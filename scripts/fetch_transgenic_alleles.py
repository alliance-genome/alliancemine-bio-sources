#!/usr/bin/env python3
"""
Fetch transgenic alleles from the Alliance API.

For each seed gene, call /gene/{id}/transgenic-alleles paginated and emit one
row per transgenic-allele record. MGI / ZFIN curate most transgenic alleles;
yeast typically returns 0.

Usage:
    python3 fetch_transgenic_alleles.py
    python3 fetch_transgenic_alleles.py --mods MGI,ZFIN
    python3 fetch_transgenic_alleles.py --ids MGI:88276
    python3 fetch_transgenic_alleles.py --limit 100

Output: data/transgenic-alleles.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    get_nested,
    iter_fms_datatype,
    join_pipe,
    normalize_taxon,
    paginate,
    resolve_gene_ids,
    run_parallel,
    selected_mods,
)

log = logging.getLogger("alliance.fetch.transgenic_alleles")

DEFAULT_MODS = ["MGI", "ZFIN", "FB", "WB", "RGD", "XBXL", "XBXT", "HUMAN"]


COLUMNS = [
    "geneId",
    "geneSymbol",
    "geneTaxon",
    "transgenicAlleleId",
    "transgenicAlleleSymbol",
    "constructIds",
    "dataProvider",
]


def _row(seed_id: str, result: dict) -> dict:
    gene = result.get("gene") or {}
    allele_doc = result.get("alleleDocument") or {}
    allele = allele_doc.get("allele") or {}
    constructs = allele_doc.get("transgenicAlleleConstructs") or []
    construct_ids = join_pipe(
        get_nested(c, "construct", "primaryExternalId", default="")
        for c in constructs
        if get_nested(c, "construct", "primaryExternalId")
    )
    return {
        "geneId": get_nested(gene, "primaryExternalId", default=seed_id),
        "geneSymbol": get_nested(gene, "geneSymbol", "displayText", default=""),
        "geneTaxon": normalize_taxon(get_nested(gene, "taxon", "curie", default="")),
        "transgenicAlleleId": get_nested(allele, "primaryExternalId", default=""),
        "transgenicAlleleSymbol": get_nested(allele, "alleleSymbol", "displayText", default=""),
        "constructIds": construct_ids,
        "dataProvider": get_nested(result, "dataProvider", "abbreviation", default=""),
    }


def _load_construct_index(mod: str) -> dict[str, dict]:
    """Load CONSTRUCT_{MOD} into a dict keyed by primaryId. The CONSTRUCT
    dumps are small (~13k records per MOD) so an in-memory dict is fine."""
    idx: dict[str, dict] = {}
    for _sub, rec in iter_fms_datatype("CONSTRUCT", sub_filter=mod):
        pid = rec.get("primaryId")
        if pid:
            idx[pid] = rec
    return idx


def _fms_rows_for_mod(mod: str) -> Iterator[dict]:
    """Yield TSV-row dicts for every transgenic allele in the given MOD's
    ALLELE FMS shard. Joins to construct components on demand."""
    constructs = _load_construct_index(mod)
    log.info("FMS %s: %d constructs indexed", mod, len(constructs))
    for _sub, allele in iter_fms_datatype("ALLELE", sub_filter=mod):
        relations = allele.get("alleleObjectRelations") or []
        construct_ids: list[str] = []
        gene_id = ""
        for rel in relations:
            rel_obj = rel.get("objectRelation") or {}
            if rel_obj.get("associationType") == "contains" and rel_obj.get("construct"):
                construct_ids.append(rel_obj["construct"])
            elif rel_obj.get("associationType") == "allele_of" and rel_obj.get("gene") and not gene_id:
                gene_id = rel_obj["gene"]
        if not construct_ids:
            continue
        # Aggregate component symbols from joined construct records.
        component_syms: list[str] = []
        for cid in construct_ids:
            comp = constructs.get(cid)
            if not comp:
                continue
            for c in (comp.get("constructComponents") or []):
                sym = c.get("componentSymbol")
                if sym:
                    component_syms.append(sym)
        yield {
            "geneId": gene_id,
            "geneSymbol": "",  # not in ALLELE shard; AGR API path fills this
            "geneTaxon": normalize_taxon(allele.get("taxonId", "")),
            "transgenicAlleleId": allele.get("primaryId", ""),
            "transgenicAlleleSymbol": allele.get("symbolText") or allele.get("symbol", ""),
            "constructIds": join_pipe(construct_ids),
            "dataProvider": mod,
        }


def _run_fms(mods: list[str], writer: TsvWriter, seen: set, limit: int | None) -> int:
    total = 0
    for mod in mods:
        log.info("FMS extraction: %s", mod)
        for row in _fms_rows_for_mod(mod):
            key = row["transgenicAlleleId"]
            if not key or key in seen:
                continue
            seen.add(key)
            writer.write_row(**row)
            total += 1
            if limit is not None and total >= limit:
                return total
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance transgenic alleles into TSV")
    parser.add_argument(
        "--source",
        choices=["fms", "api"],
        default="fms",
        help="Backend: 'fms' (bulk JSON.gz from download.alliancegenome.org; default) or "
             "'api' (per-gene /gene/{id}/transgenic-alleles; needed for fields not in FMS).",
    )
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (only with --source api)")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    add_parallel_args(parser)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    out_dir = Path(args.out_dir)
    stats = FetchStats()
    seen_alleles: set = set()

    # FMS path: bulk JSON.gz per MOD. No per-gene fan-out, no WAF.
    if args.source == "fms":
        mods = selected_mods(DEFAULT_MODS, args.mods)
        with TsvWriter(
            out_dir / "transgenic-alleles.tsv",
            COLUMNS,
            release=release,
            source=f"FMS CONSTRUCT+ALLELE ({','.join(mods)})",
        ) as writer:
            total = _run_fms(mods, writer, seen_alleles, args.limit)
        log.info("FMS Done. rows=%d", total)
        return

    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit, default_mods=DEFAULT_MODS)
    total_rows = 0

    def work(gid, cache, stats):
        return [_row(gid, r) for r in paginate(f"/gene/{gid}/transgenic-alleles", cache=cache, stats=stats)]

    def consume(_gid, rows):
        nonlocal total_rows
        for row in rows:
            key = row["transgenicAlleleId"]
            if not key or key in seen_alleles:
                continue
            seen_alleles.add(key)
            writer.write_row(**row)
            total_rows += 1

    with TsvWriter(
        out_dir / "transgenic-alleles.tsv",
        COLUMNS,
        release=release,
        source="/gene/*/transgenic-alleles",
    ) as writer:
        report = run_parallel(gene_ids, work, consume, cache_name="transgenic_alleles",
                              stats=stats, workers=args.workers, label="genes")

    log.info("transgenic_alleles: rows=%d", total_rows)
    finish_run(report, stats, out_dir, args.max_failure_rate)

if __name__ == "__main__":
    main()
