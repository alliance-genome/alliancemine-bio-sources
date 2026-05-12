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
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    enumerate_mod_genes,
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    iter_fms_datatype,
    join_pipe,
    normalize_taxon,
    open_cache,
    paginate,
)

log = logging.getLogger("alliance.fetch.transgenic_alleles")


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
        default="MGI,ZFIN,FB,WB,RGD,XBXL,XBXT,HUMAN",
        help="Comma-separated MOD prefixes for seed gene enumeration (default: non-yeast)",
    )
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=20,
        help="ThreadPool size for parallel per-gene API calls (default 20)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    out_dir = Path(args.out_dir)
    stats = FetchStats()
    seen_alleles: set = set()

    # FMS path: bulk JSON.gz per MOD. No per-gene fan-out, no WAF.
    if args.source == "fms":
        mods = [m.strip() for m in args.mods.split(",") if m.strip()]
        with TsvWriter(
            out_dir / "transgenic-alleles.tsv",
            COLUMNS,
            release=release,
            source=f"FMS CONSTRUCT+ALLELE ({','.join(mods)})",
        ) as writer:
            total = _run_fms(mods, writer, seen_alleles, args.limit)
        log.info("FMS Done. rows=%d", total)
        return

    if args.ids:
        gene_ids = [g.strip() for g in args.ids.split(",") if g.strip()]
    else:
        mods = [m.strip() for m in args.mods.split(",") if m.strip()]
        if mods == ["SGD"]:
            gene_ids = enumerate_yeast_genes()
        else:
            gene_ids = enumerate_mod_genes(mods)
    if args.limit:
        gene_ids = gene_ids[: args.limit]

    log.info("Processing %d genes", len(gene_ids))

    def fetch_one(gid: str) -> tuple[str, list[dict], Exception | None]:
        rows: list[dict] = []
        try:
            # Per-worker sqlite connection — sqlite3.Connection is not safe
            # to share across threads even in WAL mode.
            with open_cache("transgenic_alleles") as cache:
                for result in paginate(
                    f"/gene/{gid}/transgenic-alleles", cache=cache, stats=stats
                ):
                    rows.append(_row(gid, result))
            return gid, rows, None
        except Exception as e:
            return gid, [], e

    with TsvWriter(
            out_dir / "transgenic-alleles.tsv",
            COLUMNS,
            release=release,
            source="/gene/*/transgenic-alleles",
         ) as writer:
        total_rows = 0
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(fetch_one, gid): gid for gid in gene_ids}
            try:
                from tqdm import tqdm
                pbar = tqdm(total=len(futures), desc="genes", unit="gene")
            except ImportError:
                pbar = None
            for fut in as_completed(futures):
                gid, rows, err = fut.result()
                if pbar:
                    pbar.update(1)
                if err is not None:
                    log.warning("Transgenic-allele fetch failed for %s: %s", gid, err)
                    continue
                for row in rows:
                    key = row["transgenicAlleleId"]
                    if not key or key in seen_alleles:
                        continue
                    seen_alleles.add(key)
                    writer.write_row(**row)
                    total_rows += 1
            if pbar:
                pbar.close()

    log.info("Done. rows=%d  API stats: %s", total_rows, stats)


if __name__ == "__main__":
    main()
