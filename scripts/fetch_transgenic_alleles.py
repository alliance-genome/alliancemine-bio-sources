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

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    enumerate_mod_genes,
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    join_pipe,
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
        "geneTaxon": get_nested(gene, "taxon", "curie", default=""),
        "transgenicAlleleId": get_nested(allele, "primaryExternalId", default=""),
        "transgenicAlleleSymbol": get_nested(allele, "alleleSymbol", "displayText", default=""),
        "constructIds": construct_ids,
        "dataProvider": get_nested(result, "dataProvider", "abbreviation", default=""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance transgenic alleles into TSV")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs")
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
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

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
    release = get_current_release()
    out_dir = Path(args.out_dir)

    try:
        from tqdm import tqdm
        iterator = tqdm(gene_ids, desc="genes", unit="gene")
    except ImportError:
        iterator = gene_ids

    stats = FetchStats()
    seen_alleles: set = set()

    with open_cache("transgenic_alleles") as cache, \
         TsvWriter(
             out_dir / "transgenic-alleles.tsv",
             COLUMNS,
             release=release,
             source="/gene/*/transgenic-alleles",
         ) as writer:
        total_rows = 0
        for gid in iterator:
            try:
                for result in paginate(
                    f"/gene/{gid}/transgenic-alleles", cache=cache, stats=stats
                ):
                    row = _row(gid, result)
                    key = row["transgenicAlleleId"]
                    if not key or key in seen_alleles:
                        continue
                    seen_alleles.add(key)
                    writer.write_row(**row)
                    total_rows += 1
            except Exception as e:
                log.warning("Transgenic-allele fetch failed for %s: %s", gid, e)

    log.info("Done. rows=%d  API stats: %s", total_rows, stats)


if __name__ == "__main__":
    main()
