#!/usr/bin/env python3
"""
Fetch gene-level phenotype annotations from the Alliance API.

The API phenotype endpoint returns deeply nested per-annotation objects with
alleles, conditions, ontology terms, and references. For this first pass we
capture the top-level gene-phenotype statement + relation + citation list,
which is enough to populate a queryable PhenotypeAnnotation entity per gene.
Allele-level detail and condition-relation modelling is punted to a later
phase.

Usage:
    python3 fetch_phenotypes.py
    python3 fetch_phenotypes.py --ids SGD:S000004103

Output: data/phenotypes.tsv
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
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    join_pipe,
    open_cache,
    paginate,
)

log = logging.getLogger("alliance.fetch.phenotypes")


COLUMNS = [
    "subjectGeneId",
    "subjectGeneTaxon",
    "category",
    "phenotypeStatement",
    "relation",
    "uniqueId",
    "pubmedIds",
    "referenceIds",
]


def _row(seed_id: str, result: dict) -> dict:
    subject = result.get("subject") or {}
    refs = result.get("references") or []
    pubmed_ids = result.get("pubmedPubModIDs") or []
    return {
        "subjectGeneId": get_nested(subject, "primaryExternalId", default=seed_id) or seed_id,
        "subjectGeneTaxon": get_nested(subject, "taxon", "curie", default=""),
        "category": result.get("category", ""),
        "phenotypeStatement": result.get("phenotypeStatement", ""),
        "relation": get_nested(result, "relation", "name", default=""),
        "uniqueId": result.get("uniqueId", ""),
        "pubmedIds": join_pipe(pubmed_ids),
        "referenceIds": join_pipe(r.get("curie", "") for r in refs),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance phenotype annotations")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs")
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
        gene_ids = enumerate_yeast_genes()
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
    seen_uids: set = set()

    with open_cache("phenotypes") as cache, \
         TsvWriter(
             out_dir / "phenotypes.tsv",
             COLUMNS,
             release=release,
             source="/gene/*/phenotypes",
         ) as writer:
        total_rows = 0
        for gid in iterator:
            try:
                for result in paginate(f"/gene/{gid}/phenotypes", cache=cache, stats=stats):
                    row = _row(gid, result)
                    uid = row["uniqueId"]
                    if uid and uid in seen_uids:
                        continue
                    if uid:
                        seen_uids.add(uid)
                    writer.write_row(**row)
                    total_rows += 1
            except Exception as e:
                log.warning("Phenotype fetch failed for %s: %s", gid, e)

    log.info("Done. rows=%d  API stats: %s", total_rows, stats)


if __name__ == "__main__":
    main()
