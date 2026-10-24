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
    python3 fetch_phenotypes.py                    # every gene of the 9 MODs
    python3 fetch_phenotypes.py --mods MGI --workers 6
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
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    get_nested,
    join_pipe,
    paginate,
    publication_ids,
    resolve_gene_ids,
    run_parallel,
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
    pubmed_ids = publication_ids(result)
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
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--mods", help="Comma-separated MODs to seed from (default: all 9, env ALLIANCE_MODS)")
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
    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit)
    out_dir = Path(args.out_dir)
    stats = FetchStats()
    seen_uids: set = set()
    total_rows = 0

    def work(gid, cache, stats):
        return [_row(gid, r) for r in paginate(f"/gene/{gid}/phenotypes", cache=cache, stats=stats)]

    def consume(_gid, rows):
        nonlocal total_rows
        for row in rows:
            uid = row["uniqueId"]
            if uid and uid in seen_uids:
                continue
            if uid:
                seen_uids.add(uid)
            writer.write_row(**row)
            total_rows += 1

    with TsvWriter(
        out_dir / "phenotypes.tsv",
        COLUMNS,
        release=release,
        source="/gene/*/phenotypes",
    ) as writer:
        report = run_parallel(gene_ids, work, consume, cache_name="phenotypes", stats=stats,
                              workers=args.workers, label="genes")

    log.info("phenotypes: rows=%d", total_rows)
    finish_run(report, stats, out_dir, args.max_failure_rate)

if __name__ == "__main__":
    main()
