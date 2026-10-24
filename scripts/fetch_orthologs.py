#!/usr/bin/env python3
"""
Fetch the per-algorithm ortholog breakdown from the Alliance API.

The /gene/{id}/orthologs endpoint exposes data the FMS ORTHOLOGY-ALLIANCE
TSV doesn't: predictionMethodsMatched / NotMatched / NotCalled (lists of
named algorithms) and stringencyFilter. The flat aggregate counters
(AlgorithmsMatch, OutOfAlgorithms) the FMS file does have are already
populated by the alliance-orthologs converter.

This fetcher emits a companion TSV that the new alliance-ortholog-detail
bio-source converts into partial Homologue items, which InterMine merges
with the FMS-driven Homologue items via Homologue.key_pair = gene, homologue.

Usage:
    python3 fetch_orthologs.py                     # every gene of the 9 MODs
    python3 fetch_orthologs.py --mods MGI,HUMAN --workers 6
    python3 fetch_orthologs.py --ids SGD:S000004103
    python3 fetch_orthologs.py --limit 200

Output: data/orthologs.tsv
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
    resolve_gene_ids,
    run_parallel,
)

log = logging.getLogger("alliance.fetch.orthologs")


COLUMNS = [
    "subjectGeneId",
    "subjectGeneTaxon",
    "orthologGeneId",
    "orthologGeneTaxon",
    "stringencyFilter",
    "predictionMethodsMatched",
    "predictionMethodsNotMatched",
    "predictionMethodsNotCalled",
]


def _method_names(methods):
    if not methods:
        return ""
    return join_pipe(m.get("name", "") for m in methods)


def _row(seed_id: str, result: dict) -> dict:
    inner = result.get("geneToGeneOrthologyGenerated") or {}
    subj = inner.get("subjectGene") or {}
    obj = inner.get("objectGene") or {}
    if get_nested(subj, "primaryExternalId") != seed_id and get_nested(obj, "primaryExternalId") == seed_id:
        subj, obj = obj, subj
    return {
        "subjectGeneId": get_nested(subj, "primaryExternalId", default=seed_id) or seed_id,
        "subjectGeneTaxon": get_nested(subj, "taxon", "curie", default=""),
        "orthologGeneId": get_nested(obj, "primaryExternalId", default=""),
        "orthologGeneTaxon": get_nested(obj, "taxon", "curie", default=""),
        "stringencyFilter": result.get("stringencyFilter", ""),
        "predictionMethodsMatched": _method_names(inner.get("predictionMethodsMatched")),
        "predictionMethodsNotMatched": _method_names(inner.get("predictionMethodsNotMatched")),
        "predictionMethodsNotCalled": _method_names(inner.get("predictionMethodsNotCalled")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance ortholog per-algorithm detail")
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
    seen_pairs: set = set()
    total_rows = 0

    def work(gid, cache, stats):
        return [_row(gid, r) for r in paginate(f"/gene/{gid}/orthologs", cache=cache, stats=stats)]

    def consume(_gid, rows):
        nonlocal total_rows
        for row in rows:
            pair = (row["subjectGeneId"], row["orthologGeneId"])
            if not all(pair) or pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            writer.write_row(**row)
            total_rows += 1

    with TsvWriter(
        out_dir / "orthologs.tsv",
        COLUMNS,
        release=release,
        source="/gene/*/orthologs",
    ) as writer:
        report = run_parallel(gene_ids, work, consume, cache_name="orthologs", stats=stats,
                              workers=args.workers, label="genes")

    log.info("orthologs: rows=%d", total_rows)
    finish_run(report, stats, out_dir, args.max_failure_rate)

if __name__ == "__main__":
    main()
