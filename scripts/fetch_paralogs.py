#!/usr/bin/env python3
"""
Fetch paralog pairs from the Alliance API.

For each seed gene, call /gene/{id}/paralogs paginated and emit one row per
pair with the prediction-methods summary and the pairwise identity/similarity
metrics.

Usage:
    python3 fetch_paralogs.py                      # every gene of the 9 MODs
    python3 fetch_paralogs.py --mods FB,WB
    python3 fetch_paralogs.py --ids SGD:S000004103
    python3 fetch_paralogs.py --limit 200

Output: data/paralogs.tsv
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

log = logging.getLogger("alliance.fetch.paralogs")


COLUMNS = [
    "subjectGeneId",
    "subjectGeneSymbol",
    "subjectGeneTaxon",
    "paralogueGeneId",
    "paralogueGeneSymbol",
    "paralogueGeneTaxon",
    "identity",
    "similarity",
    "length",
    "rank",
    "predictionMethodsMatched",
    "predictionMethodsNotMatched",
    "predictionMethodsNotCalled",
]


def _method_names(methods):
    if not methods:
        return ""
    return join_pipe(m.get("name", "") for m in methods)


def _gene_fields(node: dict, prefix: str) -> dict:
    return {
        f"{prefix}Id": get_nested(node, "primaryExternalId", default=""),
        f"{prefix}Symbol": get_nested(node, "geneSymbol", "displayText", default=""),
        f"{prefix}Taxon": get_nested(node, "taxon", "curie", default=""),
    }


def _row(seed_id: str, result: dict) -> dict:
    p = result.get("geneToGeneParalogy") or {}
    subj = p.get("subjectGene") or {}
    obj = p.get("objectGene") or {}
    # Normalise so the seed gene is always the "subject" side.
    if subj.get("primaryExternalId") != seed_id and obj.get("primaryExternalId") == seed_id:
        subj, obj = obj, subj
    return {
        **_gene_fields(subj, "subjectGene"),
        **_gene_fields(obj, "paralogueGene"),
        "identity": p.get("identity", "") or "",
        "similarity": p.get("similarity", "") or "",
        "length": p.get("length", "") or "",
        "rank": p.get("rank", "") or "",
        "predictionMethodsMatched": _method_names(p.get("predictionMethodsMatched")),
        "predictionMethodsNotMatched": _method_names(p.get("predictionMethodsNotMatched")),
        "predictionMethodsNotCalled": _method_names(p.get("predictionMethodsNotCalled")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance paralogs into TSV")
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
    # Dedupe by unordered gene-pair so we don't emit each pair twice from
    # both ends of the traversal.
    seen_pairs: set = set()
    total_rows = 0

    def work(gid, cache, stats):
        return [_row(gid, r) for r in paginate(f"/gene/{gid}/paralogs", cache=cache, stats=stats)]

    def consume(_gid, rows):
        nonlocal total_rows
        for row in rows:
            pair_key = tuple(sorted([row["subjectGeneId"], row["paralogueGeneId"]]))
            if not all(pair_key) or pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            writer.write_row(**row)
            total_rows += 1

    with TsvWriter(
        out_dir / "paralogs.tsv",
        COLUMNS,
        release=release,
        source="/gene/*/paralogs",
    ) as writer:
        report = run_parallel(gene_ids, work, consume, cache_name="paralogs", stats=stats,
                              workers=args.workers, label="genes")
        finish_run(report, stats, out_dir, args.max_failure_rate,
                   rows_written=writer.rows_written)

    log.info("paralogs: rows=%d", total_rows)

if __name__ == "__main__":
    main()
