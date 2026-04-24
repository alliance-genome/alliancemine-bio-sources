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
    python3 fetch_orthologs.py
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
    configure_logging,
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    join_pipe,
    open_cache,
    paginate,
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
    seen_pairs: set = set()

    with open_cache("orthologs") as cache, \
         TsvWriter(
             out_dir / "orthologs.tsv",
             COLUMNS,
             release=release,
             source="/gene/*/orthologs",
         ) as writer:
        total_rows = 0
        for gid in iterator:
            try:
                for result in paginate(f"/gene/{gid}/orthologs", cache=cache, stats=stats):
                    row = _row(gid, result)
                    pair = (row["subjectGeneId"], row["orthologGeneId"])
                    if not all(pair) or pair in seen_pairs:
                        continue
                    seen_pairs.add(pair)
                    writer.write_row(**row)
                    total_rows += 1
            except Exception as e:
                log.warning("Ortholog fetch failed for %s: %s", gid, e)

    log.info("Done. rows=%d  API stats: %s", total_rows, stats)


if __name__ == "__main__":
    main()
