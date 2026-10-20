#!/usr/bin/env python3
"""
Fetch the per-disease gene-association detail from the Alliance API.

The /disease/{id}/genes endpoint exposes fields the FMS DISEASE-ALLIANCE TSV
doesn't surface: generatedRelationString (e.g. "is_implicated_in"),
diseaseQualifiers (e.g. ["onset_of"]), evidenceCodes (ECO term references),
parentSlimIds (DOID slim terms), and viaOrthologyOrder (0 when direct, >0
when inherited via ortholog chain).

Seed list is the distinct DOID column of the pinned release's FMS
DISEASE-ALLIANCE COMBINED TSV.
For each disease we paginate the genes endpoint and emit one row per
(disease, gene) annotation.

Usage:
    python3 fetch_disease_annotations.py
    python3 fetch_disease_annotations.py --ids DOID:14330
    python3 fetch_disease_annotations.py --from-tsv path/to/DISEASE-ALLIANCE.tsv
    python3 fetch_disease_annotations.py --limit 100

Output: data/disease-annotations-detail.tsv
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
    enumerate_disease_ids,
    finish_run,
    get_current_release,
    get_nested,
    join_pipe,
    paginate,
    publication_ids,
    run_parallel,
)

log = logging.getLogger("alliance.fetch.disease_annotations")


COLUMNS = [
    "diseaseId",
    "diseaseName",
    "subjectGeneId",
    "subjectGeneTaxon",
    "generatedRelationString",
    "diseaseQualifiers",
    "evidenceCodes",
    "parentSlimIds",
    "viaOrthologyOrder",
    "pubmedIds",
    "referenceIds",
    "uniqueId",
]


def _row(disease_id: str, result: dict) -> dict:
    subject = result.get("subject") or {}
    obj = result.get("object") or {}
    refs = result.get("references") or []
    pubmed_ids = publication_ids(result)
    ec = result.get("evidenceCodes") or []
    qualifiers = result.get("diseaseQualifiers") or []
    parent_slim = result.get("parentSlimIDs") or []
    via_order = result.get("viaOrthologyOrder")
    return {
        "diseaseId": obj.get("curie", "") or disease_id,
        "diseaseName": obj.get("name", ""),
        "subjectGeneId": get_nested(subject, "primaryExternalId", default=""),
        "subjectGeneTaxon": get_nested(subject, "taxon", "curie", default=""),
        "generatedRelationString": result.get("generatedRelationString", ""),
        "diseaseQualifiers": join_pipe(qualifiers),
        "evidenceCodes": join_pipe(c.get("curie", "") for c in ec),
        "parentSlimIds": join_pipe(parent_slim),
        "viaOrthologyOrder": str(via_order) if via_order is not None else "",
        "pubmedIds": join_pipe(pubmed_ids),
        "referenceIds": join_pipe(r.get("curie", "") for r in refs),
        "uniqueId": result.get("uniqueId", ""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance disease annotation detail")
    parser.add_argument("--ids", help="Comma-separated DOID CURIEs")
    parser.add_argument("--limit", type=int, help="Process at most N diseases")
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
    if args.ids:
        ids = [d.strip() for d in args.ids.split(",") if d.strip()]
    else:
        log.info("Enumerating disease IDs from FMS DISEASE-ALLIANCE COMBINED (release %s)...", release)
        ids = enumerate_disease_ids()
    if args.limit:
        ids = ids[: args.limit]

    out_dir = Path(args.out_dir)
    stats = FetchStats()
    seen_uids: set = set()
    total_rows = 0

    def work(did, cache, stats):
        return [_row(did, r) for r in paginate(f"/disease/{did}/genes", cache=cache, stats=stats)]

    def consume(_did, rows):
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
        out_dir / "disease-annotations-detail.tsv",
        COLUMNS,
        release=release,
        source="/disease/*/genes",
    ) as writer:
        report = run_parallel(ids, work, consume, cache_name="disease_annotations", stats=stats,
                              workers=args.workers, label="diseases")
        finish_run(report, stats, out_dir, args.max_failure_rate,
                   rows_written=writer.rows_written)

    log.info("disease_annotations: rows=%d", total_rows)


if __name__ == "__main__":
    main()
