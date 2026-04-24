#!/usr/bin/env python3
"""
Fetch the per-disease gene-association detail from the Alliance API.

The /disease/{id}/genes endpoint exposes fields the FMS DISEASE-ALLIANCE TSV
doesn't surface: generatedRelationString (e.g. "is_implicated_in"),
diseaseQualifiers (e.g. ["onset_of"]), evidenceCodes (ECO term references),
parentSlimIds (DOID slim terms), and viaOrthologyOrder (0 when direct, >0
when inherited via ortholog chain).

Seed list is the distinct DOID column from the existing FMS DISEASE TSV.
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
import gzip
import logging
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    get_current_release,
    get_nested,
    join_pipe,
    open_cache,
    paginate,
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
    pubmed_ids = result.get("pubmedPubModIDs") or []
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


def _enumerate_diseases_from_combined() -> list[str]:
    url = "https://fms.alliancegenome.org/download/DISEASE-ALLIANCE_COMBINED.tsv.gz"
    raw = urllib.request.urlopen(url, timeout=120).read()
    text = gzip.decompress(raw).decode("utf-8", errors="replace")
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    if not lines:
        return []
    header = lines[0].split("\t")
    try:
        idx = header.index("DOID")
    except ValueError:
        return []
    seen: set = set()
    out: list[str] = []
    for row in lines[1:]:
        cols = row.split("\t")
        if len(cols) > idx:
            d = cols[idx].strip()
            if d and d not in seen:
                seen.add(d)
                out.append(d)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance disease annotation detail")
    parser.add_argument("--ids", help="Comma-separated DOID CURIEs")
    parser.add_argument("--limit", type=int, help="Process at most N diseases")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    if args.ids:
        ids = [d.strip() for d in args.ids.split(",") if d.strip()]
    else:
        log.info("Enumerating disease IDs from FMS DISEASE-ALLIANCE_COMBINED...")
        ids = _enumerate_diseases_from_combined()
    if args.limit:
        ids = ids[: args.limit]

    log.info("Processing %d diseases", len(ids))
    release = get_current_release()
    out_dir = Path(args.out_dir)

    try:
        from tqdm import tqdm
        iterator = tqdm(ids, desc="diseases", unit="disease")
    except ImportError:
        iterator = ids

    stats = FetchStats()
    seen_uids: set = set()

    with open_cache("disease_annotations") as cache, \
         TsvWriter(
             out_dir / "disease-annotations-detail.tsv",
             COLUMNS,
             release=release,
             source="/disease/*/genes",
         ) as writer:
        total_rows = 0
        for did in iterator:
            try:
                for result in paginate(f"/disease/{did}/genes", cache=cache, stats=stats):
                    row = _row(did, result)
                    uid = row["uniqueId"]
                    if uid and uid in seen_uids:
                        continue
                    if uid:
                        seen_uids.add(uid)
                    writer.write_row(**row)
                    total_rows += 1
            except Exception as e:
                log.warning("Disease-annotation fetch failed for %s: %s", did, e)

    log.info("Done. rows=%d  API stats: %s", total_rows, stats)


if __name__ == "__main__":
    main()
