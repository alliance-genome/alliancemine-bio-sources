#!/usr/bin/env python3
"""
Build allele-detail.tsv (see fetch_allele_detail.py) from a restored Alliance
curation database instead of one /allele/{id} call per allele.

/api/allele/{id} serves the allele summary document. Its three captured
fields come from the curation DB as follows (AlleleSummaryDocumentBuilder):
    alterationType     "allele" / "allele with one variant" / "allele with
                       multiple variants" by the count of live variant links
    category           always "allele_summary"
    crossReference     the data provider cross reference's displayName

Scope: the AlleleIds of the release's FMS VARIANT-ALLELE files, as the fetcher
seeds them, limited to alleles the API serves (not obsolete, not internal).
Seed ids absent from the DB or not served are logged and skipped; the API
answers HTTP 400 for them. --all-alleles takes every summary allele instead.

Usage:
    python3 extract_allele_detail.py                  # alleles of all MODs
    python3 extract_allele_detail.py --mods SGD,WB
    python3 extract_allele_detail.py --ids SGD:S000279149
    python3 extract_allele_detail.py --all-alleles --limit 100

Output: data/allele-detail.tsv
"""

from __future__ import annotations

import argparse
import collections
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import TsvWriter, configure_logging, get_current_release, selected_mods  # noqa: E402
from curation_db import copy_rows  # noqa: E402
from extract_alleles_common import (  # noqa: E402
    VARIANT_ALLELE_MODS,
    allele_variant_links,
    alteration_type,
    mod_for,
    seed_allele_ids,
    summary_alleles,
)
from fetch_allele_detail import COLUMNS  # noqa: E402

log = logging.getLogger("alliance.extract.allele_detail")

API_CATEGORY = "allele_summary"


def detail_row(pid: str, variant_count: int, xref_display: str) -> dict:
    return {
        "alleleId": pid,
        "alterationType": alteration_type(variant_count),
        "apiCategory": API_CATEGORY,
        "apiCrossReference": xref_display,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract per-allele detail from the curation DB")
    parser.add_argument("--ids", help="Comma-separated allele CURIEs (skips the FMS seed)")
    parser.add_argument("--mods", help="Comma-separated MODs whose alleles to emit (default: all, env ALLIANCE_MODS)")
    parser.add_argument("--all-alleles", action="store_true",
                        help="Every allele the API serves, not only the FMS VARIANT-ALLELE seed")
    parser.add_argument("--limit", type=int, help="Emit at most N alleles")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    release = get_current_release()
    mods = selected_mods(None, args.mods)
    if args.ids:
        wanted = {a.strip() for a in args.ids.split(",") if a.strip()}
    else:
        wanted = seed_allele_ids([m for m in mods if m in VARIANT_ALLELE_MODS] or mods, args.all_alleles)

    links = allele_variant_links(copy_rows)
    per_mod: collections.Counter = collections.Counter()
    seen: set[str] = set()
    with TsvWriter(
        Path(args.out_dir) / "allele-detail.tsv",
        COLUMNS,
        release=release,
        source="curation DB allele summary (as /allele/{id})",
    ) as writer:
        for db_id, pid, provider, taxon, _symbol, _fmt, xref_display in summary_alleles(copy_rows, mods):
            if pid in seen or (wanted is not None and pid not in wanted):
                continue
            seen.add(pid)
            writer.write_row(**detail_row(pid, len(links.get(db_id, [])), xref_display))
            per_mod[mod_for(provider, taxon)] += 1
            if args.limit and writer.rows_written >= args.limit:
                break
    if wanted is not None and not args.limit:
        missing = wanted - seen
        log.info("%d seed alleles not served (absent, obsolete or internal in the DB)", len(missing))
    log.info("allele_detail rows per MOD: %s", dict(sorted(per_mod.items())))


if __name__ == "__main__":
    main()
