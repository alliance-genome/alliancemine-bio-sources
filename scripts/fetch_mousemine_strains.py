#!/usr/bin/env python3
"""
Fetch mouse strain data from MouseMine via InterMine PathQuery REST.

The Alliance API doesn't expose Strain (it's MOD-curated metadata that
MGI maintains in MouseMine). Pulls every Strain plus the alleles each
strain carries, emits a TSV that AllianceMine ingests via the
alliance-mouse-strains bio-source module.

Usage:
    python3 fetch_mousemine_strains.py
    python3 fetch_mousemine_strains.py --limit 1000

Output: data/mouse-strains.tsv
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
    get_current_release,
    intermine_paginate,
    open_cache,
)

log = logging.getLogger("alliance.fetch.mousemine_strains")

MOUSEMINE = "https://www.mousemine.org/mousemine"

# Two-pass query: pull strain metadata + each strain's allele list.
# InterMine returns nulls when a strain has no alleles, so we get a
# row per (strain, allele) plus rows for strain-with-no-allele.
QUERY = """<query name="" model="genomic" view="
    Strain.primaryIdentifier
    Strain.name
    Strain.attributeString
    Strain.carries.primaryIdentifier
"/>"""


COLUMNS = [
    "strainId",
    "strainName",
    "attributeString",
    "alleleId",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch MouseMine strain data into TSV")
    parser.add_argument("--limit", type=int, help="Process at most N rows")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    log.info("Fetching mouse strains from MouseMine...")
    release = get_current_release()
    out_dir = Path(args.out_dir)
    stats = FetchStats()

    with open_cache("mousemine_strains") as cache, \
         TsvWriter(
             out_dir / "mouse-strains.tsv",
             COLUMNS,
             release=release,
             source="MouseMine PathQuery /service/query/results",
         ) as writer:
        n = 0
        for row in intermine_paginate(MOUSEMINE, QUERY, page_size=2000, cache=cache, stats=stats):
            # row is positional: [strainId, strainName, attributeString, alleleId]
            writer.write_row(
                strainId=row[0] or "",
                strainName=row[1] or "",
                attributeString=row[2] or "",
                alleleId=row[3] or "",
            )
            n += 1
            if args.limit and n >= args.limit:
                break

    log.info("Done. rows=%d  API stats: %s", n, stats)


if __name__ == "__main__":
    main()
