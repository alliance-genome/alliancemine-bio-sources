#!/usr/bin/env python3
"""
Fetch RNAi screen data from WormMine via InterMine PathQuery REST.

Pulls every RNAi record and the gene it targets, emitting a TSV that
AllianceMine ingests via the alliance-worm-rnai bio-source module.
WormMine is hosted at the same Alliance server domain.

Usage:
    python3 fetch_wormmine_rnai.py
    python3 fetch_wormmine_rnai.py --limit 5000

Output: data/worm-rnai.tsv
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

log = logging.getLogger("alliance.fetch.wormmine_rnai")

WORMMINE = "https://wormmine.alliancegenome.org/wormmine"

# Pull RNAi entry plus the targeted gene plus key metadata.
# WormMine RNAi has a Gene.RNAi reverse-ref so we can join.
QUERY = """<query name="" model="genomic" view="
    RNAi.primaryIdentifier
    RNAi.method
    RNAi.treatment
    RNAi.temperature
    RNAi.genotype
    RNAi.deliveredBy
    RNAi.phenotypeRemark
    RNAi.remark
    RNAi.reference.pubMedId
    RNAi.strain.primaryIdentifier
    RNAi.inhibitsGene.primaryIdentifier
"/>"""


COLUMNS = [
    "rnaiId",
    "method",
    "treatment",
    "temperature",
    "genotype",
    "deliveredBy",
    "phenotypeRemark",
    "remark",
    "referenceId",
    "strainId",
    "geneId",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch WormMine RNAi data into TSV")
    parser.add_argument("--limit", type=int, help="Process at most N rows")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    log.info("Fetching RNAi data from WormMine...")
    release = get_current_release()
    out_dir = Path(args.out_dir)
    stats = FetchStats()

    with open_cache("wormmine_rnai") as cache, \
         TsvWriter(
             out_dir / "worm-rnai.tsv",
             COLUMNS,
             release=release,
             source="WormMine PathQuery /service/query/results",
         ) as writer:
        n = 0
        for row in intermine_paginate(WORMMINE, QUERY, page_size=2000, cache=cache, stats=stats):
            writer.write_row(**{
                col: (row[i] if i < len(row) and row[i] is not None else "")
                for i, col in enumerate(COLUMNS)
            })
            n += 1
            if args.limit and n >= args.limit:
                break

    log.info("Done. rows=%d  API stats: %s", n, stats)


if __name__ == "__main__":
    main()
