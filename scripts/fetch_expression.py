#!/usr/bin/env python3
"""
Fetch the Alliance EXPRESSION-ALLIANCE combined TSV shards from FMS and
concatenate them into a single file the existing AllianceExpressionConverter
can read.

The converter expects a 23-column TSV: row 0 leads with NCBI taxon and
row 1 with species; columns 2-22 cover geneId, location, stageterm,
assayID, cellular component, anatomy, sourceUrl, source, reference, etc.

Output: data/expression-alliance.tsv
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
    configure_logging,
    fms_download_gz,
    fms_list_files,
    get_current_release,
)

log = logging.getLogger("alliance.fetch.expression")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance EXPRESSION-ALLIANCE combined TSV (FMS bulk)")
    parser.add_argument("--sub", default="COMBINED", help="dataSubType filter (default COMBINED)")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    out_path = Path(args.out_dir) / "expression-alliance.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")

    files = fms_list_files("EXPRESSION-ALLIANCE", sub_filter=args.sub)
    log.info("EXPRESSION-ALLIANCE shards: %d", len(files))
    if not files:
        log.error("No FMS files for EXPRESSION-ALLIANCE / %s", args.sub)
        sys.exit(1)

    header_written = False
    row_count = 0
    with tmp.open("w", encoding="utf-8") as out:
        for fi in files:
            local = fms_download_gz(fi["url"])
            with gzip.open(local, "rt", encoding="utf-8", errors="replace") as fh:
                seen_header = False
                for line in fh:
                    if line.startswith("#"):
                        if not header_written:
                            out.write(line)
                        continue
                    if not seen_header:
                        seen_header = True
                        if not header_written:
                            out.write(line)
                            header_written = True
                        continue
                    out.write(line)
                    row_count += 1
    tmp.replace(out_path)
    log.info("Wrote %s (%d data rows across %d shards)", out_path, row_count, len(files))


if __name__ == "__main__":
    main()
