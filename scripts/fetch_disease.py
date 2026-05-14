#!/usr/bin/env python3
"""
Fetch the Alliance DISEASE-ALLIANCE combined TSV shards from FMS and
concatenate them into a single file the existing AllianceDiseaseConverter
can read.

The converter expects an 18-column TSV: taxon, species, productType,
productId, productSymbol, qualifier, DOid, evidenceCode, withText,
inferredFrom, ..., publication, dateAssigned, dataSource. FMS publishes
one shard per ~500k rows; this fetcher pulls all 14 9.0.0 shards
sequentially and emits a single concatenated TSV, dropping the
per-shard comment headers but preserving the column header line.

Output: data/disease-alliance.tsv
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
    FMS_BASE,
    configure_logging,
    fms_download_gz,
    fms_list_files,
    get_current_release,
)

log = logging.getLogger("alliance.fetch.disease")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance DISEASE-ALLIANCE combined TSV (FMS bulk)")
    parser.add_argument("--sub", default="COMBINED", help="dataSubType filter (default COMBINED)")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    out_path = Path(args.out_dir) / "disease-alliance.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")

    files = fms_list_files("DISEASE-ALLIANCE", sub_filter=args.sub)
    log.info("DISEASE-ALLIANCE shards: %d", len(files))
    if not files:
        log.error("No FMS files for DISEASE-ALLIANCE / %s", args.sub)
        sys.exit(1)

    header_written = False
    row_count = 0
    with tmp.open("w", encoding="utf-8") as out:
        for fi in files:
            local = fms_download_gz(fi["url"])
            with gzip.open(local, "rt", encoding="utf-8", errors="replace") as fh:
                # Each shard begins with several "#" comment rows then a
                # column header line. We keep the header from the first
                # shard only; subsequent shards' headers are skipped.
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
