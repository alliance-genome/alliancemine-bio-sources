#!/usr/bin/env python3
"""
Fetch Affected Genomic Models (AGMs) from the Alliance FMS AGM bulk files.
Maps to DiseaseModel items in the genomic model. An AGM is a strain,
genotype, or fish-style model that disease/phenotype annotations reference
as their experimental subject.

Usage:
    python3 fetch_agms.py
    python3 fetch_agms.py --mods MGI,RGD --limit 50

Output: data/agms.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    get_current_release,
    iter_fms_datatype,
    normalize_taxon,
)

log = logging.getLogger("alliance.fetch.agms")


DEFAULT_MODS = ["MGI", "ZFIN", "FB", "WB", "RGD", "SGD", "XBXL", "XBXT"]


COLUMNS = [
    "modelId",
    "modelName",
    "modelSubtype",
    "taxon",
    "dataProvider",
]


def _row(rec: dict, mod: str) -> dict | None:
    model_id = rec.get("primaryID") or rec.get("primaryId")
    if not model_id:
        return None
    return {
        "modelId": model_id,
        "modelName": rec.get("name") or "",
        "modelSubtype": rec.get("subtype") or "",
        "taxon": normalize_taxon(rec.get("taxonId") or ""),
        "dataProvider": mod,
    }


def _iter_mod(mod: str) -> Iterator[dict]:
    for _sub, rec in iter_fms_datatype("AGM", sub_filter=mod):
        row = _row(rec, mod)
        if row is not None:
            yield row


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance AGMs (FMS bulk)")
    parser.add_argument(
        "--mods",
        default=",".join(DEFAULT_MODS),
        help="Comma-separated MOD prefixes",
    )
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    mods = [m.strip() for m in args.mods.split(",") if m.strip()]
    release = get_current_release()
    out_dir = Path(args.out_dir)
    seen: set = set()
    total = 0

    with TsvWriter(
        out_dir / "agms.tsv",
        COLUMNS,
        release=release,
        source=f"FMS AGM ({','.join(mods)})",
    ) as writer:
        for mod in mods:
            log.info("AGMs FMS: %s", mod)
            for row in _iter_mod(mod):
                if row["modelId"] in seen:
                    continue
                seen.add(row["modelId"])
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    log.info("Done. rows=%d", total)
                    return

    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
