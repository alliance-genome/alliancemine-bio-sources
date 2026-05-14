#!/usr/bin/env python3
"""
Fetch small-molecule catalog from the Alliance FMS MOLECULE datatype.
Currently only WB ships MOLECULE shards (stale at 8.3.0; included for
schema completeness ahead of any future fresh upstream).

Schema per record: id (CHEBI/MOD curie), name, iupac, smiles, inchi,
inchikey, crossReferences[].

Output: data/molecules.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import TsvWriter, configure_logging, get_current_release, iter_fms_datatype  # noqa: E402

log = logging.getLogger("alliance.fetch.molecules")


DEFAULT_MODS = ["WB", "FB", "MGI", "ZFIN", "RGD", "SGD", "HUMAN", "XBXL", "XBXT"]


COLUMNS = [
    "primaryIdentifier",
    "chebiId",
    "name",
    "iupacName",
    "smiles",
    "inchi",
    "inchiKey",
    "dataProvider",
]


def _row(rec: dict, mod: str) -> dict | None:
    pid = rec.get("id")
    if not pid:
        return None
    chebi = ""
    for x in (rec.get("crossReferences") or []):
        xid = x.get("id", "") if isinstance(x, dict) else ""
        if xid.startswith("CHEBI:"):
            chebi = xid
            break
    return {
        "primaryIdentifier": pid,
        "chebiId": chebi,
        "name": rec.get("name") or "",
        "iupacName": rec.get("iupac") or "",
        "smiles": rec.get("smiles") or "",
        "inchi": rec.get("inchi") or "",
        "inchiKey": rec.get("inchikey") or "",
        "dataProvider": mod,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance small molecules (FMS bulk)")
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
        out_dir / "molecules.tsv",
        COLUMNS,
        release=release,
        source=f"FMS MOLECULE ({','.join(mods)})",
    ) as writer:
        for mod in mods:
            log.info("Molecules FMS: %s", mod)
            for _sub, rec in iter_fms_datatype("MOLECULE", sub_filter=mod):
                row = _row(rec, mod)
                if row is None or row["primaryIdentifier"] in seen:
                    continue
                seen.add(row["primaryIdentifier"])
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    log.info("Done. rows=%d", total)
                    return
    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
