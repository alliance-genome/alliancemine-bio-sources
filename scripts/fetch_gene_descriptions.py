#!/usr/bin/env python3
"""
Fetch curated gene descriptions from the Alliance FMS GENE-DESCRIPTION-JSON
bulk files. One record per gene with auto-generated descriptions for the
overall summary plus GO function/process/component, disease ontology,
expression, and orthology slices.

The fields populate enrichment columns on the existing Gene class
(autoDescription, goDescription, goFunctionDescription, etc.) added in
alliancemine-global_additions.xml.

Usage:
    python3 fetch_gene_descriptions.py
    python3 fetch_gene_descriptions.py --mods MGI,SGD --limit 100

Output: data/gene-descriptions.tsv
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
)

log = logging.getLogger("alliance.fetch.gene_descriptions")


DEFAULT_MODS = ["MGI", "ZFIN", "FB", "WB", "RGD", "SGD", "HUMAN", "XBXL", "XBXT"]


COLUMNS = [
    "geneId",
    "geneName",
    "autoDescription",
    "goDescription",
    "goFunctionDescription",
    "goProcessDescription",
    "goComponentDescription",
    "doDescription",
    "expressionDescription",
    "orthologyDescription",
    "dataProvider",
]


def _row(rec: dict, mod: str) -> dict | None:
    gene_id = rec.get("gene_id")
    if not gene_id:
        return None
    return {
        "geneId": gene_id,
        "geneName": rec.get("gene_name") or "",
        "autoDescription": rec.get("description") or "",
        "goDescription": rec.get("go_description") or "",
        "goFunctionDescription": rec.get("go_function_description") or "",
        "goProcessDescription": rec.get("go_process_description") or "",
        "goComponentDescription": rec.get("go_component_description") or "",
        "doDescription": rec.get("do_description") or "",
        "expressionDescription": rec.get("tissue_expression_description") or "",
        "orthologyDescription": rec.get("orthology_description") or "",
        "dataProvider": mod,
    }


def _iter_mod(mod: str) -> Iterator[dict]:
    for sub, rec in iter_fms_datatype("GENE-DESCRIPTION-JSON", sub_filter=mod):
        row = _row(rec, mod)
        if row is not None:
            yield row


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance gene descriptions (FMS bulk)")
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
        out_dir / "gene-descriptions.tsv",
        COLUMNS,
        release=release,
        source=f"FMS GENE-DESCRIPTION-JSON ({','.join(mods)})",
    ) as writer:
        for mod in mods:
            log.info("Gene descriptions FMS: %s", mod)
            for row in _iter_mod(mod):
                if row["geneId"] in seen:
                    continue
                seen.add(row["geneId"])
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    log.info("Hit limit; stopping")
                    log.info("Done. rows=%d", total)
                    return

    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
