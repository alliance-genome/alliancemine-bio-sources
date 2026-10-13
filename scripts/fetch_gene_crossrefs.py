#!/usr/bin/env python3
"""
Fetch gene cross-references from the Alliance FMS GENECROSSREFERENCEJSON
bulk files. Each record carries (MOD gene id) → (external db id + URL).
Emits one CrossReference row per record.

Output: data/gene-crossrefs.tsv
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

log = logging.getLogger("alliance.fetch.gene_crossrefs")


MOD_PREFIXES = {"MGI", "SGD", "ZFIN", "FB", "WB", "RGD", "HGNC", "XBXL", "XBXT", "XB"}


COLUMNS = [
    "geneId",
    "geneTaxon",
    "xrefId",
    "xrefUrl",
    "xrefType",
]


def _iter_files() -> Iterator[dict]:
    for _sub, rec in iter_fms_datatype("GENECROSSREFERENCEJSON"):
        gene_id = rec.get("GeneID", "")
        if not gene_id:
            continue
        prefix = gene_id.split(":", 1)[0]
        if prefix not in MOD_PREFIXES:
            continue
        xref_id = rec.get("GlobalCrossReferenceID", "")
        if not xref_id:
            continue
        yield {
            "geneId": gene_id,
            "geneTaxon": normalize_taxon(rec.get("TaxonID", "")),
            "xrefId": xref_id,
            "xrefUrl": rec.get("CrossReferenceCompleteURL", ""),
            "xrefType": rec.get("ResourceDescriptorPage", ""),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance gene cross-references (FMS bulk)")
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    out_dir = Path(args.out_dir)
    seen: set = set()
    total = 0

    with TsvWriter(
        out_dir / "gene-crossrefs.tsv",
        COLUMNS,
        release=release,
        source="FMS GENECROSSREFERENCEJSON",
    ) as writer:
        for row in _iter_files():
            key = (row["geneId"], row["xrefId"])
            if key in seen:
                continue
            seen.add(key)
            writer.write_row(**row)
            total += 1
            if args.limit is not None and total >= args.limit:
                log.info("Done. rows=%d", total)
                return

    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
