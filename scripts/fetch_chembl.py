#!/usr/bin/env python3
"""
Fetch ChEMBL target ↔ UniProt mapping. Single tiny file (~700KB,
~17k rows) from EBI FTP. Each row maps a ChEMBL target ID to a UniProt
accession with a human-readable name and target type.

For drug-target enrichment beyond the bare mapping (i.e. drug→target
activity / mechanism / max_phase) the ChEMBL SQLite dump is required
(5.6GB). This fetcher stays focused on the mapping so AllianceMine
proteins gain a ChEMBL cross-reference; downstream tools can pivot from
the ChEMBL ID via ChEMBL's web service or local SQLite when needed.

Input: chembl_uniprot_mapping.txt (no filter — we keep all 17k rows;
       the integration engine drops xrefs that don't bind to a Protein)
Output: data/chembl-targets.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import TsvWriter, configure_logging, get_current_release  # noqa: E402

log = logging.getLogger("alliance.fetch.chembl")

CHEMBL_MAPPING_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_uniprot_mapping.txt"
CHEMBL_TARGET_BASE = "https://www.ebi.ac.uk/chembl/target_report_card/"


COLUMNS = [
    "uniProtAccession",
    "chemblTargetId",
    "targetName",
    "targetType",
    "url",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch ChEMBL target → UniProt mapping")
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)
    out_dir = Path(args.out_dir)
    release = get_current_release()
    log.info("Downloading %s", CHEMBL_MAPPING_URL)
    req = urllib.request.Request(CHEMBL_MAPPING_URL, headers={"User-Agent": "AllianceMine-fetcher/1.0"})
    body = urllib.request.urlopen(req, timeout=120).read().decode("utf-8", errors="replace")
    total = 0
    with TsvWriter(
        out_dir / "chembl-targets.tsv",
        COLUMNS,
        release=release,
        source=CHEMBL_MAPPING_URL,
    ) as writer:
        for line in body.splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            acc, chembl_id, name, ttype = cols[0], cols[1], cols[2], cols[3]
            if not acc or not chembl_id:
                continue
            writer.write_row(
                uniProtAccession=acc,
                chemblTargetId=chembl_id,
                targetName=name,
                targetType=ttype,
                url=f"{CHEMBL_TARGET_BASE}{chembl_id}",
            )
            total += 1
            if args.limit is not None and total >= args.limit:
                break
    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
