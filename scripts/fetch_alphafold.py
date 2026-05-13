#!/usr/bin/env python3
"""
Synthesize AlphaFold structure URLs for every UniProt accession surfaced
by the Alliance gene-crossrefs feed. No HTTP — pure URL templating
against the public EBI AlphaFold viewer
(https://alphafold.ebi.ac.uk/entry/{accession}). Coverage is near-100%
for AGR-target taxa since AlphaFoldDB v4 covers the entire UniProt
proteome.

Input: data/gene-crossrefs.tsv (must exist; produced by
       scripts/fetch_gene_crossrefs.py)
Output: data/alphafold.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import TsvWriter, configure_logging, get_current_release  # noqa: E402

log = logging.getLogger("alliance.fetch.alphafold")

ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/entry/"


COLUMNS = [
    "uniProtAccession",
    "geneId",
    "geneTaxon",
    "alphaFoldUrl",
]


def _rows_from_crossrefs(path: Path):
    """Yield (acc, geneId, taxon, url) per UniProtKB xref."""
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0] == "geneId":
                continue
            if len(cols) < 5:
                continue
            gene_id, taxon, xref_id, xref_url, xref_type = cols[:5]
            if not xref_id.startswith("UniProtKB:"):
                continue
            acc = xref_id.split(":", 1)[1]
            yield {
                "uniProtAccession": acc,
                "geneId": gene_id,
                "geneTaxon": taxon,
                "alphaFoldUrl": f"{ALPHAFOLD_URL}{acc}",
            }


def main() -> None:
    parser = argparse.ArgumentParser(description="Emit AlphaFold URLs from gene-crossrefs UniProtKB rows")
    parser.add_argument(
        "--crossrefs",
        default=str(Path(__file__).parent.parent / "data" / "gene-crossrefs.tsv"),
        help="Path to gene-crossrefs.tsv (defaults to ../data/gene-crossrefs.tsv)",
    )
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)
    crossrefs = Path(args.crossrefs)
    if not crossrefs.exists():
        log.error("Missing %s — run scripts/fetch_gene_crossrefs.py first", crossrefs)
        sys.exit(1)

    out_dir = Path(args.out_dir)
    release = get_current_release()
    seen: set = set()
    total = 0
    with TsvWriter(
        out_dir / "alphafold.tsv",
        COLUMNS,
        release=release,
        source="alphafold.ebi.ac.uk (URL template; no HTTP fetch)",
    ) as writer:
        for row in _rows_from_crossrefs(crossrefs):
            if row["uniProtAccession"] in seen:
                continue
            seen.add(row["uniProtAccession"])
            writer.write_row(**row)
            total += 1
            if args.limit is not None and total >= args.limit:
                break
    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
