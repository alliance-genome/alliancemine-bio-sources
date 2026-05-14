#!/usr/bin/env python3
"""
Fetch CRISPR/RNAi screen hits from the Alliance FMS BIOGRID-ORCS bulk
tarballs. Each tarball contains hundreds of per-screen tab files; the
HUMAN tarball is ~750MB compressed. To avoid blowing disk, the fetcher
streams the tarball over HTTP via tarfile.open(fileobj=...) and
emits only HIT=YES rows to the output TSV.

Per-row schema (after the leading #SCREEN_ID comment header):
    SCREEN_ID IDENTIFIER_ID IDENTIFIER_TYPE OFFICIAL_SYMBOL ALIASES
    ORGANISM_ID ORGANISM_OFFICIAL SCORE.1 SCORE.2 SCORE.3 SCORE.4 SCORE.5
    HIT SOURCE

Usage:
    python3 fetch_crispr_screens.py
    python3 fetch_crispr_screens.py --mods HUMAN --limit 1000
    python3 fetch_crispr_screens.py --keep-all   # drop HIT-only filter

Output: data/crispr-screens.tsv
"""

from __future__ import annotations

import argparse
import io
import logging
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    fms_list_files,
    get_current_release,
    normalize_taxon,
)

log = logging.getLogger("alliance.fetch.crispr_screens")


DEFAULT_MODS = ["HUMAN", "MGI", "FB"]
DATATYPE = "BIOGRID-ORCS"


COLUMNS = [
    "screenId",
    "geneId",          # OFFICIAL_SYMBOL (BIOGRID-ORCS uses symbol, not curie)
    "geneSymbol",
    "score1",
    "score2",
    "hit",
    "taxon",
    "source",
    "dataProvider",
]


def _stream_open(url: str) -> tarfile.TarFile:
    """Open a tarball over HTTP without writing it to disk. tarfile.open
    accepts a streaming file-like object as long as it's read sequentially."""
    log.info("Opening %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": "AllianceMine-fetcher/1.0"})
    response = urllib.request.urlopen(req, timeout=600)
    # tarfile needs a seekable-or-stream object; passing the response works
    # because we read members in archive order without seeking.
    return tarfile.open(fileobj=response, mode="r|gz")


def _iter_orcs_rows(tar: tarfile.TarFile, keep_all: bool) -> Iterator[list[str]]:
    for member in tar:
        if not member.isfile() or not member.name.endswith(".screen.tab.txt"):
            continue
        fh = tar.extractfile(member)
        if fh is None:
            continue
        for raw in fh:
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 14:
                continue
            if not keep_all and cols[12].strip().upper() != "YES":
                continue
            yield cols


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch CRISPR/RNAi screen hits from FMS BIOGRID-ORCS")
    parser.add_argument(
        "--mods",
        default=",".join(DEFAULT_MODS),
        help="Comma-separated MOD prefixes (HUMAN, MGI, FB available)",
    )
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="Keep every row (default: HIT=YES filter only)",
    )
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
        out_dir / "crispr-screens.tsv",
        COLUMNS,
        release=release,
        source=f"FMS BIOGRID-ORCS ({','.join(mods)})",
    ) as writer:
        for mod in mods:
            files = fms_list_files(DATATYPE, sub_filter=mod)
            for fi in files:
                tar = _stream_open(fi["url"])
                try:
                    for cols in _iter_orcs_rows(tar, args.keep_all):
                        screen_id = cols[0].strip()
                        gene_sym = cols[3].strip()
                        if not screen_id or not gene_sym:
                            continue
                        key = (screen_id, gene_sym)
                        if key in seen:
                            continue
                        seen.add(key)
                        writer.write_row(
                            screenId=screen_id,
                            geneId=gene_sym,
                            geneSymbol=gene_sym,
                            score1=cols[7].strip(),
                            score2=cols[8].strip(),
                            hit=cols[12].strip().lower() if cols[12].strip() != "-" else "",
                            taxon=normalize_taxon(f"NCBITaxon:{cols[5].strip()}") if cols[5].strip() and cols[5].strip() != "-" else "",
                            source=cols[13].strip() if len(cols) > 13 else "",
                            dataProvider=mod,
                        )
                        total += 1
                        if args.limit is not None and total >= args.limit:
                            log.info("Hit limit; stopping")
                            return
                finally:
                    tar.close()
    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
