#!/usr/bin/env python3
"""
Fetch STRING DB v12 per-organism interaction network. Downloads two
files per taxon:

    {taxid}.protein.aliases.v12.0.txt.gz  — STRING ID → UniProt mapping
    {taxid}.protein.links.v12.0.txt.gz    — score-weighted protein pairs

Emits one row per (proteinA, proteinB) interaction whose combined score
crosses the threshold. The default 700 corresponds to STRING's
"high confidence" tier.

Usage:
    python3 fetch_string.py
    python3 fetch_string.py --taxa 4932,9606 --min-score 700

Output: data/string-interactions.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import urllib.request
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import TsvWriter, configure_logging, get_current_release  # noqa: E402

log = logging.getLogger("alliance.fetch.string")


STRING_BASE = "https://stringdb-downloads.org/download"
STRING_VERSION = "v12.0"

# Taxa shared with AllianceMine. STRING uses bare NCBI tax IDs.
DEFAULT_TAXA = ["4932", "9606", "10090", "10116", "6239", "7227", "7955"]
LEGACY_REMAP = {"4932": "559292"}  # match D4 normalization at emit


COLUMNS = [
    "stringIdA",
    "stringIdB",
    "uniProtAccA",
    "uniProtAccB",
    "combinedScore",
    "taxon",
]


def _cache_path(filename: str) -> Path:
    base = Path(__file__).parent / ".string_cache"
    base.mkdir(parents=True, exist_ok=True)
    return base / filename


def _download(url: str) -> Path:
    fname = url.rsplit("/", 1)[-1]
    target = _cache_path(fname)
    if target.exists() and target.stat().st_size > 0:
        return target
    log.info("STRING download %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": "AllianceMine-fetcher/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp, target.open("wb") as fh:
        while True:
            chunk = resp.read(1 << 16)
            if not chunk:
                break
            fh.write(chunk)
    return target


def _alias_map(taxid: str) -> dict[str, str]:
    """STRING ID → first-seen UniProt accession. Some STRING entries have
    multiple UniProt mappings (isoforms); first is canonical enough for
    cross-source merge via Protein.primaryAccession."""
    url = f"{STRING_BASE}/protein.aliases.{STRING_VERSION}/{taxid}.protein.aliases.{STRING_VERSION}.txt.gz"
    path = _download(url)
    out: dict[str, str] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sid, alias, source = parts[0], parts[1], parts[2]
            if source == "UniProt_AC" and sid not in out:
                out[sid] = alias
    return out


def _iter_links(taxid: str, min_score: int) -> Iterator[tuple[str, str, int]]:
    url = f"{STRING_BASE}/protein.links.{STRING_VERSION}/{taxid}.protein.links.{STRING_VERSION}.txt.gz"
    path = _download(url)
    with gzip.open(path, "rt") as fh:
        next(fh)  # header
        for line in fh:
            cols = line.rstrip("\n").split(" ")
            if len(cols) != 3:
                continue
            try:
                score = int(cols[2])
            except ValueError:
                continue
            if score < min_score:
                continue
            yield cols[0], cols[1], score


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch STRING DB v12 interactions")
    parser.add_argument(
        "--taxa",
        default=",".join(DEFAULT_TAXA),
        help="Comma-separated NCBI tax IDs (STRING-style, e.g. 4932,9606)",
    )
    parser.add_argument("--min-score", type=int, default=700, help="Minimum combined_score (0-1000); default 700")
    parser.add_argument("--limit", type=int, help="Cap rows emitted (across all taxa)")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)
    taxa = [t.strip() for t in args.taxa.split(",") if t.strip()]
    release = get_current_release()
    out_dir = Path(args.out_dir)
    total = 0
    with TsvWriter(
        out_dir / "string-interactions.tsv",
        COLUMNS,
        release=release,
        source=f"STRING DB {STRING_VERSION} ({','.join(taxa)})",
    ) as writer:
        for taxid in taxa:
            log.info("STRING %s: building alias map", taxid)
            aliases = _alias_map(taxid)
            log.info("STRING %s: %d STRING→UniProt mappings", taxid, len(aliases))
            canonical_taxon = LEGACY_REMAP.get(taxid, taxid)
            emitted = 0
            for a, b, score in _iter_links(taxid, args.min_score):
                up_a = aliases.get(a, "")
                up_b = aliases.get(b, "")
                writer.write_row(
                    stringIdA=a,
                    stringIdB=b,
                    uniProtAccA=up_a,
                    uniProtAccB=up_b,
                    combinedScore=str(score),
                    taxon=f"NCBITaxon:{canonical_taxon}",
                )
                emitted += 1
                total += 1
                if args.limit is not None and total >= args.limit:
                    log.info("Done. rows=%d (taxon %s contributed %d)", total, taxid, emitted)
                    return
            log.info("STRING %s: emitted %d rows", taxid, emitted)

    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
