#!/usr/bin/env python3
"""
Fetch per-allele detail from the Alliance API and emit a companion TSV that
the existing alliance-alleles converter joins to its FMS VARIANT-ALLELE rows.

The /allele/{id} endpoint is single-resource (no pagination) and returns
slim metadata for yeast alleles plus richer fields for mouse/fly alleles.
We capture three fields as Allele attribute additions: alterationType,
apiCategory, apiCrossReference.

Usage:
    python3 fetch_allele_detail.py                    # all yeast alleles
    python3 fetch_allele_detail.py --ids SGD:S000279149
    python3 fetch_allele_detail.py --from-tsv path/to/VARIANT-ALLELE.tsv
    python3 fetch_allele_detail.py --limit 100

Output: data/allele-detail.tsv
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
    FetchStats,
    TsvWriter,
    configure_logging,
    get_current_release,
    get_nested,
    http_get_json,
    open_cache,
)

log = logging.getLogger("alliance.fetch.allele_detail")


COLUMNS = [
    "alleleId",
    "alterationType",
    "apiCategory",
    "apiCrossReference",
]


def _row(allele_id: str, payload: dict) -> dict:
    return {
        "alleleId": allele_id,
        "alterationType": payload.get("alterationType", ""),
        "apiCategory": payload.get("category", ""),
        "apiCrossReference": get_nested(payload, "crossReference", "displayName", default=""),
    }


def _enumerate_alleles_from_combined_tsv(url: str = "https://fms.alliancegenome.org/download/VARIANT-ALLELE_COMBINED.tsv.gz") -> list[str]:
    """Pull AlleleId column from the FMS COMBINED TSV (while FMS lasts)."""
    raw = urllib.request.urlopen(url, timeout=120).read()
    text = gzip.decompress(raw).decode("utf-8", errors="replace")
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    if not lines:
        return []
    header = lines[0].split("\t")
    try:
        allele_idx = header.index("AlleleId")
    except ValueError:
        log.error("VARIANT-ALLELE TSV missing AlleleId column")
        return []
    seen: set = set()
    out: list[str] = []
    for row in lines[1:]:
        cols = row.split("\t")
        if len(cols) > allele_idx:
            aid = cols[allele_idx].strip()
            if aid and aid not in seen:
                seen.add(aid)
                out.append(aid)
    return out


def _enumerate_alleles_from_local_tsv(path: Path) -> list[str]:
    seen: set = set()
    out: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0] == "AlleleId":  # header
                continue
            if len(cols) >= 3:
                aid = cols[2].strip()
                if aid and aid not in seen:
                    seen.add(aid)
                    out.append(aid)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance per-allele detail")
    parser.add_argument("--ids", help="Comma-separated allele CURIEs (skips enumeration)")
    parser.add_argument("--from-tsv", help="Path to a local VARIANT-ALLELE TSV to seed allele IDs from")
    parser.add_argument("--limit", type=int, help="Process at most N alleles")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    if args.ids:
        allele_ids = [a.strip() for a in args.ids.split(",") if a.strip()]
    elif args.from_tsv:
        allele_ids = _enumerate_alleles_from_local_tsv(Path(args.from_tsv))
    else:
        log.info("Enumerating allele IDs from FMS VARIANT-ALLELE_COMBINED...")
        allele_ids = _enumerate_alleles_from_combined_tsv()
    if args.limit:
        allele_ids = allele_ids[: args.limit]

    log.info("Processing %d alleles", len(allele_ids))
    release = get_current_release()
    out_dir = Path(args.out_dir)

    try:
        from tqdm import tqdm
        iterator = tqdm(allele_ids, desc="alleles", unit="allele")
    except ImportError:
        iterator = allele_ids

    stats = FetchStats()
    with open_cache("allele_detail") as cache, \
         TsvWriter(
             out_dir / "allele-detail.tsv",
             COLUMNS,
             release=release,
             source="/allele/{id}",
         ) as writer:
        for aid in iterator:
            try:
                # URL-encode the colon to avoid path issues
                payload = http_get_json(f"/allele/{aid}", cache=cache, stats=stats)
            except Exception as e:
                log.warning("Allele detail fetch failed for %s: %s", aid, e)
                continue
            writer.write_row(**_row(aid, payload))

    log.info("Done. API stats: %s", stats)


if __name__ == "__main__":
    main()
