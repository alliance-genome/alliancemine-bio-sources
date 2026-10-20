#!/usr/bin/env python3
"""
Fetch per-allele detail from the Alliance API and emit a companion TSV that
the existing alliance-alleles converter joins to its FMS VARIANT-ALLELE rows.

Seed list: AlleleId of the pinned release's per-taxon VARIANT-ALLELE TSVs
(the COMBINED file is frozen at 4.0.0; in 9.x the HUMAN shard is also 4.0.0
and XBXT 8.1.0, which is what the release ships).

The /allele/{id} endpoint is single-resource (no pagination) and returns
slim metadata for yeast alleles plus richer fields for mouse/fly alleles.
We capture three fields as Allele attribute additions: alterationType,
apiCategory, apiCrossReference.

Usage:
    python3 fetch_allele_detail.py                    # alleles of all 9 MODs
    python3 fetch_allele_detail.py --mods SGD,WB
    python3 fetch_allele_detail.py --ids SGD:S000279149
    python3 fetch_allele_detail.py --from-tsv path/to/VARIANT-ALLELE.tsv
    python3 fetch_allele_detail.py --limit 100

Output: data/allele-detail.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    add_parallel_args,
    configure_logging,
    enumerate_allele_ids,
    finish_run,
    get_current_release,
    get_nested,
    http_get_json,
    run_parallel,
    selected_mods,
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
    parser.add_argument("--mods", help="Comma-separated MODs whose alleles to seed (default: all 9, env ALLIANCE_MODS)")
    parser.add_argument("--limit", type=int, help="Process at most N alleles")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    add_parallel_args(parser)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    if args.ids:
        allele_ids = [a.strip() for a in args.ids.split(",") if a.strip()]
    elif args.from_tsv:
        allele_ids = _enumerate_alleles_from_local_tsv(Path(args.from_tsv))
    else:
        log.info("Enumerating allele IDs from FMS VARIANT-ALLELE (release %s)...", release)
        allele_ids = enumerate_allele_ids(selected_mods(None, args.mods))
    if args.limit:
        allele_ids = allele_ids[: args.limit]

    out_dir = Path(args.out_dir)
    stats = FetchStats()

    def work(aid, cache, stats):
        return _row(aid, http_get_json(f"/allele/{aid}", cache=cache, stats=stats))

    with TsvWriter(
        out_dir / "allele-detail.tsv",
        COLUMNS,
        release=release,
        source="/allele/{id}",
    ) as writer:
        report = run_parallel(
            allele_ids, work, lambda _aid, row: writer.write_row(**row),
            cache_name="allele_detail", stats=stats, workers=args.workers, label="alleles",
        )
        finish_run(report, stats, out_dir, args.max_failure_rate,
                   rows_written=writer.rows_written)


if __name__ == "__main__":
    main()
