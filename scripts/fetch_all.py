#!/usr/bin/env python3
"""
Orchestrator for every Alliance API fetcher in this directory.

Runs each fetch_*.py in the canonical order that the downstream Gradle build
will ingest them. Designed to be called from the AllianceMine Docker pipeline
as a single pre-build step replacing the old FMS pull.

Usage:
    python3 fetch_all.py                  # full run, all fetchers
    python3 fetch_all.py --only genes,paralogs
    python3 fetch_all.py --skip phenotypes
    python3 fetch_all.py --limit 50       # smoke test — applies to every fetcher

Environment variables:
    ALLIANCE_API_BASE   - API root (default https://www.alliancegenome.org/api)
    ALLIANCE_FETCH_CACHE - SQLite cache dir (default scripts/.cache/)
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import configure_logging  # noqa: E402

log = logging.getLogger("alliance.fetch.orchestrator")


# Run order matters: genes first, since everything else joins to gene primary
# identifiers; interactions next (biggest volume, longest runtime); then the
# smaller feeds and the *-detail enrichment passes.
FETCHERS = [
    ("genes",               "fetch_genes.py"),
    ("interactions",        "fetch_interactions.py"),
    ("orthologs",           "fetch_orthologs.py"),
    ("paralogs",            "fetch_paralogs.py"),
    ("allele_detail",       "fetch_allele_detail.py"),
    ("disease_annotations", "fetch_disease_annotations.py"),
    ("disease_models",      "fetch_disease_models.py"),
    ("phenotypes",          "fetch_phenotypes.py"),
    ("transgenic_alleles",  "fetch_transgenic_alleles.py"),
    ("disease_alleles",     "fetch_disease_alleles.py"),
    ("experimental_disease", "fetch_experimental_disease.py"),
    # Phase 6d: cross-mine federation via InterMine PathQuery REST.
    ("mousemine_strains",   "fetch_mousemine_strains.py"),
    ("wormmine_rnai",       "fetch_wormmine_rnai.py"),
]


def run_one(name: str, script: str, extra_args: list[str]) -> tuple[str, int, float]:
    script_path = Path(__file__).parent / script
    cmd = [sys.executable, str(script_path)] + extra_args
    log.info("=> %s  (%s)", name, " ".join(cmd))
    start = time.time()
    proc = subprocess.run(cmd)
    elapsed = time.time() - start
    status = "ok" if proc.returncode == 0 else f"FAILED rc={proc.returncode}"
    log.info("<= %s  %s  in %.1fs", name, status, elapsed)
    return name, proc.returncode, elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description="Run every Alliance fetcher in order")
    parser.add_argument("--only", help="Comma-separated subset of fetcher names to run")
    parser.add_argument("--skip", help="Comma-separated subset to skip")
    parser.add_argument("--limit", type=int, help="Pass --limit N through to each fetcher")
    parser.add_argument("--out-dir", help="Pass --out-dir through to each fetcher")
    parser.add_argument("--ids", help="Pass --ids through to each fetcher (smoke tests)")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true",
                        help="Abort immediately on any fetcher failure")
    args = parser.parse_args()

    configure_logging(args.verbose)

    only = set(n.strip() for n in args.only.split(",")) if args.only else None
    skip = set(n.strip() for n in args.skip.split(",")) if args.skip else set()

    extra: list[str] = []
    if args.limit:
        extra += ["--limit", str(args.limit)]
    if args.out_dir:
        extra += ["--out-dir", args.out_dir]
    if args.ids:
        extra += ["--ids", args.ids]
    if args.verbose:
        extra += ["--verbose"]

    results = []
    for name, script in FETCHERS:
        if only and name not in only:
            continue
        if name in skip:
            log.info("-- skipping %s", name)
            continue
        result = run_one(name, script, extra)
        results.append(result)
        if args.stop_on_error and result[1] != 0:
            log.error("Stopping: %s failed", name)
            break

    # Summary
    log.info("")
    log.info("=== Summary ===")
    any_failed = False
    for name, rc, elapsed in results:
        status = "ok" if rc == 0 else f"FAILED (rc={rc})"
        log.info("  %-15s  %s  %.1fs", name, status, elapsed)
        if rc != 0:
            any_failed = True
    sys.exit(1 if any_failed else 0)


if __name__ == "__main__":
    main()
