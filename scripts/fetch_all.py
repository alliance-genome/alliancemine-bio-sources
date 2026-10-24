#!/usr/bin/env python3
"""
Orchestrator for every Alliance API fetcher in this directory.

Runs each fetch_*.py in the canonical order that the downstream Gradle build
will ingest them. Designed to be called from the AllianceMine Docker pipeline
as a single pre-build step replacing the old FMS pull.

Usage:
    ALLIANCE_RELEASE=9.2.0 python3 fetch_all.py --out-dir /root/data/api
    python3 fetch_all.py --release 9.2.0 --workers 6 --out-dir /root/data/api
    python3 fetch_all.py --only genes,paralogs
    python3 fetch_all.py --skip phenotypes
    python3 fetch_all.py --mods SGD,MGI --limit 50 --allow-release-mismatch   # smoke test

Before running anything the pinned release (--release or ALLIANCE_RELEASE)
is compared with what the API serves (/api/releaseInfo); a mismatch exits 2
unless --allow-release-mismatch. Unpinned runs pin to the API's release.
The pin is exported to every fetcher so FMS file selection, gene seeding and
TSV provenance all use the same release.

Pass-through flags (--limit, --out-dir, --workers, --max-failure-rate,
--verbose) go only to fetchers whose --help lists them; --ids goes only to
the fetchers seeded by gene CURIEs. --mods is exported as ALLIANCE_MODS,
which narrows each fetcher's own default MOD list.

Environment variables:
    ALLIANCE_RELEASE     - pinned Alliance release (e.g. 9.2.0)
    ALLIANCE_MODS        - comma-separated MOD subset (default: all 9)
    ALLIANCE_API_BASE    - API root (default https://www.alliancegenome.org/api)
    ALLIANCE_FETCH_CACHE - SQLite cache + gene-list dir (default scripts/.cache/)
    ALLIANCE_FMS_CACHE   - downloaded FMS files (default scripts/.fms_cache/)
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import configure_logging, get_api_release, parse_mods  # noqa: E402

log = logging.getLogger("alliance.fetch.orchestrator")


# Run order:
#   1. Primary identifiers (genes) first — all downstream feeds join on Gene.
#   2. Per-gene API enrichment (interactions, orthologs, paralogs, phenotypes,
#      disease/models/annotations) — these still use the live API (per-entity).
#   3. FMS-bulk loads (allele/transgenic/disease-alleles/experimental-disease/
#      variants/gene-descriptions/agms/allele-phenotypes/htp). Default source
#      is FMS where the fetcher accepts --source.
#   4. Derived / external (gene_crossrefs → alphafold; string; crispr_screens;
#      molecules; chembl). gene_crossrefs must run BEFORE alphafold because
#      alphafold synthesises URLs from the UniProtKB rows it produces.
#   5. Cross-mine PathQuery federation (mousemine, wormmine).
FETCHERS = [
    # Primary
    ("genes",                "fetch_genes.py"),

    # API-driven per-gene enrichment
    ("interactions",         "fetch_interactions.py"),
    ("orthologs",            "fetch_orthologs.py"),
    ("paralogs",             "fetch_paralogs.py"),
    ("allele_detail",        "fetch_allele_detail.py"),
    ("disease_annotations",  "fetch_disease_annotations.py"),
    ("disease_models",       "fetch_disease_models.py"),
    ("phenotypes",           "fetch_phenotypes.py"),

    # FMS bulk loads (default --source fms; some take a long time)
    ("alleles",              "fetch_alleles.py"),
    ("disease",              "fetch_disease.py"),
    ("expression",           "fetch_expression.py"),
    ("transgenic_alleles",   "fetch_transgenic_alleles.py"),
    ("disease_alleles",      "fetch_disease_alleles.py"),
    ("experimental_disease", "fetch_experimental_disease.py"),
    ("variants",             "fetch_variants.py"),
    ("gene_descriptions",    "fetch_gene_descriptions.py"),
    ("agms",                 "fetch_agms.py"),
    ("allele_phenotypes",    "fetch_allele_phenotypes.py"),
    ("htp",                  "fetch_htp.py"),

    # Cross-reference + external (gene_crossrefs feeds alphafold)
    ("gene_crossrefs",       "fetch_gene_crossrefs.py"),
    ("alphafold",            "fetch_alphafold.py"),

    # External data sources (independent of MOD pipeline)
    ("string",               "fetch_string.py"),
    ("crispr_screens",       "fetch_crispr_screens.py"),
    ("molecules",            "fetch_molecules.py"),
    ("chembl",               "fetch_chembl.py"),

    # Cross-mine PathQuery REST federation
    ("mousemine_strains",    "fetch_mousemine_strains.py"),
    ("wormmine_rnai",        "fetch_wormmine_rnai.py"),
]

# With --source db these names run extract_<name>.py (SQL against a restored
# curation dump, see curation_db.py) in place of fetch_<name>.py. Names absent
# here, or whose extractor file does not exist yet, keep their fetcher.
# htp keeps FMS: 39,176 RGD samples have no dataset link in the 9.2.0 curation DB.
KEEP_FETCHER = {"htp"}


def script_for(name: str, script: str, source: str) -> str:
    extractor = f"extract_{name}.py"
    if source == "db" and name not in KEEP_FETCHER and (Path(__file__).parent / extractor).exists():
        return extractor
    return script


# Fetchers whose --ids takes gene CURIEs (allele_detail takes allele ids,
# disease_annotations / disease_alleles take DOIDs).
GENE_ID_FETCHERS = {
    "genes", "interactions", "orthologs", "paralogs", "disease_models", "phenotypes",
    "transgenic_alleles", "experimental_disease",
}

_OPTION_RE = re.compile(r"^\s+(?:-\w, )?(--[\w-]+)", re.MULTILINE)


def accepted_flags(script: str) -> set[str]:
    """Options a fetcher declares, read from its own --help."""
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).parent / script), "--help"],
        capture_output=True, text=True, timeout=120,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{script} --help failed: {proc.stderr.strip()}")
    return set(_OPTION_RE.findall(proc.stdout))


def args_for(name: str, script: str, passthrough: dict[str, str | None]) -> list[str]:
    flags = accepted_flags(script)
    out: list[str] = []
    for flag, value in passthrough.items():
        if flag not in flags:
            continue
        if flag == "--ids" and name not in GENE_ID_FETCHERS:
            continue
        out += [flag] if value is None else [flag, value]
    return out


def check_release(pinned: str | None, allow_mismatch: bool) -> str:
    """Return the release to pin every fetcher to, or exit when it does not
    match the release the API is serving."""
    api_release = get_api_release()
    if not pinned:
        if not api_release:
            log.error("No ALLIANCE_RELEASE / --release given and /api/releaseInfo is unreachable")
            sys.exit(2)
        log.warning("No release pinned; pinning to the API's release %s", api_release)
        return api_release
    if api_release == pinned:
        log.info("Release %s matches /api/releaseInfo", pinned)
        return pinned
    message = (f"Pinned release {pinned} but the API serves {api_release or 'unknown (unreachable)'}; "
               "API-fetched data would not match the pinned FMS files")
    if not allow_mismatch:
        log.error("%s. Wait for the API switch-over or pass --allow-release-mismatch.", message)
        sys.exit(2)
    log.warning("%s (--allow-release-mismatch given)", message)
    return pinned


def run_one(name: str, script: str, extra_args: list[str], env: dict[str, str]) -> tuple[str, int, float]:
    script_path = Path(__file__).parent / script
    cmd = [sys.executable, str(script_path)] + extra_args
    log.info("=> %s  (%s)", name, " ".join(cmd))
    start = time.time()
    proc = subprocess.run(cmd, env=env)
    elapsed = time.time() - start
    status = "ok" if proc.returncode == 0 else f"FAILED rc={proc.returncode}"
    log.info("<= %s  %s  in %.1fs", name, status, elapsed)
    return name, proc.returncode, elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description="Run every Alliance fetcher in order")
    parser.add_argument("--only", help="Comma-separated subset of fetcher names to run")
    parser.add_argument("--skip", help="Comma-separated subset to skip")
    parser.add_argument("--release", help="Alliance release to pin (default: env ALLIANCE_RELEASE)")
    parser.add_argument("--allow-release-mismatch", action="store_true",
                        help="Run even when the pinned release differs from /api/releaseInfo (smoke tests)")
    parser.add_argument("--mods", help="Comma-separated MOD subset, exported as ALLIANCE_MODS")
    parser.add_argument("--limit", type=int, help="Pass --limit N through to each fetcher")
    parser.add_argument("--out-dir", help="Pass --out-dir through to each fetcher")
    parser.add_argument("--ids", help="Pass --ids through to the gene-seeded fetchers (smoke tests)")
    parser.add_argument("--workers", type=int, help="Pass --workers N through to the API fan-out fetchers")
    parser.add_argument("--max-failure-rate", type=float,
                        help="Pass --max-failure-rate through to the API fan-out fetchers")
    parser.add_argument("--source", choices=("api", "db"), default="api",
                        help="db: run extract_<name>.py against the curation dump where one exists")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true",
                        help="Abort immediately on any fetcher failure")
    args = parser.parse_args()

    configure_logging(args.verbose)

    only = set(n.strip() for n in args.only.split(",")) if args.only else None
    skip = set(n.strip() for n in args.skip.split(",")) if args.skip else set()
    known = {name for name, _ in FETCHERS}
    unknown = ((only or set()) | skip) - known
    if unknown:
        parser.error(f"unknown fetcher name(s): {', '.join(sorted(unknown))}")

    env = dict(os.environ)
    pinned = args.release or env.get("ALLIANCE_RELEASE", "").strip() or None
    env["ALLIANCE_RELEASE"] = check_release(pinned, args.allow_release_mismatch)
    if args.mods:
        env["ALLIANCE_MODS"] = ",".join(parse_mods(args.mods))
    log.info("Release %s, MODs %s", env["ALLIANCE_RELEASE"], env.get("ALLIANCE_MODS") or "all 9")

    passthrough: dict[str, str | None] = {}
    if args.limit:
        passthrough["--limit"] = str(args.limit)
    if args.out_dir:
        passthrough["--out-dir"] = args.out_dir
    if args.ids:
        passthrough["--ids"] = args.ids
    if args.workers:
        passthrough["--workers"] = str(args.workers)
    if args.max_failure_rate is not None:
        passthrough["--max-failure-rate"] = str(args.max_failure_rate)
    if args.verbose:
        passthrough["--verbose"] = None

    results = []
    for name, script in FETCHERS:
        if only and name not in only:
            continue
        if name in skip:
            log.info("-- skipping %s", name)
            continue
        script = script_for(name, script, args.source)
        result = run_one(name, script, args_for(name, script, passthrough), env)
        results.append(result)
        if args.stop_on_error and result[1] != 0:
            log.error("Stopping: %s failed", name)
            break

    log.info("")
    log.info("=== Summary (release %s) ===", env["ALLIANCE_RELEASE"])
    any_failed = False
    for name, rc, elapsed in results:
        status = "ok" if rc == 0 else f"FAILED (rc={rc})"
        log.info("  %-20s  %s  %.1fs", name, status, elapsed)
        if rc != 0:
            any_failed = True
    sys.exit(1 if any_failed else 0)


if __name__ == "__main__":
    main()
