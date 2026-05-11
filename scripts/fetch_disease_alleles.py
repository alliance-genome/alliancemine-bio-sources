#!/usr/bin/env python3
"""
Fetch per-disease allele annotations from the Alliance API.

For each disease (DOID), call /disease/{id}/alleles paginated. Emits partial
DiseaseAnnotation items keyed on (alleleSubject, ontologyTerm) that the
integration engine merges into the AGM/gene-shape annotations from sibling
sources.

Seed list: same as fetch_disease_annotations.py — DISEASE-ALLIANCE_COMBINED
FMS export DOID column.

Usage:
    python3 fetch_disease_alleles.py
    python3 fetch_disease_alleles.py --ids DOID:2841
    python3 fetch_disease_alleles.py --limit 100

Output: data/disease-alleles.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    get_current_release,
    get_nested,
    join_pipe,
    normalize_taxon,
    open_cache,
    paginate,
)

log = logging.getLogger("alliance.fetch.disease_alleles")


COLUMNS = [
    "diseaseId",
    "diseaseName",
    "alleleId",
    "alleleSymbol",
    "alleleTaxon",
    "relationName",
    "evidenceCodes",
    "evidencePmids",
    "evidenceCurie",
    "negated",
    "dataProvider",
]


def _enumerate_diseases_from_combined() -> list[str]:
    url = "https://fms.alliancegenome.org/download/DISEASE-ALLIANCE_COMBINED.tsv.gz"
    raw = urllib.request.urlopen(url, timeout=120).read()
    text = gzip.decompress(raw).decode("utf-8", errors="replace")
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    if not lines:
        return []
    header = lines[0].split("\t")
    try:
        idx = header.index("DOID")
    except ValueError:
        return []
    seen: set = set()
    out: list[str] = []
    for row in lines[1:]:
        cols = row.split("\t")
        if len(cols) > idx:
            d = cols[idx].strip()
            if d and d not in seen:
                seen.add(d)
                out.append(d)
    return out


def _row(disease_id: str, result: dict) -> dict | None:
    obj = result.get("object") or {}
    primary = (result.get("primaryAnnotations") or [{}])[0]
    subj = primary.get("diseaseAnnotationSubject") or {}
    allele_id = get_nested(subj, "primaryExternalId", default="")
    if not allele_id:
        return None
    eco_codes = join_pipe(
        get_nested(c, "curie", default="")
        for c in (result.get("evidenceCodes") or [])
        if get_nested(c, "curie")
    )
    pmids = join_pipe(
        get_nested(r, "referenceID", default="")
        for r in (result.get("references") or [])
        if get_nested(r, "referenceID")
    )
    return {
        "diseaseId": get_nested(obj, "curie", default=disease_id),
        "diseaseName": get_nested(obj, "name", default=""),
        "alleleId": allele_id,
        "alleleSymbol": get_nested(subj, "alleleSymbol", "displayText", default=""),
        "alleleTaxon": normalize_taxon(get_nested(subj, "taxon", "curie", default="")),
        "relationName": get_nested(result, "relation", "name", default=""),
        "evidenceCodes": eco_codes,
        "evidencePmids": pmids,
        "evidenceCurie": get_nested(primary, "evidenceItem", "curie", default=""),
        "negated": str(primary.get("negated", "")).lower() if primary.get("negated") is not None else "",
        "dataProvider": get_nested(primary, "dataProvider", "abbreviation", default=""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance per-disease alleles into TSV")
    parser.add_argument("--ids", help="Comma-separated DOID CURIEs")
    parser.add_argument("--limit", type=int, help="Process at most N diseases")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="ThreadPool size (default 4; keep low to avoid AGR WAF rate limit)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    if args.ids:
        disease_ids = [d.strip() for d in args.ids.split(",") if d.strip()]
    else:
        disease_ids = _enumerate_diseases_from_combined()
    if args.limit:
        disease_ids = disease_ids[: args.limit]

    log.info("Processing %d diseases", len(disease_ids))
    release = get_current_release()
    out_dir = Path(args.out_dir)

    stats = FetchStats()
    seen_pairs: set = set()

    def fetch_one(did: str) -> tuple[str, list[dict], Exception | None]:
        rows: list[dict] = []
        try:
            with open_cache("disease_alleles") as cache:
                for result in paginate(
                    f"/disease/{did}/alleles", cache=cache, stats=stats
                ):
                    row = _row(did, result)
                    if row is not None:
                        rows.append(row)
            return did, rows, None
        except Exception as e:
            return did, [], e

    with TsvWriter(
            out_dir / "disease-alleles.tsv",
            COLUMNS,
            release=release,
            source="/disease/*/alleles",
         ) as writer:
        total_rows = 0
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(fetch_one, did): did for did in disease_ids}
            try:
                from tqdm import tqdm
                pbar = tqdm(total=len(futures), desc="diseases", unit="dz")
            except ImportError:
                pbar = None
            for fut in as_completed(futures):
                did, rows, err = fut.result()
                if pbar:
                    pbar.update(1)
                if err is not None:
                    log.warning("Disease-alleles fetch failed for %s: %s", did, err)
                    continue
                for row in rows:
                    pair_key = (row["alleleId"], row["diseaseId"])
                    if pair_key in seen_pairs:
                        continue
                    seen_pairs.add(pair_key)
                    writer.write_row(**row)
                    total_rows += 1
            if pbar:
                pbar.close()

    log.info("Done. rows=%d  API stats: %s", total_rows, stats)


if __name__ == "__main__":
    main()
