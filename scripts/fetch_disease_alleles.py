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
import logging
import sys
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    MOD_TAXON,
    FetchStats,
    TsvWriter,
    add_parallel_args,
    configure_logging,
    enumerate_disease_ids,
    finish_run,
    get_current_release,
    get_nested,
    iter_fms_datatype,
    join_pipe,
    normalize_taxon,
    paginate,
    run_parallel,
    selected_mods,
)

DEFAULT_MODS = ["MGI", "ZFIN", "FB", "WB", "RGD", "SGD", "HUMAN", "XBXL", "XBXT"]


def _fms_rows(mods: list[str]) -> Iterator[dict]:
    """Yield TSV-row dicts for each DAF record with objectType=='allele'.
    Disease + evidence + PMID + ECO are all in-line on the DAF record."""
    for mod in mods:
        log.info("FMS DAF extraction: %s", mod)
        for _sub, rec in iter_fms_datatype("DAF", sub_filter=mod):
            rel = rec.get("objectRelation") or {}
            if rel.get("objectType") != "allele":
                continue
            allele_id = rec.get("objectId", "")
            disease_id = rec.get("DOid", "")
            if not allele_id or not disease_id:
                continue
            evidence = rec.get("evidence") or {}
            eco_codes = join_pipe(evidence.get("evidenceCodes") or [])
            pub = evidence.get("publication") or {}
            pmid = pub.get("publicationId", "")
            yield {
                "diseaseId": disease_id,
                "diseaseName": "",                          # DAF lacks DOID name; converter keys by DOid alone
                "alleleId": allele_id,
                "alleleSymbol": rec.get("objectName", ""),  # DAF objectName == subject's symbol (allele here)
                "alleleTaxon": normalize_taxon(MOD_TAXON[mod]),
                "relationName": rel.get("associationType", ""),
                "evidenceCodes": eco_codes,
                "evidencePmids": pmid,
                "evidenceCurie": "",
                "negated": "",
                "dataProvider": mod,
            }

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
    parser.add_argument(
        "--source",
        choices=["fms", "api"],
        default="fms",
        help="Backend: 'fms' (bulk DAF JSON.gz from download.alliancegenome.org; default) or "
             "'api' (per-disease /disease/{id}/alleles).",
    )
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes, FMS path (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--ids", help="Comma-separated DOID CURIEs (api path)")
    parser.add_argument("--limit", type=int, help="Process at most N diseases / rows")
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
    out_dir = Path(args.out_dir)
    stats = FetchStats()
    seen_pairs: set = set()

    if args.source == "fms":
        mods = selected_mods(DEFAULT_MODS, args.mods)
        total = 0
        with TsvWriter(
            out_dir / "disease-alleles.tsv",
            COLUMNS,
            release=release,
            source=f"FMS DAF ({','.join(mods)})",
        ) as writer:
            for row in _fms_rows(mods):
                key = (row["alleleId"], row["diseaseId"])
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    break
        log.info("FMS Done. rows=%d", total)
        return

    if args.ids:
        disease_ids = [d.strip() for d in args.ids.split(",") if d.strip()]
    else:
        disease_ids = enumerate_disease_ids()
    if args.limit:
        disease_ids = disease_ids[: args.limit]
    total_rows = 0

    def work(did, cache, stats):
        results = paginate(f"/disease/{did}/alleles", cache=cache, stats=stats)
        return [row for row in (_row(did, r) for r in results) if row is not None]

    def consume(_did, rows):
        nonlocal total_rows
        for row in rows:
            pair_key = (row["alleleId"], row["diseaseId"])
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            writer.write_row(**row)
            total_rows += 1

    with TsvWriter(
        out_dir / "disease-alleles.tsv",
        COLUMNS,
        release=release,
        source="/disease/*/alleles",
    ) as writer:
        report = run_parallel(disease_ids, work, consume, cache_name="disease_alleles",
                              stats=stats, workers=args.workers, label="diseases")
        finish_run(report, stats, out_dir, args.max_failure_rate,
                   rows_written=writer.rows_written)

    log.info("disease_alleles: rows=%d", total_rows)

if __name__ == "__main__":
    main()
