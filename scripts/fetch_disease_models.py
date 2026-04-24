#!/usr/bin/env python3
"""
Fetch disease-model associations from the Alliance API.

This data is MOD-curated (mouse, rat, zebrafish, fly, worm, xenopus, human)
and has no FMS equivalent. Yeast returns empty — the endpoint is populated
for mammalian/vertebrate models. Each hit links a gene to an "Affected
Genomic Model" (AGM, an allele-background strain/genotype) and to one or
more DO terms with is_model_of / is_marker_for associations.

Usage:
    python3 fetch_disease_models.py                         # all non-yeast MODs
    python3 fetch_disease_models.py --mods MGI,RGD          # subset
    python3 fetch_disease_models.py --ids MGI:98834         # smoke test
    python3 fetch_disease_models.py --limit 200

Output: data/disease-models.tsv
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
    configure_logging,
    enumerate_mod_genes,
    get_current_release,
    get_nested,
    join_pipe,
    open_cache,
    paginate,
)

log = logging.getLogger("alliance.fetch.disease_models")


# MODs where /gene/{id}/models returns data. SGD excluded (always empty).
DEFAULT_MODS = ["MGI", "RGD", "ZFIN", "FB", "WB", "XBXL", "XBXT", "HUMAN"]


COLUMNS = [
    "subjectGeneId",
    "subjectGeneTaxon",
    "modelId",
    "modelName",
    "modelSubtype",
    "dataProvider",
    "diseaseId",
    "diseaseName",
    "associationType",
    "associatedPhenotypes",
    "modifierRelationshipTypes",
    "hasDiseaseAnnotations",
    "hasPhenotypeAnnotations",
]


def _rows_for_result(seed_id: str, result: dict) -> list[dict]:
    """One API result wraps a single model but can carry multiple disease
    linkages via `diseaseModels`. Emit one row per (model, disease) pair."""
    gene = result.get("gene") or {}
    model = result.get("model") or {}
    gene_id = get_nested(gene, "primaryExternalId", default=seed_id) or seed_id
    gene_taxon = get_nested(gene, "taxon", "curie", default="")
    model_id = get_nested(model, "primaryExternalId", default="")
    model_name = get_nested(model, "agmFullName", "displayText", default="")
    model_subtype = get_nested(model, "subtype", "name", default="")
    dp = result.get("dataProvider", "") or get_nested(model, "dataProvider", "abbreviation", default="")

    phenotypes = result.get("associatedPhenotype") or []
    modifiers = result.get("modifierRelationshipTypes") or []
    has_disease = "true" if result.get("hasDiseaseAnnotations") else "false"
    has_pheno = "true" if result.get("hasPhenotypeAnnotations") else "false"

    disease_models = result.get("diseaseModels") or []
    rows: list[dict] = []
    if not disease_models:
        # Occasionally a model appears with phenotype-only annotation.
        rows.append({
            "subjectGeneId": gene_id,
            "subjectGeneTaxon": gene_taxon,
            "modelId": model_id,
            "modelName": model_name,
            "modelSubtype": model_subtype,
            "dataProvider": dp,
            "diseaseId": "",
            "diseaseName": "",
            "associationType": "",
            "associatedPhenotypes": join_pipe(phenotypes),
            "modifierRelationshipTypes": join_pipe(modifiers),
            "hasDiseaseAnnotations": has_disease,
            "hasPhenotypeAnnotations": has_pheno,
        })
        return rows

    for dm in disease_models:
        disease = dm.get("disease") or {}
        rows.append({
            "subjectGeneId": gene_id,
            "subjectGeneTaxon": gene_taxon,
            "modelId": model_id,
            "modelName": model_name,
            "modelSubtype": model_subtype,
            "dataProvider": dp,
            "diseaseId": disease.get("curie", ""),
            "diseaseName": disease.get("name", "") or dm.get("diseaseModel", ""),
            "associationType": dm.get("associationType", ""),
            "associatedPhenotypes": join_pipe(phenotypes),
            "modifierRelationshipTypes": join_pipe(modifiers),
            "hasDiseaseAnnotations": has_disease,
            "hasPhenotypeAnnotations": has_pheno,
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance disease-model associations")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--mods", help=f"Comma-separated MOD list (default: {','.join(DEFAULT_MODS)})")
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    if args.ids:
        gene_ids = [g.strip() for g in args.ids.split(",") if g.strip()]
    else:
        mods = args.mods.split(",") if args.mods else DEFAULT_MODS
        gene_ids = enumerate_mod_genes([m.strip() for m in mods])
    if args.limit:
        gene_ids = gene_ids[: args.limit]

    log.info("Processing %d genes", len(gene_ids))
    release = get_current_release()
    out_dir = Path(args.out_dir)

    try:
        from tqdm import tqdm
        iterator = tqdm(gene_ids, desc="genes", unit="gene")
    except ImportError:
        iterator = gene_ids

    stats = FetchStats()
    seen_keys: set = set()

    with open_cache("disease_models") as cache, \
         TsvWriter(
             out_dir / "disease-models.tsv",
             COLUMNS,
             release=release,
             source="/gene/*/models",
         ) as writer:
        total_rows = 0
        skipped_empty = 0
        for gid in iterator:
            try:
                page_had_results = False
                for result in paginate(f"/gene/{gid}/models", cache=cache, stats=stats):
                    page_had_results = True
                    for row in _rows_for_result(gid, result):
                        key = (row["modelId"], row["diseaseId"], row["subjectGeneId"])
                        if key in seen_keys:
                            continue
                        seen_keys.add(key)
                        writer.write_row(**row)
                        total_rows += 1
                if not page_had_results:
                    skipped_empty += 1
            except Exception as e:
                log.warning("Disease-model fetch failed for %s: %s", gid, e)

    log.info("Done. rows=%d  genes_without_models=%d  API stats: %s",
             total_rows, skipped_empty, stats)


if __name__ == "__main__":
    main()
