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
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    get_nested,
    join_pipe,
    mod_gene_taxa,
    paginate,
    resolve_gene_ids,
    run_parallel,
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


NOT_MODEL_PREFIX = "does not model "
MAX_UNMAPPED_RATE = 0.01


def load_do_names(obo_path: Path) -> dict[str, str]:
    """Lower-cased DO term name and EXACT synonym -> DOID. The 9.2.0
    /gene/{id}/models payload carries the disease name only (disease: {})."""
    names: dict[str, str] = {}
    synonyms: dict[str, str] = {}
    term_id = None
    with obo_path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                term_id = None
            elif line.startswith("id: DOID:"):
                term_id = line[4:]
            elif term_id and line.startswith("name: "):
                names.setdefault(line[6:].lower(), term_id)
            elif term_id and line.startswith('synonym: "') and " EXACT " in line:
                synonyms.setdefault(line.split('"')[1].lower(), term_id)
    if not names:
        raise RuntimeError(f"No DOID terms read from {obo_path}")
    for synonym, doid in synonyms.items():
        names.setdefault(synonym, doid)
    return names


def dedupe_key(row: dict) -> tuple:
    """Rows of one model differ by gene, disease and association type. An
    unmapped disease falls back to its name so two of them do not collapse."""
    return (row["modelId"], row["diseaseId"] or row["diseaseName"], row["subjectGeneId"],
            row["associationType"])


def unmapped_rate_exceeded(unmapped: int, named: int, max_rate: float = MAX_UNMAPPED_RATE) -> bool:
    return named > 0 and unmapped / named > max_rate


def disease_curie(dm: dict, do_names: dict[str, str]) -> str:
    curie = get_nested(dm, "disease", "curie", default="")
    if curie:
        return curie
    name = (get_nested(dm, "disease", "name", default="") or dm.get("diseaseModel", "")).lower()
    if name.startswith(NOT_MODEL_PREFIX):
        name = name[len(NOT_MODEL_PREFIX):]
    return do_names.get(name, "")


def _rows_for_result(seed_id: str, result: dict, gene_taxa: dict[str, str],
                     do_names: dict[str, str]) -> list[dict]:
    """One API result wraps a single model but can carry multiple disease
    linkages via `diseaseModels`. Emit one row per (model, disease) pair."""
    gene = result.get("gene") or {}
    model = result.get("model") or {}
    gene_id = get_nested(gene, "primaryExternalId", default=seed_id) or seed_id
    gene_taxon = (get_nested(gene, "taxon", "curie", default="")
                  or gene_taxa.get(gene_id) or gene_taxa.get(seed_id, ""))
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
            "diseaseId": disease_curie(dm, do_names),
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
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MODs to seed from (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--do-obo",
        help="DO ontology OBO for disease name lookup (default: <out-dir>/../fms/ONTOLOGY_DOID.obo)",
    )
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
    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit, default_mods=DEFAULT_MODS)
    out_dir = Path(args.out_dir)
    do_names = load_do_names(Path(args.do_obo) if args.do_obo else out_dir.parent / "fms" / "ONTOLOGY_DOID.obo")
    gene_taxa = mod_gene_taxa()
    stats = FetchStats()
    seen_keys: set = set()
    total_rows = 0
    genes_without_models = 0
    named_rows = 0
    unmapped_rows = 0

    def work(gid, cache, stats):
        return list(paginate(f"/gene/{gid}/models", cache=cache, stats=stats))

    def consume(gid, results):
        nonlocal total_rows, genes_without_models, named_rows, unmapped_rows
        if not results:
            genes_without_models += 1
        for result in results:
            for row in _rows_for_result(gid, result, gene_taxa, do_names):
                key = dedupe_key(row)
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                if row["diseaseName"]:
                    named_rows += 1
                    if not row["diseaseId"]:
                        unmapped_rows += 1
                writer.write_row(**row)
                total_rows += 1

    with TsvWriter(
        out_dir / "disease-models.tsv",
        COLUMNS,
        release=release,
        source="/gene/*/models",
    ) as writer:
        report = run_parallel(gene_ids, work, consume, cache_name="disease_models", stats=stats,
                              workers=args.workers, label="genes")
        log.info("disease_models: rows=%d  genes_without_models=%d  disease_names_unmapped=%d/%d",
                 total_rows, genes_without_models, unmapped_rows, named_rows)
        unmapped_exceeded = unmapped_rate_exceeded(unmapped_rows, named_rows)
        if unmapped_exceeded:
            log.error("disease_models: %d of %d disease names have no DOID (limit %.0f%%)",
                      unmapped_rows, named_rows, MAX_UNMAPPED_RATE * 100)
        try:
            finish_run(report, stats, out_dir, args.max_failure_rate,
                       rows_written=writer.rows_written)
        finally:
            if unmapped_exceeded:
                sys.exit(1)


if __name__ == "__main__":
    main()
