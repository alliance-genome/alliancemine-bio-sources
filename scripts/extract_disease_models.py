#!/usr/bin/env python3
"""
Build disease-models.tsv (fetch_disease_models.py's TSV) from a restored
curation database.

/gene/{id}/models lists the affected genomic models (AGMs) of a gene: AGMs
holding an allele associated with the gene (any allele-gene relation),
AGMs holding a sequence targeting reagent that targets it, and AGMs whose
disease or phenotype annotations infer or assert the gene. Construct-only
links are not models. Obsolete or internal AGMs, alleles and associations
are left out. This reproduces the API's model set for the 28 genes checked.

Each (gene, model) gives one row per disease the model's annotations name
(associationType IS_MODEL_OF, IS_NOT_MODEL_OF when negated), or one row
without disease when it has none. associatedPhenotypes are the model's
phenotype statements, sorted; modifierRelationshipTypes is the fixed list the
API reports.

Output: data/disease-models.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    get_current_release,
    join_pipe,
    resolve_gene_ids,
)
from curation_db import copy_rows  # noqa: E402
from disease_db import split_pipe  # noqa: E402
from fetch_disease_models import COLUMNS, DEFAULT_MODS, dedupe_key  # noqa: E402

log = logging.getLogger("alliance.extract.disease_models")

MODIFIER_RELATIONSHIP_TYPES = "ameliorated_by|exacerbated_by"
NOT_MODEL_PREFIX = "does not model "

GENE_MODELS_SQL = """
    WITH live AS (SELECT id FROM biologicalentity WHERE NOT obsolete AND NOT internal),
    links AS (
        SELECT ag.allelegeneassociationobject_id AS gene, aa.agmassociationsubject_id AS agm
        FROM allelegeneassociation ag
        JOIN agmalleleassociation aa ON aa.agmalleleassociationobject_id = ag.alleleassociationsubject_id
             AND NOT aa.obsolete AND NOT aa.internal
        WHERE NOT ag.obsolete AND NOT ag.internal AND ag.alleleassociationsubject_id IN (SELECT id FROM live)
        UNION
        SELECT sg.sequencetargetingreagentgeneassociationobject_id, sa.agmassociationsubject_id
        FROM sqtrgeneassociation sg
        JOIN agmsequencetargetingreagentassociation sa
             ON sa.agmsequencetargetingreagentassociationobject_id = sg.sequencetargetingreagentassociationsubject_id
             AND NOT sa.obsolete AND NOT sa.internal
        WHERE NOT sg.obsolete AND NOT sg.internal
        UNION
        SELECT m.inferredgene_id, m.diseaseannotationsubject_id FROM agmdiseaseannotation m
        JOIN diseaseannotation d ON d.id = m.id AND NOT d.obsolete AND NOT d.internal
        UNION
        SELECT x.assertedgenes_id, m.diseaseannotationsubject_id FROM agmdiseaseannotation_gene x
        JOIN agmdiseaseannotation m ON m.id = x.agmdiseaseannotation_id
        JOIN diseaseannotation d ON d.id = m.id AND NOT d.obsolete AND NOT d.internal
        UNION
        SELECT m.inferredgene_id, m.phenotypeannotationsubject_id FROM agmphenotypeannotation m
        JOIN phenotypeannotation p ON p.id = m.id AND NOT p.obsolete AND NOT p.internal
        UNION
        SELECT x.assertedgenes_id, m.phenotypeannotationsubject_id FROM agmphenotypeannotation_gene x
        JOIN agmphenotypeannotation m ON m.id = x.agmphenotypeannotation_id
        JOIN phenotypeannotation p ON p.id = m.id AND NOT p.obsolete AND NOT p.internal
    )
    SELECT g.primaryexternalid, gt.curie, a.primaryexternalid
    FROM links
    JOIN biologicalentity g ON g.id = links.gene AND NOT g.obsolete AND NOT g.internal
    JOIN gene ON gene.id = g.id
    JOIN biologicalentity a ON a.id = links.agm AND NOT a.obsolete AND NOT a.internal
    JOIN affectedgenomicmodel m ON m.id = a.id
    LEFT JOIN ontologyterm gt ON gt.id = g.taxon_id
"""

MODELS_SQL = """
    WITH dis AS (
        SELECT m.diseaseannotationsubject_id AS agm,
               string_agg(DISTINCT d.negated::text || '=' || o.curie || '=' || o.name, '|') AS diseases
        FROM agmdiseaseannotation m
        JOIN diseaseannotation d ON d.id = m.id AND NOT d.obsolete AND NOT d.internal
        JOIN ontologyterm o ON o.id = d.diseaseannotationobject_id
        GROUP BY 1
    ),
    phe AS (
        SELECT m.phenotypeannotationsubject_id AS agm,
               string_agg(DISTINCT p.phenotypeannotationobject, '|') AS phenotypes, count(*) AS n
        FROM agmphenotypeannotation m
        JOIN phenotypeannotation p ON p.id = m.id AND NOT p.obsolete AND NOT p.internal
        GROUP BY 1
    ),
    name AS (
        SELECT s.singleagm_id AS agm, min(s.displaytext) AS name FROM slotannotation s
        WHERE s.slotannotationtype = 'AgmFullNameSlotAnnotation' AND NOT s.obsolete
        GROUP BY 1
    )
    SELECT b.primaryexternalid, name.name, st.name, org.abbreviation, dis.diseases, phe.phenotypes,
           (phe.n IS NOT NULL)::text
    FROM affectedgenomicmodel m
    JOIN biologicalentity b ON b.id = m.id AND NOT b.obsolete AND NOT b.internal
    LEFT JOIN vocabularyterm st ON st.id = m.subtype_id
    LEFT JOIN organization org ON org.id = b.dataprovider_id
    LEFT JOIN name ON name.agm = m.id
    LEFT JOIN dis ON dis.agm = m.id
    LEFT JOIN phe ON phe.agm = m.id
"""


def disease_models(diseases: str) -> list[dict]:
    """'false=DOID:1=name|true=DOID:2=name' -> the API's diseaseModels, sorted by name."""
    out = []
    for item in split_pipe(diseases):
        negated, doid, name = item.split("=", 2)
        negated_flag = negated == "true"
        out.append({
            "diseaseId": doid,
            "diseaseName": (NOT_MODEL_PREFIX + name) if negated_flag else name,
            "associationType": "IS_NOT_MODEL_OF" if negated_flag else "IS_MODEL_OF",
        })
    out.sort(key=lambda d: (d["diseaseName"].lower(), d["associationType"]))
    return out


def sorted_phenotypes(phenotypes: str) -> list[str]:
    return sorted(split_pipe(phenotypes), key=str.lower)


def model_rows(gene_id: str, gene_taxon: str, model: dict) -> list[dict]:
    base = {
        "subjectGeneId": gene_id,
        "subjectGeneTaxon": gene_taxon,
        "modelId": model["id"],
        "modelName": model["name"],
        "modelSubtype": model["subtype"],
        "dataProvider": model["provider"],
        "associatedPhenotypes": join_pipe(model["phenotypes"]),
        "modifierRelationshipTypes": MODIFIER_RELATIONSHIP_TYPES,
        "hasDiseaseAnnotations": "true" if model["diseases"] else "false",
        "hasPhenotypeAnnotations": "true" if model["has_phenotypes"] else "false",
    }
    if not model["diseases"]:
        return [dict(base, diseaseId="", diseaseName="", associationType="")]
    return [dict(base, **dm) for dm in model["diseases"]]


def load_models() -> dict[str, dict]:
    models = {}
    for agm_id, name, subtype, provider, diseases, phenotypes, has_pheno in copy_rows(MODELS_SQL):
        models[agm_id] = {
            "id": agm_id, "name": name, "subtype": subtype, "provider": provider,
            "diseases": disease_models(diseases), "phenotypes": sorted_phenotypes(phenotypes),
            "has_phenotypes": has_pheno == "true",
        }
    log.info("models: %d", len(models))
    return models


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance disease-model associations from the curation DB")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MODs (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    release = get_current_release()
    # The fetcher's seeds: --ids, else the release's gene list of the MODs.
    seeds = resolve_gene_ids(args.ids, args.mods, args.limit, default_mods=DEFAULT_MODS)
    wanted = set(seeds)

    gene_models: dict[str, list[str]] = {}
    gene_taxa: dict[str, str] = {}
    for gene_id, taxon, agm_id in copy_rows(GENE_MODELS_SQL):
        if gene_id in wanted:
            gene_models.setdefault(gene_id, []).append(agm_id)
            gene_taxa[gene_id] = taxon
    genes = [g for g in seeds if g in gene_models]
    models = load_models()

    seen: set = set()
    with TsvWriter(Path(args.out_dir) / "disease-models.tsv", COLUMNS, release=release,
                   source="curation DB (AGM components and annotations)") as writer:
        for gene_id in genes:
            for agm_id in sorted(set(gene_models[gene_id])):
                model = models.get(agm_id)
                if model is None:
                    continue
                for row in model_rows(gene_id, gene_taxa[gene_id], model):
                    key = dedupe_key(row)
                    if key in seen:
                        continue
                    seen.add(key)
                    writer.write_row(**row)
    log.info("disease_models: rows=%d genes_with_models=%d of %d", writer.rows_written, len(genes), len(seeds))


if __name__ == "__main__":
    main()
