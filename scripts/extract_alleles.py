#!/usr/bin/env python3
"""
Build alliance-alleles.tsv (the 27-column VARIANT-ALLELE shape, see
fetch_alleles.py) from a restored Alliance curation database.

One row per allele and live variant, one row for an allele without variants,
as the VARIANT-ALLELE generator writes them. Columns follow the allele summary
document (AlleleDAO) and the generator:
    AlleleSymbol / AlleleSynonyms   symbol displayText / sorted synonym displayTexts
    AlleleAssociatedGene*           first live is_allele_of gene
    VariantAffectedGene*            first gene of the variant's gene-level predicted consequences
    Category                        "allele" / "allele with N known variant(s)"
    MostSevereConsequenceName       per gene-level predicted consequence, its VEP
                                    terms comma-joined, most severe first
    HasDisease/PhenotypeAnnotations "yes" or "-" (AlleleDAO existence queries)
    VariantInformationReference     one curie per variant reference: PMID, else
                                    the allele MOD's paper id, else any paper id

Scope: the AlleleIds of the release's FMS VARIANT-ALLELE files (the
fetcher's input) that the curation DB still serves; --all-alleles takes every
summary allele of the selected MODs.

Usage:
    python3 extract_alleles.py
    python3 extract_alleles.py --mods MGI --limit 1000

Output: data/alliance-alleles.tsv
"""

from __future__ import annotations

import argparse
import collections
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import TsvWriter, configure_logging, get_current_release, normalize_taxon, selected_mods  # noqa: E402
from curation_db import copy_rows  # noqa: E402
from extract_alleles_common import (  # noqa: E402
    VARIANT_ALLELE_MODS,
    allele_variant_links,
    first_genes,
    fms_category,
    information_references,
    load_variants,
    mod_for,
    seed_allele_ids,
    summary_alleles,
    yes_dash,
)
from fetch_alleles import COLUMNS  # noqa: E402

log = logging.getLogger("alliance.extract.alleles")


SYNONYMS_SQL = """
    SELECT singleallele_id, displaytext FROM slotannotation
    WHERE slotannotationtype = 'AlleleSynonymSlotAnnotation' AND displaytext IS NOT NULL
"""

# AlleleDAO phenotype / disease existence queries.
PHENOTYPE_SQL = """
    SELECT phenotypeannotationsubject_id FROM allelephenotypeannotation
    UNION SELECT inferredallele_id FROM agmphenotypeannotation
    UNION SELECT assertedalleles_id FROM agmphenotypeannotation_allele
"""
DISEASE_SQL = """
    SELECT ada.diseaseannotationsubject_id
    FROM allelediseaseannotation ada JOIN diseaseannotation da ON da.id = ada.id
    WHERE da.negated = false
    UNION SELECT agmda.inferredallele_id
    FROM agmdiseaseannotation agmda JOIN diseaseannotation da ON da.id = agmda.id
    WHERE da.negated = false
    UNION SELECT agmdaa.assertedalleles_id
    FROM agmdiseaseannotation_allele agmdaa
    JOIN diseaseannotation da ON da.id = agmdaa.agmdiseaseannotation_id
    WHERE da.negated = false
"""
SPECIES_SQL = """
    SELECT ot.curie, org.abbreviation, sp.fullname
    FROM species sp
    LEFT JOIN ontologyterm ot ON ot.id = sp.taxon_id
    LEFT JOIN organization org ON org.id = sp.dataprovider_id
    ORDER BY sp.id
"""


def species_names(rows) -> tuple[dict[str, str], dict[str, str]]:
    """(taxon curie -> species name, provider -> first species name). SGD
    strain taxa have no species row and resolve through the provider."""
    by_taxon: dict[str, str] = {}
    by_provider: dict[str, str] = {}
    for taxon, provider, name in rows:
        by_taxon.setdefault(taxon, name)
        by_provider.setdefault(provider, name)
    return by_taxon, by_provider


def allele_fields(pid, taxon, species, symbol, synonyms, gene_id, gene_symbol,
                  has_disease: bool, has_phenotype: bool, variant_count: int) -> dict:
    return {
        "Taxon": normalize_taxon(taxon),
        "SpeciesName": species,
        "AlleleId": pid,
        "AlleleSymbol": symbol,
        "AlleleSynonyms": "|".join(sorted(synonyms)),
        "AlleleAssociatedGeneId": gene_id,
        "AlleleAssociatedGeneSymbol": gene_symbol,
        "Category": fms_category(variant_count),
        "HasDiseaseAnnotations": yes_dash(has_disease),
        "HasPhenotypeAnnotations": yes_dash(has_phenotype),
    }


def variant_fields(v: dict, mod_prefix: str) -> dict:
    return {
        "VariantId": v["hgvs"],
        "VariantSymbol": v["symbol"],
        "VariantSynonyms": v["hgvs"],
        "VariantCrossReferences": "|".join(v["xrefs"]),
        "VariantAffectedGeneId": v["gene_ids"][0] if v["gene_ids"] else "",
        "VariantAffectedGeneSymbol": v["gene_symbols"][0] if v["gene_symbols"] else "",
        "VariantsTypeId": v["type_id"],
        "VariantsTypeName": v["type_name"],
        "VariantsHgvsNames": v["hgvs"],
        "Assembly": v["assembly"],
        "Chromosome": v["chromosome"],
        "StartPosition": v["start"],
        "EndPosition": v["end"],
        "SequenceOfReference": v["ref"],
        "SequenceOfVariant": v["var"],
        "MostSevereConsequenceName": "|".join(v["consequences"]),
        "VariantInformationReference": "|".join(information_references(v["references"], mod_prefix)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance alleles into the legacy 27-col TSV")
    parser.add_argument("--mods", help="Comma-separated MODs (default: all, env ALLIANCE_MODS)")
    parser.add_argument("--all-alleles", action="store_true",
                        help="Every summary allele, not only the FMS VARIANT-ALLELE seed")
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    release = get_current_release()
    mods = [m for m in selected_mods(None, args.mods) if m in VARIANT_ALLELE_MODS]
    wanted = seed_allele_ids(mods, args.all_alleles)
    alleles = summary_alleles(copy_rows, mods)
    if wanted is not None:
        alleles = [a for a in alleles if a[1] in wanted]
    log.info("%d alleles in scope", len(alleles))
    in_scope = {a[0] for a in alleles}
    synonyms: dict[str, list[str]] = collections.defaultdict(list)
    for allele_id, name in copy_rows(SYNONYMS_SQL):
        if allele_id in in_scope:
            synonyms[allele_id].append(name)
    genes = first_genes(copy_rows)
    with_phenotype = {r[0] for r in copy_rows(PHENOTYPE_SQL) if r}
    with_disease = {r[0] for r in copy_rows(DISEASE_SQL) if r}
    by_taxon, by_provider = species_names(copy_rows(SPECIES_SQL))
    variants = load_variants(copy_rows)
    links = allele_variant_links(copy_rows)

    per_mod: collections.Counter = collections.Counter()
    seen: set[str] = set()
    with TsvWriter(
        Path(args.out_dir) / "alliance-alleles.tsv",
        COLUMNS,
        release=release,
        source=f"curation DB allele summary -> legacy 27-col TSV ({','.join(mods)})",
    ) as writer:
        for db_id, pid, provider, taxon, symbol, _fmt, _xref in alleles:
            if pid in seen:
                continue
            seen.add(pid)
            gene_id, gene_symbol = genes.get(db_id, ("", ""))
            allele_variants = [variants[v] for v in links.get(db_id, []) if v in variants]
            base = allele_fields(pid, taxon, by_taxon.get(taxon) or by_provider.get(provider, ""),
                                 symbol, synonyms.get(db_id, []), gene_id, gene_symbol,
                                 db_id in with_disease, db_id in with_phenotype, len(allele_variants))
            for v in allele_variants or [None]:
                writer.write_row(**base, **(variant_fields(v, provider) if v else {}))
                per_mod[mod_for(provider, taxon)] += 1
                if args.limit and writer.rows_written >= args.limit:
                    break
            if args.limit and writer.rows_written >= args.limit:
                break
    if wanted is not None and not args.limit:
        log.info("%d seed alleles not served (absent, obsolete or internal in the DB)", len(wanted - seen))
    log.info("alliance-alleles rows per MOD: %s", dict(sorted(per_mod.items())))


if __name__ == "__main__":
    main()
