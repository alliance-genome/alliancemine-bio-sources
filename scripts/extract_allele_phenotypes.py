#!/usr/bin/env python3
"""
Build allele-phenotypes.tsv (same columns as fetch_allele_phenotypes.py) from
a restored Alliance curation database, with the FMS PHENOTYPE files kept for
the MODs whose curation annotations carry no phenotype terms.

Two annotation kinds name an allele:
    AllelePhenotypeAnnotation  allele = subject, geneId = inferred gene, else
                               the allele's only gene (FMS: gene rows of SGD
                               listing alleles)
    AGMPhenotypeAnnotation     allele = inferred allele, geneId empty
                               (FMS: MGI rows keyed on the allele)
One row per (allele, phenotype term); the first row of a pair wins.

The 9.2.0 curation DB has no phenotypeannotation_ontologyterm rows for FB,
WB or ZFIN annotations (statements only), and termId is half of the
integration key, so those MODs still come from FMS through the fetcher's
_iter_mod (FMS_MODS below).

evidencePmid is the reference's PMID, else its MOD publication id.

Usage:
    python3 extract_allele_phenotypes.py
    python3 extract_allele_phenotypes.py --mods MGI,SGD --limit 100

Output: data/allele-phenotypes.tsv
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
    TsvWriter,
    configure_logging,
    enumerate_allele_ids,
    get_current_release,
    mod_gene_taxa,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402
from extract_phenotypes import reference_publication, reference_publications  # noqa: E402
from fetch_allele_phenotypes import COLUMNS, DEFAULT_MODS  # noqa: E402
from fetch_allele_phenotypes import _iter_mod as iter_fms_mod  # noqa: E402

log = logging.getLogger("alliance.extract.allele_phenotypes")

# Curation data-provider abbreviation -> MOD for providers that are not
# already MOD names. Xenbase covers two taxa.
XENBASE_BY_TAXON = {MOD_TAXON["XBXL"]: "XBXL", MOD_TAXON["XBXT"]: "XBXT"}

# MODs without phenotype term ids in the curation DB; read from FMS instead.
FMS_MODS = ("ZFIN", "FB", "WB")

ALLELE_PHENOTYPES_SQL = """
WITH sole_gene AS (
    SELECT alleleassociationsubject_id AS allele_id, min(allelegeneassociationobject_id) AS gene_id
    FROM allelegeneassociation WHERE NOT obsolete
    GROUP BY 1 HAVING count(DISTINCT allelegeneassociationobject_id) = 1
), ann AS (
    SELECT a.id, a.phenotypeannotationsubject_id AS allele_id,
           coalesce(a.inferredgene_id, sg.gene_id) AS gene_id
    FROM allelephenotypeannotation a
    LEFT JOIN sole_gene sg ON sg.allele_id = a.phenotypeannotationsubject_id
    UNION ALL
    SELECT id, inferredallele_id, NULL::bigint
    FROM agmphenotypeannotation WHERE inferredallele_id IS NOT NULL
)
SELECT org.abbreviation, atx.curie, al.primaryexternalid, coalesce(g.primaryexternalid, ''),
       t.curie, pa.phenotypeannotationobject, coalesce(ref.curie, '')
FROM ann
JOIN phenotypeannotation pa ON pa.id = ann.id
JOIN biologicalentity al ON al.id = ann.allele_id
LEFT JOIN biologicalentity g ON g.id = ann.gene_id AND NOT g.obsolete AND NOT g.internal
LEFT JOIN ontologyterm atx ON atx.id = al.taxon_id
LEFT JOIN organization org ON org.id = pa.dataprovider_id
JOIN phenotypeannotation_ontologyterm pot ON pot.phenotypeannotation_id = pa.id
JOIN ontologyterm t ON t.id = pot.phenotypeterms_id
LEFT JOIN informationcontententity ref ON ref.id = pa.evidenceitem_id
WHERE NOT pa.obsolete AND NOT pa.internal
  AND NOT al.obsolete AND NOT al.internal
  AND org.abbreviation = ANY(ARRAY[{providers}])
ORDER BY array_position(ARRAY[{providers}], org.abbreviation::text), atx.curie, pa.id, t.curie
"""


def mod_for(provider: str, taxon: str) -> str:
    """MOD name of a curation data provider; Xenbase splits by taxon."""
    if provider == "XB":
        return XENBASE_BY_TAXON.get(taxon, "")
    return provider


def evidence_publication(pubs: list[str], mod: str) -> str:
    """The FMS publicationId: PMID, else the publication id at the MOD."""
    return reference_publication(pubs, "Xenbase" if mod.startswith("XB") else mod)


def providers_sql(mods: list[str]) -> str:
    """Quoted curation provider names of the MODs, in MOD order."""
    names = dict.fromkeys("XB" if m.startswith("XB") else m for m in mods)
    return ",".join(f"'{n}'" for n in names)


def split_mods(mods: list[str]) -> tuple[list[str], list[str]]:
    """(curation DB MODs, FMS MODs), each in the given order."""
    return [m for m in mods if m not in FMS_MODS], [m for m in mods if m in FMS_MODS]


def iter_db_rows(mods: list[str], genes: set[str], ref_pubs: dict[str, list[str]]) -> Iterator[dict]:
    """Curation DB rows of the MODs, MOD by MOD in the given order."""
    if not mods:
        return
    sql = ALLELE_PHENOTYPES_SQL.format(providers=providers_sql(mods))
    for provider, taxon, allele, gene, term, statement, ref in copy_rows(sql):
        mod = mod_for(provider, taxon)
        if mod not in mods:
            continue
        yield {
            "geneId": gene if gene in genes else "",
            "alleleId": allele,
            "termId": term,
            "phenotypeStatement": statement,
            "evidencePmid": evidence_publication(ref_pubs.get(ref, []), mod),
            "dataProvider": mod,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance allele phenotypes from the curation DB")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    mods = selected_mods(DEFAULT_MODS, args.mods)
    db_mods, fms_mods = split_mods(mods)
    release = get_current_release()
    genes = set(mod_gene_taxa(mods))
    ref_pubs = reference_publications() if db_mods else {}
    alleles = set(enumerate_allele_ids(fms_mods)) if fms_mods else set()
    seen: set = set()

    def rows() -> Iterator[dict]:
        yield from iter_db_rows(db_mods, genes, ref_pubs)
        for mod in fms_mods:
            log.info("Allele phenotypes FMS: %s", mod)
            yield from iter_fms_mod(mod, genes, alleles)

    with TsvWriter(
        Path(args.out_dir) / "allele-phenotypes.tsv",
        COLUMNS,
        release=release,
        source=f"curation DB allele/AGM phenotype annotations ({','.join(db_mods)}); "
               f"FMS PHENOTYPE ({','.join(fms_mods)})",
    ) as writer:
        for row in rows():
            key = (row["alleleId"], row["termId"])
            if key in seen:
                continue
            seen.add(key)
            writer.write_row(**row)
            if args.limit is not None and writer.rows_written >= args.limit:
                break
        if writer.rows_written == 0:
            raise SystemExit("allele phenotypes: no rows extracted")

    log.info("Done. rows=%d", writer.rows_written)


if __name__ == "__main__":
    main()
