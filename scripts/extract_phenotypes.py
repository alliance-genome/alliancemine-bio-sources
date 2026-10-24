#!/usr/bin/env python3
"""
Build phenotypes.tsv (same columns as fetch_phenotypes.py) from a restored
Alliance curation database instead of the /gene/{id}/phenotypes API.

The API serves a per-gene rollup: one record per (gene, phenotype statement)
whose primary annotations are the gene's own GenePhenotypeAnnotations plus
every AllelePhenotypeAnnotation and AGMPhenotypeAnnotation whose inferred
gene is that gene. This script rebuilds the same grouping in SQL:

    relation    has_phenotype when a gene-level annotation is in the group,
                otherwise is_implicated_in (the API's choice)
    uniqueId    "<gene>|<statement>"; the API value is a random search-index
                document id, so no source can reproduce it
    pubmedIds   per reference its PMID, else its MOD publication id
    referenceIds AGRKB reference curies, sorted and distinct

PMIDs come from the curation reference_crossreference table; the script stops
when that table is missing from the restore.

Usage:
    python3 extract_phenotypes.py                    # every gene of the 9 MODs
    python3 extract_phenotypes.py --mods SGD,FB
    python3 extract_phenotypes.py --ids SGD:S000004103

Output: data/phenotypes.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    get_current_release,
    join_pipe,
    resolve_gene_ids,
)
from curation_db import copy_rows, scalar  # noqa: E402
from fetch_phenotypes import COLUMNS  # noqa: E402

log = logging.getLogger("alliance.extract.phenotypes")

CATEGORY = "gene_phenotype_annotation"
# Gene sets up to this size are pushed into the SQL as a VALUES join.
SQL_ID_FILTER_MAX = 5000

# One row per (gene, statement, reference) with a flag for gene-level
# annotations. {gene_filter} is empty or a join on the selected gene ids.
ROLLUP_SQL = """
WITH prim AS (
    SELECT phenotypeannotationsubject_id AS gene_id, id, true AS direct FROM genephenotypeannotation
    UNION ALL
    SELECT inferredgene_id, id, false FROM allelephenotypeannotation WHERE inferredgene_id IS NOT NULL
    UNION ALL
    SELECT inferredgene_id, id, false FROM agmphenotypeannotation WHERE inferredgene_id IS NOT NULL
)
SELECT g.primaryexternalid, tx.curie, pa.phenotypeannotationobject,
       bool_or(p.direct), coalesce(ref.curie, '')
FROM prim p
JOIN phenotypeannotation pa ON pa.id = p.id
JOIN biologicalentity g ON g.id = p.gene_id
{gene_filter}
LEFT JOIN ontologyterm tx ON tx.id = g.taxon_id
LEFT JOIN informationcontententity ref ON ref.id = pa.evidenceitem_id
WHERE NOT pa.obsolete AND NOT pa.internal
  AND NOT g.obsolete AND NOT g.internal
  AND coalesce(pa.phenotypeannotationobject, '') <> ''
GROUP BY 1, 2, 3, 5
"""

# Reference curie -> publication cross-reference curies (PMIDs and MOD ids).
REFERENCE_PUBS_SQL = """
SELECT ref.curie, x.referencedcurie
FROM reference_crossreference rx
JOIN informationcontententity ref ON ref.id = rx.reference_id
JOIN crossreference x ON x.id = rx.crossreferences_id
WHERE x.referencedcurie ~ '^(PMID|MGI|RGD|ZFIN|FB|WB|SGD|Xenbase):'
  AND NOT coalesce(x.obsolete, false)
"""


def relation_for(has_gene_annotation: bool) -> str:
    """The API names the rollup has_phenotype when the gene itself is the
    subject of one of its annotations, is_implicated_in when every annotation
    is on an allele or a model of the gene."""
    return "has_phenotype" if has_gene_annotation else "is_implicated_in"


def unique_id(gene_id: str, statement: str) -> str:
    return f"{gene_id}|{statement}"


def reference_publication(pubs: list[str], mod_prefix: str) -> str:
    """One publication id per reference, the API's referenceID rule: the PMID
    when there is one, else the id at the MOD (prefix such as 'MGI' or
    'Xenbase')."""
    pmids = sorted(p for p in pubs if p.startswith("PMID:"))
    if pmids:
        return pmids[0]
    mod_ids = sorted(p for p in pubs if p.startswith(f"{mod_prefix}:"))
    return mod_ids[0] if mod_ids else ""


def publication_ids(references: Iterable[str], ref_pubs: dict[str, list[str]], gene_id: str) -> list[str]:
    """Distinct sorted publication ids of the references, as the API's
    pubmedPublications list: a reference without a PMID contributes its id at
    the gene's MOD, and a non-AGRKB evidence item (MIM:, ORPHA: on human
    genes) contributes its own curie."""
    prefix = gene_id.split(":", 1)[0]
    pubs = (reference_publication(ref_pubs.get(r, []), prefix) or ("" if r.startswith("AGRKB:") else r)
            for r in references)
    return sorted({p for p in pubs if p})


def gene_filter_sql(gene_ids: list[str] | None) -> str:
    if not gene_ids or len(gene_ids) > SQL_ID_FILTER_MAX:
        return ""
    values = ",".join("('" + g.replace("'", "''") + "')" for g in gene_ids)
    return f"JOIN (VALUES {values}) AS sel(pid) ON sel.pid = g.primaryexternalid"


def reference_publications() -> dict[str, list[str]]:
    """AGRKB reference curie -> publication curies, from reference_crossreference."""
    if scalar("SELECT to_regclass('public.reference_crossreference') IS NOT NULL") != "t":
        raise SystemExit(
            "reference_crossreference is missing from the curation restore; PMIDs cannot be "
            "resolved. Restore it: pg_restore -t reference_crossreference <dump>"
        )
    pubs: dict[str, list[str]] = defaultdict(list)
    for ref, xref in copy_rows(REFERENCE_PUBS_SQL):
        pubs[ref].append(xref)
    log.info("References with publication ids: %d", len(pubs))
    return pubs


def iter_rollups(gene_ids: list[str]) -> Iterator[dict]:
    """Group the per-reference SQL rows into one dict per (gene, statement)."""
    wanted = set(gene_ids)
    groups: dict[tuple[str, str], dict] = {}
    for gene, taxon, statement, direct, ref in copy_rows(
            ROLLUP_SQL.format(gene_filter=gene_filter_sql(gene_ids))):
        if gene not in wanted:
            continue
        group = groups.setdefault((gene, statement), {
            "gene": gene, "taxon": taxon, "statement": statement, "direct": False, "refs": set(),
        })
        group["direct"] = group["direct"] or direct == "t"
        if ref:
            group["refs"].add(ref)
    order = {g: i for i, g in enumerate(gene_ids)}
    yield from sorted(groups.values(), key=lambda g: (order[g["gene"]], g["statement"]))


def to_row(group: dict, ref_pubs: dict[str, list[str]]) -> dict:
    refs = sorted(group["refs"])
    return {
        "subjectGeneId": group["gene"],
        "subjectGeneTaxon": group["taxon"],
        "category": CATEGORY,
        "phenotypeStatement": group["statement"],
        "relation": relation_for(group["direct"]),
        "uniqueId": unique_id(group["gene"], group["statement"]),
        "pubmedIds": join_pipe(publication_ids(refs, ref_pubs, group["gene"])),
        "referenceIds": join_pipe(refs),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance phenotype annotations from the curation DB")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs")
    parser.add_argument("--mods", help="Comma-separated MODs (default: all 9, env ALLIANCE_MODS)")
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
    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit)
    ref_pubs = reference_publications()
    out_dir = Path(args.out_dir)

    with TsvWriter(
        out_dir / "phenotypes.tsv",
        COLUMNS,
        release=release,
        source="curation DB gene/allele/AGM phenotype annotations",
    ) as writer:
        for group in iter_rollups(gene_ids):
            writer.write_row(**to_row(group, ref_pubs))
        if writer.rows_written == 0:
            raise SystemExit("phenotypes: no rows extracted")

    log.info("phenotypes: rows=%d", writer.rows_written)


if __name__ == "__main__":
    main()
