#!/usr/bin/env python3
"""
Build molecular-interactions.tsv and genetic-interactions.tsv from a restored
Alliance curation database, in the column shape fetch_interactions.py writes.

The fetcher asked /gene/{id}/*-interactions for every seeded gene and kept
each uniqueId once, oriented with the smaller gene id first when that gene
was seeded, else from the seeded partner's side. This script streams every
interaction from SQL and applies the same rule against the same seed list.

Publications: referenceId is the AGRKB curie of the evidence row and
pubmedId its PMID (else MOD reference id) from reference_crossreference. The
dump has no citation text, so shortCitation comes from the API's
/reference/{curie} (one cached call per distinct publication);
--no-reference-lookup skips it and leaves shortCitation empty.

Usage:
    python3 extract_interactions.py                  # every gene of the 9 MODs
    python3 extract_interactions.py --mods SGD,HUMAN
    python3 extract_interactions.py --ids SGD:S000004103,SGD:S000002429
    python3 extract_interactions.py --limit 50       # smoke test

Output: data/molecular-interactions.tsv, data/genetic-interactions.tsv
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import urllib.error
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    http_get_json,
    resolve_gene_ids,
    run_parallel,
)
from curation_db import copy_rows  # noqa: E402
from fetch_interactions import GENETIC_COLUMNS, MOLECULAR_COLUMNS  # noqa: E402

log = logging.getLogger("alliance.extract.interactions")

# Term columns shared by both tables: (TSV curie column, interaction FK column).
_TERM_FIELDS = [
    ("interactionType", "interactiontype_id"),
    ("interactorAType", "interactoratype_id"),
    ("interactorBType", "interactorbtype_id"),
    ("interactorARole", "interactorarole_id"),
    ("interactorBRole", "interactorbrole_id"),
    ("interactionSource", "interactionsource_id"),
]
_KINDS = {
    "molecular": {
        "table": "genemolecularinteraction",
        "terms": _TERM_FIELDS + [
            ("detectionMethod", "detectionmethod_id"),
            ("aggregationDatabase", "aggregationdatabase_id"),
        ],
        "columns": MOLECULAR_COLUMNS,
    },
    "genetic": {
        "table": "genegeneticinteraction",
        "terms": _TERM_FIELDS,
        "columns": GENETIC_COLUMNS,
    },
}
# SQL result columns before the term pairs; extra columns follow the terms.
_BASE_FIELDS = [
    "subjectId", "subjectSymbol", "subjectTaxon",
    "objectId", "objectSymbol", "objectTaxon",
    "relation", "interactionId", "uniqueId", "referenceId", "pubmedId", "crossRefs",
]
# A publication without a PubMed id is cited by its MOD reference id.
_INTERACTOR_A_COLUMNS = ["interactorAType", "interactorATypeName", "interactorARole", "interactorARoleName"]
_PUBLICATION_PREFIXES = "'PMID', 'FB', 'WB', 'MGI', 'ZFIN', 'RGD', 'SGD', 'Xenbase'"
_CITATION_HEAD = re.compile(r"^(?P<author>.*?)(?: et al\.)? \((?P<year>[^)]*)\) ?(?P<tail>.*)$")


def interactions_sql(kind: str) -> str:
    """One row per interaction of `kind`, participants and terms resolved."""
    spec = _KINDS[kind]
    table = spec["table"]
    term_ids = " UNION ".join(f"SELECT {fk} FROM {table}" for _, fk in spec["terms"])
    term_joins = "\n".join(
        f"LEFT JOIN terms t{n} ON t{n}.id = i.{fk}" for n, (_, fk) in enumerate(spec["terms"])
    )
    term_cols = ", ".join(f"t{n}.curie, t{n}.name" for n in range(len(spec["terms"])))
    extra_cte = extra_col = extra_join = ""
    if kind == "genetic":
        extra_cte = f""",
phenotypes AS (
    SELECT {table}_id AS id, string_agg(phenotypesortraits, '|' ORDER BY phenotypesortraits COLLATE "C") AS value
    FROM {table}_phenotypesortraits
    WHERE phenotypesortraits <> ''
    GROUP BY 1
)"""
        extra_col = ", ph.value"
        extra_join = "LEFT JOIN phenotypes ph ON ph.id = i.id"
    return f"""
WITH participants AS (
    SELECT geneassociationsubject_id AS id FROM {table}
    UNION SELECT genegeneassociationobject_id FROM {table}
),
genes AS (
    SELECT be.id, be.primaryexternalid, tx.curie AS taxon
    FROM participants p
    JOIN biologicalentity be ON be.id = p.id
    LEFT JOIN ontologyterm tx ON tx.id = be.taxon_id
    WHERE NOT be.obsolete AND NOT be.internal
),
symbols AS (
    SELECT sa.singlegene_id AS id, min(sa.displaytext) AS symbol
    FROM slotannotation sa
    JOIN participants p ON p.id = sa.singlegene_id
    WHERE sa.slotannotationtype = 'GeneSymbolSlotAnnotation' AND NOT sa.obsolete
    GROUP BY 1
),
terms AS (
    SELECT id, curie, name FROM ontologyterm WHERE id IN ({term_ids})
),
evidence AS (
    SELECT e.association_id AS id, min(ice.curie) AS curie, min(ice.id) AS reference_id
    FROM {table}_informationcontententity e
    JOIN informationcontententity ice ON ice.id = e.evidence_id
    GROUP BY 1
),
publication_ids AS (
    SELECT rx.reference_id AS id,
           (array_agg(x.referencedcurie
                      ORDER BY x.referencedcurie NOT LIKE 'PMID:%', x.referencedcurie))[1] AS curie
    FROM reference_crossreference rx
    JOIN crossreference x ON x.id = rx.crossreferences_id
    WHERE rx.reference_id IN (SELECT reference_id FROM evidence)
      AND NOT x.obsolete
      AND split_part(x.referencedcurie, ':', 1) IN ({_PUBLICATION_PREFIXES})
    GROUP BY 1
),
xrefs AS (
    SELECT l.geneinteraction_id AS id, string_agg(x.referencedcurie, '|' ORDER BY x.id) AS curies
    FROM {table}_crossreference l
    JOIN crossreference x ON x.id = l.crossreferences_id
    WHERE x.referencedcurie <> ''
    GROUP BY 1
){extra_cte}
SELECT sg.primaryexternalid, ss.symbol, sg.taxon,
       og.primaryexternalid, os.symbol, og.taxon,
       rel.name, i.interactionid, i.uniqueid, ev.curie, pub.curie, xr.curies,
       {term_cols}{extra_col}
FROM {table} i
JOIN genes sg ON sg.id = i.geneassociationsubject_id
JOIN genes og ON og.id = i.genegeneassociationobject_id
LEFT JOIN symbols ss ON ss.id = i.geneassociationsubject_id
LEFT JOIN symbols os ON os.id = i.genegeneassociationobject_id
LEFT JOIN vocabularyterm rel ON rel.id = i.relation_id
LEFT JOIN evidence ev ON ev.id = i.id
LEFT JOIN publication_ids pub ON pub.id = ev.reference_id
LEFT JOIN xrefs xr ON xr.id = i.id
{term_joins}
{extra_join}
WHERE NOT i.obsolete AND NOT i.internal
ORDER BY i.id
"""


def references_sql() -> str:
    """Distinct publication curies cited by either interaction table."""
    return """
SELECT DISTINCT ice.curie
FROM (SELECT evidence_id FROM genemolecularinteraction_informationcontententity
      UNION SELECT evidence_id FROM genegeneticinteraction_informationcontententity) e
JOIN informationcontententity ice ON ice.id = e.evidence_id
ORDER BY 1
"""


def queried_gene(subject_id: str, object_id: str, seeded: set[str]) -> str | None:
    """The gene whose API call the fetcher kept this interaction from.

    fetch_interactions.py writes one row per uniqueId, preferring the copy
    whose queried gene sorts first; the larger id wins only when the smaller
    one was not seeded. None when neither participant was seeded."""
    low, high = sorted((subject_id, object_id))
    if low in seeded:
        return low
    if high in seeded:
        return high
    return None


def sql_row_to_record(kind: str, values: list[str]) -> dict:
    """Map one interactions_sql() result row to named fields."""
    spec = _KINDS[kind]
    rec = dict(zip(_BASE_FIELDS, values))
    pos = len(_BASE_FIELDS)
    for column, _ in spec["terms"]:
        rec[column], rec[column + "Name"] = values[pos], values[pos + 1]
        pos += 2
    if kind == "genetic":
        rec["phenotypesOrTraits"] = values[pos]
    return rec


def tsv_row(rec: dict, gene_id: str, citation: str) -> dict:
    """Orient rec on gene_id like the API: when gene_id is the object, gene
    and partner swap and so do the interactor A and B type/role columns, so
    interactor A always describes geneId."""
    swap = rec["subjectId"] != gene_id
    side, other = ("object", "subject") if swap else ("subject", "object")
    row = {k: v for k, v in rec.items() if not k.startswith(("subject", "object"))}
    if swap:
        for a_col in _INTERACTOR_A_COLUMNS:
            b_col = a_col.replace("interactorA", "interactorB")
            row[a_col], row[b_col] = rec[b_col], rec[a_col]
    row.update({
        "geneId": rec[f"{side}Id"],
        "geneSymbol": rec[f"{side}Symbol"],
        "geneTaxon": rec[f"{side}Taxon"],
        "partnerGeneId": rec[f"{other}Id"],
        "partnerSymbol": rec[f"{other}Symbol"],
        "partnerTaxon": rec[f"{other}Taxon"],
        "shortCitation": citation,
    })
    return row


def short_citation(summary: dict) -> str:
    """'Author AB (year) Journal volume(issue):pages', every part present
    even when empty. The literature service's short_citation adds 'et al.'
    and drops empty parts; the journal is blank when it lists no resource.
    Input is an API /reference literatureSummary."""
    text = summary.get("short_citation") or ""
    m = _CITATION_HEAD.match(text)
    if not m:
        return text
    volume = summary.get("volume") or ""
    issue = summary.get("issue_name") or ""
    pages = summary.get("page_range") or ""
    journal = m["tail"]
    if not summary.get("resource_title"):
        journal = " "
    elif volume:
        head = re.match(rf"(.*?) {re.escape(volume)}(?:\(|:|$)", journal)
        journal = head.group(1) if head else journal
    elif pages and journal.endswith(pages):
        journal = journal[: -len(pages)].rstrip(" :")
    return f"{m['author']} ({m['year']}) {journal} {volume}({issue}):{pages}"


def fetch_citation(curie: str, cache, stats: FetchStats) -> str:
    """Short citation of one publication; '' when the API has no record
    (HTTP 400 for 317 of 75.5K curies in 9.2.0)."""
    try:
        data = http_get_json(f"/reference/{curie}", cache=cache, stats=stats)
    except urllib.error.HTTPError as e:
        if e.code in (400, 404):
            log.debug("no /reference record for %s (HTTP %s)", curie, e.code)
            return ""
        raise
    return short_citation(data.get("literatureSummary") or {})


def lookup_citations(curies: list[str], workers: int, stats: FetchStats):
    """Map each publication curie to its short citation via the API."""
    found: dict[str, str] = {}

    def consume(curie, result):
        found[curie] = result

    report = run_parallel(curies, fetch_citation, consume, cache_name="references",
                          stats=stats, workers=workers, label="references")
    return found, report


def write_kind(kind: str, writer: TsvWriter, seeded: set[str], citations: dict) -> dict:
    counts = {"rows": 0, "unseeded": 0, "no_pubmed_id": 0, "no_citation": 0}
    for values in copy_rows(interactions_sql(kind)):
        rec = sql_row_to_record(kind, values)
        gene_id = queried_gene(rec["subjectId"], rec["objectId"], seeded)
        if gene_id is None:
            counts["unseeded"] += 1
            continue
        citation = citations.get(rec["referenceId"], "")
        counts["no_pubmed_id"] += not rec["pubmedId"]
        counts["no_citation"] += not citation
        writer.write_row(**tsv_row(rec, gene_id, citation))
        counts["rows"] += 1
    log.info("%s: %d rows, %d interactions without a seeded participant, "
             "%d rows without pubmedId, %d without shortCitation", kind, counts["rows"],
             counts["unseeded"], counts["no_pubmed_id"], counts["no_citation"])
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance interactions from the curation DB into TSVs")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--mods", help="Comma-separated MODs to seed from (default: all 9, env ALLIANCE_MODS)")
    parser.add_argument("--limit", type=int, help="Seed at most N genes (for smoke tests)")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSVs (default: ../data)",
    )
    parser.add_argument("--no-reference-lookup", action="store_true",
                        help="Skip /reference/{curie}; leave shortCitation empty")
    add_parallel_args(parser)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    seeded = set(resolve_gene_ids(args.ids, args.mods, args.limit))
    out_dir = Path(args.out_dir)
    stats = FetchStats()

    citations: dict[str, str] = {}
    report = None
    if not args.no_reference_lookup:
        curies = [r[0] for r in copy_rows(references_sql())]
        citations, report = lookup_citations(curies, args.workers, stats)

    with TsvWriter(
        out_dir / "molecular-interactions.tsv",
        MOLECULAR_COLUMNS,
        release=release,
        source="curation DB genemolecularinteraction",
    ) as mol_writer, TsvWriter(
        out_dir / "genetic-interactions.tsv",
        GENETIC_COLUMNS,
        release=release,
        source="curation DB genegeneticinteraction",
    ) as gen_writer:
        mol = write_kind("molecular", mol_writer, seeded, citations)
        gen = write_kind("genetic", gen_writer, seeded, citations)
        if report is not None:
            finish_run(report, stats, out_dir, args.max_failure_rate,
                       rows_written=mol_writer.rows_written + gen_writer.rows_written)
        elif mol["rows"] + gen["rows"] == 0:
            log.error("interactions: no rows written")
            sys.exit(1)

    log.info("interactions: molecular rows=%d  genetic rows=%d", mol["rows"], gen["rows"])


if __name__ == "__main__":
    main()
