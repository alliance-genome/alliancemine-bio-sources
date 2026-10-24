#!/usr/bin/env python3
"""
Build data/expression-alliance.tsv from a restored Alliance curation
database, in the layout of the FMS EXPRESSION-ALLIANCE COMBINED file that
fetch_expression.py downloads.

AllianceExpressionConverter skips the first 15 lines of the file, so the
output keeps the FMS block of 14 comment lines plus the column header.

One row per distinct (gene, location, stage, assay, anatomical site) over
the non-obsolete, non-internal gene expression annotations. SourceURL holds
the comma-joined annotation (ZFIN figure) and expression experiment (MGI
assay, WB expression pattern) cross references, each filled into its
resourcedescriptorpage urlTemplate; a stage of "N/A" is
written blank, as FMS does. Reference holds the
comma-joined publications, a PMID where the reference has one, else the MOD
publication id.

Publications resolve through reference_crossreference.

Usage:
    python3 extract_expression.py
    python3 extract_expression.py --mods ZFIN --limit 1000
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    MOD_TAXON,
    configure_logging,
    get_current_release,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402

log = logging.getLogger("alliance.extract.expression")

DEFAULT_MODS = ["SGD", "MGI", "RGD", "ZFIN", "FB", "WB", "XBXL", "XBXT"]

COLUMNS = [
    "Species", "SpeciesID", "GeneID", "GeneSymbol", "Location", "StageTerm",
    "AssayID", "AssayTermName",
    "CellularComponentID", "CellularComponentTerm",
    "CellularComponentQualifierIDs", "CellularComponentQualifierTermNames",
    "SubStructureID", "SubStructureName",
    "SubStructureQualifierIDs", "SubStructureQualifierTermNames",
    "AnatomyTermID", "AnatomyTermName",
    "AnatomyTermQualifierIDs", "AnatomyTermQualifierTermNames",
    "SourceURL", "Source", "Reference",
]

# Species names as the FMS file prints them (the curation DB names the
# yeast taxon "Saccharomyces cerevisiae S288C").
SPECIES_NAME = {
    "SGD": "Saccharomyces cerevisiae",
    "MGI": "Mus musculus",
    "RGD": "Rattus norvegicus",
    "ZFIN": "Danio rerio",
    "FB": "Drosophila melanogaster",
    "WB": "Caenorhabditis elegans",
    "XBXL": "Xenopus laevis",
    "XBXT": "Xenopus tropicalis",
}

# Publication id prefix each data provider uses for papers without a PMID.
PUB_PREFIX = {
    "SGD": "SGD:", "MGI": "MGI:", "RGD": "RGD:", "ZFIN": "ZFIN:", "FB": "FB:",
    "WB": "WB:", "XBXL": "Xenbase:", "XBXT": "Xenbase:", "HUMAN": "PMID:",
}
PUB_PREFIX_PATTERN = "^(PMID|SGD|MGI|RGD|ZFIN|FB|WB|Xenbase):"


def pick_reference(curies: Iterable[str], mod: str) -> str:
    """The publication id of one reference: its PMID, else the providing
    MOD's own id, else any other MOD publication id."""
    curies = sorted(c for c in curies if c)
    for c in curies:
        if c.startswith("PMID:"):
            return c
    own = PUB_PREFIX.get(mod, "")
    for c in curies:
        if own and c.startswith(own):
            return c
    return curies[0] if curies else ""


def reference_curies(reference_ids_sql: str) -> dict[str, list[str]]:
    """PMID and MOD publication curies per reference id, for the reference
    ids that reference_ids_sql selects."""
    sql = f"""
        SELECT rx.reference_id, c.referencedcurie
        FROM reference_crossreference rx
        JOIN crossreference c ON c.id = rx.crossreferences_id
        WHERE rx.reference_id IN ({reference_ids_sql})
          AND c.referencedcurie ~ '{PUB_PREFIX_PATTERN}'
          AND NOT coalesce(c.obsolete, false)
    """
    out: dict[str, list[str]] = {}
    for ref_id, curie in copy_rows(sql):
        out.setdefault(ref_id, []).append(curie)
    return out


def taxon_mods(mods: list[str]) -> dict[str, str]:
    return {MOD_TAXON[m]: m for m in mods}


def _qualifiers(table: str, column: str) -> str:
    return f"""
        SELECT q.anatomicalsite_id AS sid,
               string_agg(t.curie, '|' ORDER BY t.curie) AS ids,
               string_agg(t.name, '|' ORDER BY t.curie) AS names
        FROM {table} q JOIN ontologyterm t ON t.id = q.{column}
        GROUP BY q.anatomicalsite_id"""


def annotation_filter_sql(taxa: list[str]) -> str:
    """Live annotations whose subject is a live gene of one of the taxa."""
    taxon_list = ", ".join(f"'{t}'" for t in taxa)
    return f"""
        SELECT a.*, g.primaryexternalid AS gene_id, tx.curie AS taxon
        FROM geneexpressionannotation a
        JOIN biologicalentity g ON g.id = a.expressionannotationsubject_id
        JOIN ontologyterm tx ON tx.id = g.taxon_id
        WHERE NOT coalesce(a.obsolete, false) AND NOT coalesce(a.internal, false)
          AND NOT coalesce(g.obsolete, false) AND NOT coalesce(g.internal, false)
          AND tx.curie IN ({taxon_list})"""


def expression_sql(taxa: list[str], limit: int | None) -> str:
    limit_sql = f"LIMIT {int(limit)}" if limit else ""
    return f"""
        WITH ann AS ({annotation_filter_sql(taxa)}),
        xref AS (
            SELECT ax.geneexpressionannotation_id AS aid, ax.crossreferences_id AS xid
            FROM geneexpressionannotation_crossreference ax
            UNION
            SELECT a.id, ex.crossreferences_id
            FROM ann a
            JOIN geneexpressionexperiment_crossreference ex
              ON ex.geneexpressionexperiment_id = a.expressionexperiment_id
        ),
        url AS (
            SELECT x.aid,
                   replace(p.urltemplate, '[%s]', regexp_replace(c.referencedcurie, '^[^:]*:', '')) AS url
            FROM xref x
            JOIN crossreference c ON c.id = x.xid
            JOIN resourcedescriptorpage p ON p.id = c.resourcedescriptorpage_id
            WHERE p.urltemplate IS NOT NULL
        ),
        symbol AS (
            SELECT singlegene_id AS gid, min(displaytext) AS symbol
            FROM slotannotation
            WHERE slotannotationtype = 'GeneSymbolSlotAnnotation' AND NOT coalesce(obsolete, false)
            GROUP BY singlegene_id
        ),
        ccq AS ({_qualifiers('anatomicalsite_cellularcomponentqualifiers', 'cellularcomponentqualifiers_id')}),
        subq AS ({_qualifiers('anatomicalsite_anatomicalsubstructurequalifiers', 'anatomicalsubstructurequalifiers_id')}),
        anq AS ({_qualifiers('anatomicalsite_anatomicalstructurequalifiers', 'anatomicalstructurequalifiers_id')})
        SELECT a.taxon, a.gene_id, coalesce(sy.symbol, ''),
               coalesce(a.whereexpressedstatement, ''),
               coalesce(nullif(a.whenexpressedstagename, 'N/A'), ''),
               coalesce(mmo.curie, ''), coalesce(mmo.name, ''),
               coalesce(cc.curie, ''), coalesce(cc.name, ''), coalesce(ccq.ids, ''), coalesce(ccq.names, ''),
               coalesce(sub.curie, ''), coalesce(sub.name, ''), coalesce(subq.ids, ''), coalesce(subq.names, ''),
               coalesce(an.curie, ''), coalesce(an.name, ''), coalesce(anq.ids, ''), coalesce(anq.names, ''),
               string_agg(DISTINCT u.url, '|'),
               string_agg(DISTINCT a.evidenceitem_id::text, '|')
        FROM ann a
        LEFT JOIN symbol sy ON sy.gid = a.expressionannotationsubject_id
        LEFT JOIN ontologyterm mmo ON mmo.id = a.expressionassayused_id
        LEFT JOIN expressionpattern ep ON ep.id = a.expressionpattern_id
        LEFT JOIN anatomicalsite s ON s.id = ep.whereexpressed_id
        LEFT JOIN ontologyterm cc ON cc.id = s.cellularcomponentterm_id
        LEFT JOIN ontologyterm sub ON sub.id = s.anatomicalsubstructure_id
        LEFT JOIN ontologyterm an ON an.id = s.anatomicalstructure_id
        LEFT JOIN ccq ON ccq.sid = s.id
        LEFT JOIN subq ON subq.sid = s.id
        LEFT JOIN anq ON anq.sid = s.id
        LEFT JOIN url u ON u.aid = a.id
        GROUP BY 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
        ORDER BY 1, 2
        {limit_sql}"""


def build_row(rec: list[str], mod: str, refs: dict[str, list[str]]) -> list[str]:
    """One TSV row from an expression_sql result row."""
    taxon, gene_id, symbol, *middle, urls, ref_ids = rec
    pubs: list[str] = []
    for ref_id in sorted(r for r in ref_ids.split("|") if r):
        pub = pick_reference(refs.get(ref_id, []), mod)
        if pub and pub not in pubs:
            pubs.append(pub)
    return [SPECIES_NAME.get(mod, ""), taxon, gene_id, symbol, *middle,
            ",".join(sorted(u for u in urls.split("|") if u)), mod, ",".join(sorted(pubs))]


def header_lines(mods: list[str], release: str) -> list[str]:
    """The 14 comment lines that precede the column header in the FMS file."""
    rule = "#" * 74
    return [
        rule,
        "#",
        "# Data type: Expression",
        "# Data format: tsv",
        "# README: ",
        "# Source: Alliance curation database (geneexpressionannotation)",
        "# Source URL: http://alliancegenome.org/downloads",
        "# Help Desk: help@alliancegenome.org",
        "# Taxon IDs: " + ", ".join(MOD_TAXON[m] for m in mods),
        "# Species: " + ", ".join(SPECIES_NAME[m] for m in mods),
        f"# Alliance Database Version: {release}",
        "# Date file generated (UTC): " + time.strftime("%Y-%m-%d %H:%M", time.gmtime()),
        "#",
        rule,
    ]


def clean(value: str) -> str:
    return value.replace("\t", " ").replace("\n", " ").replace("\r", " ")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance expression TSV from the curation DB")
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

    mods = [m for m in selected_mods(DEFAULT_MODS, args.mods) if m in SPECIES_NAME]
    if not mods:
        log.error("No expression MODs selected (expression covers %s)", ",".join(DEFAULT_MODS))
        sys.exit(1)
    by_taxon = taxon_mods(mods)
    release = get_current_release()
    taxa = list(by_taxon)

    started = time.monotonic()
    refs = reference_curies(
        f"SELECT evidenceitem_id FROM ({annotation_filter_sql(taxa)}) a")
    log.info("Publication ids for %d references (%.0fs)", len(refs), time.monotonic() - started)

    out_path = Path(args.out_dir) / "expression-alliance.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    rows = 0
    try:
        with tmp.open("w", encoding="utf-8") as out:
            for line in header_lines(mods, release):
                out.write(line + "\n")
            out.write("\t".join(COLUMNS) + "\n")
            for rec in copy_rows(expression_sql(taxa, args.limit)):
                row = build_row(rec, by_taxon[rec[0]], refs)
                out.write("\t".join(clean(v) for v in row) + "\n")
                rows += 1
        os.replace(tmp, out_path)
    finally:
        tmp.unlink(missing_ok=True)
    log.info("Wrote %s (%d data rows, %.0fs)", out_path, rows, time.monotonic() - started)


if __name__ == "__main__":
    main()
