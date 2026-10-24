#!/usr/bin/env python3
"""
Build alliance-genes.tsv from a restored Alliance curation database dump.

Same file, columns and value formats as fetch_genes.py (see its COLUMNS), so
AllianceGenesConverter reads either. One set-based query for all MODs
replaces 394k /gene/{id} calls. Differences from the API path:

    automatedDescription  The dump holds no automated_gene_description notes
                          (the API merges them in), so the column comes from
                          the release's FMS GENE-DESCRIPTION-JSON "description".
    dateProduced          The API omits dateUpdated, so fetch_genes.py leaves
                          it empty; here it is dateUpdated (else dateCreated),
                          read as UTC, in the converter's
                          yyyy-MM-dd'T'HH:mm:ss.SSSXXX form.
    synonyms, crossrefs   Same members; order follows slot annotation and
                          crossreference ids, the API's order is unspecified.

The gene set is the fetcher's: the per-MOD gene lists (resolve_gene_ids).
Listed ids the database does not yield are fetched from the API with
fetch_genes.py's own mapping and appended, so no gene is lost.

Usage:
    python3 extract_genes.py                     # every gene of the 9 MODs
    python3 extract_genes.py --mods SGD,MGI
    python3 extract_genes.py --ids SGD:S000004103,HGNC:1100
    python3 extract_genes.py --limit 100         # smoke test

Environment: CURATION_PSQL or CURATION_DB_URL (see curation_db.py).
Output: data/alliance-genes.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    ALL_MODS,
    MOD_TAXON,
    FetchStats,
    TsvWriter,
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    http_get_json,
    iter_fms_datatype,
    resolve_gene_ids,
    run_parallel,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402
from fetch_genes import COLUMNS, _list_literal, _row_for_gene  # noqa: E402

log = logging.getLogger("alliance.extract.genes")

# string_agg separator for multi-valued columns; never occurs in curated text.
SEP = "\x1f"

# Query result columns, in SELECT order.
DB_FIELDS = [
    "primaryExternalId", "taxon", "dataProvider", "geneType", "dateUpdated", "dateCreated",
    "symbol", "fullName", "systematicName", "synonyms", "secondaryIds", "crossrefs",
    "modDescription", "chromosome", "start", "end", "strand", "modXrefCurie", "modXrefTemplate",
]

# Sequence and regulatory features the gene table also holds (1.3M rows in
# 9.2.0) that no MOD gene list carries. Dropping them early keeps the joins
# small; the gene-list filter in main() stays the authority either way.
# biological_region stays: 8 listed HGNC genes carry it.
FEATURE_TYPES = (
    "ARS", "CTCF_binding_site", "CpG_island", "E_box_motif", "LTR_retrotransposon",
    "TF_binding_site", "TSS_region", "base_call_error_correction", "binding_site",
    "centromere", "chromosome_part", "endogenous_retroviral_sequence",
    "enhancer", "gene_group", "histone_binding_site", "imprinting_control_region", "insulator",
    "intein_encoding_region", "interchromosomal_duplication", "intrachromosomal_duplication",
    "intronic_regulatory_region", "locus_control_region", "long_terminal_repeat",
    "mating_type_region", "matrix_attachment_site", "miRNA_target_site", "minisatellite",
    "open_chromatin_region", "origin_of_replication", "polyA_signal_sequence", "polyA_site",
    "possible_base_call_error", "promoter", "protein_binding_site", "recombination_enhancer",
    "regulatory_region", "response_element", "retrotransposon", "sequence_feature",
    "silencer", "silent_mating_type_cassette_array", "telomere", "three_prime_UTR",
    "trans_splice_acceptor_site", "transcription_start_cluster",
    "transcriptional_cis_regulatory_region",
)

# Live genes of the requested taxa; extract_gene_crossrefs.py reuses it.
GENE_BASE_CTE = """g AS (
  SELECT b.id, b.primaryexternalid, t.curie AS taxon, o.abbreviation AS provider,
         gt.name AS genetype, b.dateupdated, b.datecreated,
         b.dataprovidercrossreference_id, ge.gcrpcrossreference_id
  FROM gene ge
  JOIN biologicalentity b ON b.id = ge.id
  JOIN ontologyterm t ON t.id = b.taxon_id
  LEFT JOIN organization o ON o.id = b.dataprovider_id
  LEFT JOIN ontologyterm gt ON gt.id = ge.genetype_id
  WHERE t.curie IN ({taxa}) {id_filter}
    AND NOT coalesce(b.obsolete, false) AND NOT coalesce(b.internal, false)
    AND coalesce(gt.name, '') NOT IN ({feature_types})
)"""

GENE_SQL = """
WITH {gene_base},
slots AS (
  SELECT s.singlegene_id AS id,
    min(s.displaytext) FILTER (WHERE s.slotannotationtype = 'GeneSymbolSlotAnnotation') AS symbol,
    min(s.displaytext) FILTER (WHERE s.slotannotationtype = 'GeneFullNameSlotAnnotation') AS fullname,
    min(s.displaytext) FILTER (WHERE s.slotannotationtype = 'GeneSystematicNameSlotAnnotation') AS sysname,
    string_agg(s.displaytext, E'\\x1f' ORDER BY s.id)
      FILTER (WHERE s.slotannotationtype = 'GeneSynonymSlotAnnotation') AS synonyms,
    string_agg(s.secondaryid, E'\\x1f' ORDER BY s.id)
      FILTER (WHERE s.slotannotationtype = 'GeneSecondaryIdSlotAnnotation') AS secondary
  FROM slotannotation s JOIN g ON g.id = s.singlegene_id
  WHERE s.slotannotationtype IN ('GeneSymbolSlotAnnotation', 'GeneFullNameSlotAnnotation',
        'GeneSystematicNameSlotAnnotation', 'GeneSynonymSlotAnnotation',
        'GeneSecondaryIdSlotAnnotation')
    AND NOT coalesce(s.obsolete, false) AND NOT coalesce(s.internal, false)
  GROUP BY s.singlegene_id
),
xrefs AS (
  SELECT gc.genomicentity_id AS id,
         string_agg(c.referencedcurie, E'\\x1f' ORDER BY c.id) AS curies
  FROM genomicentity_crossreference gc
  JOIN g ON g.id = gc.genomicentity_id
  JOIN crossreference c ON c.id = gc.crossreferences_id
  WHERE NOT coalesce(c.obsolete, false) AND NOT coalesce(c.internal, false)
  GROUP BY gc.genomicentity_id
),
notes AS (
  SELECT bn.submittedobject_id AS id, min(n.freetext ORDER BY n.id) AS mod_description
  FROM biologicalentity_note bn
  JOIN g ON g.id = bn.submittedobject_id
  JOIN note n ON n.id = bn.relatednotes_id
  JOIN vocabularyterm vt ON vt.id = n.notetype_id
  WHERE vt.name = 'MOD_provided_gene_description' AND NOT coalesce(n.obsolete, false)
  GROUP BY bn.submittedobject_id
),
loc AS (
  SELECT DISTINCT ON (l.geneassociationsubject_id)
         l.geneassociationsubject_id AS id, ac.name AS chromosome, l.start, l."end", l.strand
  FROM genegenomiclocationassociation l
  JOIN g ON g.id = l.geneassociationsubject_id
  LEFT JOIN assemblycomponent ac ON ac.id = l.genegenomiclocationassociationobject_id
  WHERE NOT coalesce(l.obsolete, false) AND NOT coalesce(l.internal, false)
  ORDER BY l.geneassociationsubject_id, l.id
)
SELECT g.primaryexternalid, g.taxon, g.provider, g.genetype,
       to_char(g.dateupdated, 'YYYY-MM-DD"T"HH24:MI:SS.MS"Z"'),
       to_char(g.datecreated, 'YYYY-MM-DD"T"HH24:MI:SS.MS"Z"'),
       s.symbol, s.fullname, s.sysname, s.synonyms, s.secondary, x.curies,
       n.mod_description, loc.chromosome, loc.start, loc."end", loc.strand,
       dx.referencedcurie, rp.urltemplate
FROM g
LEFT JOIN slots s ON s.id = g.id
LEFT JOIN xrefs x ON x.id = g.id
LEFT JOIN notes n ON n.id = g.id
LEFT JOIN loc ON loc.id = g.id
LEFT JOIN crossreference dx ON dx.id = g.dataprovidercrossreference_id
LEFT JOIN resourcedescriptorpage rp ON rp.id = dx.resourcedescriptorpage_id
ORDER BY g.primaryexternalid
"""


def _sql_list(values: list[str]) -> str:
    return ", ".join("'" + v.replace("'", "''") + "'" for v in values)


def gene_base_cte(taxa: list[str], ids: list[str] | None = None) -> str:
    id_filter = f"AND b.primaryexternalid IN ({_sql_list(ids)})" if ids else ""
    return GENE_BASE_CTE.format(taxa=_sql_list(taxa), id_filter=id_filter,
                                 feature_types=_sql_list(list(FEATURE_TYPES)))


def gene_sql(taxa: list[str], ids: list[str] | None = None) -> str:
    return GENE_SQL.format(gene_base=gene_base_cte(taxa, ids))


def _split(value: str) -> list[str]:
    return [v for v in value.split(SEP) if v] if value else []


def _mod_url(curie: str, template: str) -> str:
    """Same rule as fetch_genes._mod_url: [%s] takes the curie's local part."""
    if not curie or "[%s]" not in template:
        return ""
    return template.replace("[%s]", curie.split(":", 1)[-1])


def row_from_db(rec: dict, auto_descriptions: dict[str, str]) -> dict:
    """Map one GENE_SQL result (keyed by DB_FIELDS) onto fetch_genes COLUMNS,
    following fetch_genes._row_for_gene rule for rule."""
    secondary_ids = _split(rec["secondaryIds"])
    synonyms = _split(rec["synonyms"]) + [rec["systematicName"]] + secondary_ids[1:]
    synonyms = list(dict.fromkeys(v for v in synonyms if v))
    crossrefs = list(dict.fromkeys(_split(rec["crossrefs"])))

    def pos(value: str) -> str:
        return value if value else "null"

    gene_id = rec["primaryExternalId"]
    return {
        "primaryIdentifier": gene_id,
        "secondaryIdentifier": secondary_ids[0] if secondary_ids else "",
        "synonyms": _list_literal(synonyms),
        "crossrefs": _list_literal(crossrefs),
        "name": rec["fullName"],
        "symbol": rec["symbol"],
        "modDescription": rec["modDescription"],
        "automatedDescription": auto_descriptions.get(gene_id, ""),
        "species": rec["taxon"],
        "chromosome": pos(rec["chromosome"]),
        "start": pos(rec["start"]),
        "end": pos(rec["end"]),
        "strand": pos(rec["strand"]),
        "featureType": rec["geneType"],
        "dateProduced": rec["dateUpdated"] or rec["dateCreated"],
        "dataProvider": rec["dataProvider"],
        "modCrossRefCompleteUrl": _mod_url(rec["modXrefCurie"], rec["modXrefTemplate"]),
    }


def iter_db_genes(taxa: list[str], ids: list[str] | None = None) -> Iterator[dict]:
    for values in copy_rows(gene_sql(taxa, ids)):
        yield dict(zip(DB_FIELDS, values))


def load_auto_descriptions(mods: list[str]) -> dict[str, str]:
    """gene id -> automated description from FMS GENE-DESCRIPTION-JSON."""
    descriptions: dict[str, str] = {}
    for mod in mods:
        for _sub, rec in iter_fms_datatype("GENE-DESCRIPTION-JSON", sub_filter=mod):
            gene_id, text = rec.get("gene_id"), rec.get("description")
            if gene_id and text:
                descriptions.setdefault(gene_id, text)
    log.info("Automated descriptions: %d genes from FMS GENE-DESCRIPTION-JSON", len(descriptions))
    return descriptions


def _mods_for_ids(ids: list[str]) -> list[str]:
    """MODs whose taxon the --ids may belong to (prefix alone cannot split XBXL/XBXT)."""
    prefixes = {i.split(":", 1)[0] for i in ids}
    by_prefix = {"SGD": ["SGD"], "MGI": ["MGI"], "RGD": ["RGD"], "ZFIN": ["ZFIN"], "FB": ["FB"],
                 "WB": ["WB"], "Xenbase": ["XBXL", "XBXT"], "HGNC": ["HUMAN"]}
    return [m for m in ALL_MODS if any(m in by_prefix.get(p, []) for p in prefixes)]


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance gene metadata from the curation DB into TSV")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--mods", help="Comma-separated MODs (default: all 9, env ALLIANCE_MODS)")
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV (default: ../data)",
    )
    add_parallel_args(parser)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit)
    wanted = set(gene_ids)
    mods = _mods_for_ids(gene_ids) if args.ids else selected_mods(None, args.mods)
    # Small id sets go into the SQL; full runs filter the per-taxon result in Python.
    sql_ids = gene_ids if len(gene_ids) <= 5000 else None
    auto = load_auto_descriptions(mods)
    out_dir = Path(args.out_dir)
    start = time.time()

    with TsvWriter(
        out_dir / "alliance-genes.tsv",
        COLUMNS,
        release=release,
        source="curation DB dump (gene tables) + FMS GENE-DESCRIPTION-JSON",
    ) as writer:
        found: set[str] = set()
        per_taxon: dict[str, int] = {}
        unlisted = 0
        for rec in iter_db_genes([MOD_TAXON[m] for m in mods], sql_ids):
            gene_id = rec["primaryExternalId"]
            if gene_id not in wanted:
                unlisted += 1
                continue
            if gene_id in found:
                continue
            found.add(gene_id)
            writer.write_row(**row_from_db(rec, auto))
            per_taxon[rec["taxon"]] = per_taxon.get(rec["taxon"], 0) + 1
        for mod in mods:
            log.info("%s: %d genes from the database", mod, per_taxon.get(MOD_TAXON[mod], 0))
        log.info("Database pass: %d genes in %.0fs; skipped %d live genes outside the gene lists",
                 len(found), time.time() - start, unlisted)

        if not found:
            log.error("the database yielded none of the %d listed genes", len(gene_ids))
            sys.exit(1)
        missing = [g for g in gene_ids if g not in found]
        log.info("API fallback: %d listed genes the database did not yield", len(missing))
        stats = FetchStats()

        def work(gid, cache, stats):
            row = _row_for_gene(http_get_json(f"/gene/{gid}", cache=cache, stats=stats))
            if not row["primaryIdentifier"]:
                raise ValueError("response has no gene.primaryExternalId")
            return row

        if missing:
            report = run_parallel(
                missing, work, lambda _gid, row: writer.write_row(**row),
                cache_name="genes", stats=stats, workers=args.workers, label="genes-api-fallback",
            )
            finish_run(report, stats, out_dir, args.max_failure_rate,
                       rows_written=writer.rows_written)
        log.info("Done. rows=%d (db=%d, api=%d) in %.0fs",
                 writer.rows_written, len(found), writer.rows_written - len(found),
                 time.time() - start)


if __name__ == "__main__":
    main()
