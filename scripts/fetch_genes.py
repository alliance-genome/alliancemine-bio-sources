#!/usr/bin/env python3
"""
Fetch per-gene metadata from the Alliance API, emit an enriched TSV.

The TSV layout keeps the 14-column schema that AllianceGenesConverter already
understands (so an in-place swap of data sources is possible) and appends three
new trailing columns for attributes that FMS never surfaced:

    14  dateProduced            ISO 8601 UTC timestamp from the API
    15  dataProvider            MOD abbreviation, e.g. "SGD"
    16  modCrossRefCompleteUrl  direct link to the MOD gene page

Usage:
    python3 fetch_genes.py                       # every gene of the 9 MODs
    python3 fetch_genes.py --mods SGD,MGI --workers 6
    python3 fetch_genes.py --ids SGD:S000004103,HGNC:1100
    python3 fetch_genes.py --limit 100           # smoke test

Output: data/alliance-genes.tsv
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
    http_get_json,
    resolve_gene_ids,
    run_parallel,
)

log = logging.getLogger("alliance.fetch.genes")


# Column order must match AllianceGenesConverter COL_* constants so the existing
# Java code doesn't have to move. Trailing three columns are the new additions.
COLUMNS = [
    "primaryIdentifier",        # 0
    "secondaryIdentifier",      # 1
    "synonyms",                 # 2  (Python-list format: "[a, b, c]" or "[]")
    "crossrefs",                # 3  (same format)
    "name",                     # 4
    "symbol",                   # 5
    "modDescription",           # 6
    "automatedDescription",     # 7
    "species",                  # 8  "NCBITaxon:559292"
    "chromosome",               # 9
    "start",                    # 10
    "end",                      # 11
    "strand",                   # 12
    "featureType",              # 13
    "dateProduced",             # 14  (new)
    "dataProvider",             # 15  (new)
    "modCrossRefCompleteUrl",   # 16  (new)
]


def _list_literal(values) -> str:
    """Format a Python list like "[a, b, c]" that the existing
    AllianceGenesConverter.getSynonyms() / getCrossReference() parsers expect."""
    if not values:
        return "[]"
    return "[" + ", ".join(str(v) for v in values if v) + "]"


def _display(node: dict | None) -> str:
    return get_nested(node, "displayText", default="") or ""


def _note(gene: dict, note_type: str) -> str:
    for note in gene.get("relatedNotes") or []:
        if get_nested(note, "noteType", "name") == note_type:
            return note.get("freeText", "") or ""
    return ""


def _mod_url(xref: dict | None) -> str:
    """Fill the resource descriptor's [%s] with the local part of the curie."""
    curie = get_nested(xref, "referencedCurie", default="") or ""
    template = get_nested(xref, "resourceDescriptorPage", "urlTemplate", default="") or ""
    if not curie or "[%s]" not in template:
        return ""
    return template.replace("[%s]", curie.split(":", 1)[-1])


def _row_for_gene(payload: dict) -> dict:
    """Map the /gene/{id} payload ({"gene": {...curation Gene...}}) onto the
    converter's columns. Missing location fields are the literal "null" the
    converter tests for, as in bgi_to_genes_tsv.py."""
    gene = payload.get("gene") or {}
    secondary_ids = [s.get("secondaryId", "") for s in gene.get("geneSecondaryIds") or []]
    # Systematic name (SGD ORF, WormBase sequence name) and the secondary IDs after the first
    # have no Gene attribute of their own; they become synonyms so lookups find the gene.
    synonyms = [_display(s) for s in gene.get("geneSynonyms") or []]
    synonyms += [_display(gene.get("geneSystematicName"))] + secondary_ids[1:]
    synonyms = list(dict.fromkeys(v for v in synonyms if v))
    crossrefs: list[str] = []
    for x in gene.get("crossReferences") or []:
        curie = x.get("referencedCurie", "")
        if curie and curie not in crossrefs:
            crossrefs.append(curie)
    locs = gene.get("geneGenomicLocationAssociations") or []
    loc = locs[0] if locs else {}

    def pos(value) -> str:
        return "null" if value in (None, "") else str(value)

    return {
        "primaryIdentifier": gene.get("primaryExternalId", ""),
        "secondaryIdentifier": secondary_ids[0] if secondary_ids else "",
        "synonyms": _list_literal(synonyms),
        "crossrefs": _list_literal(crossrefs),
        "name": _display(gene.get("geneFullName")),
        "symbol": _display(gene.get("geneSymbol")),
        "modDescription": _note(gene, "MOD_provided_gene_description"),
        "automatedDescription": _note(gene, "automated_gene_description"),
        "species": get_nested(gene, "taxon", "curie", default=""),
        "chromosome": pos(get_nested(loc, "geneGenomicLocationAssociationObject", "name")),
        "start": pos(loc.get("start")),
        "end": pos(loc.get("end")),
        "strand": pos(loc.get("strand")),
        "featureType": get_nested(gene, "geneType", "name", default=""),
        "dateProduced": gene.get("dateUpdated", "") or gene.get("dateCreated", "") or "",
        "dataProvider": get_nested(gene, "dataProvider", "abbreviation", default=""),
        "modCrossRefCompleteUrl": _mod_url(gene.get("dataProviderCrossReference")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance gene metadata into TSV")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--mods", help="Comma-separated MODs to seed from (default: all 9, env ALLIANCE_MODS)")
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
    out_dir = Path(args.out_dir)
    stats = FetchStats()

    def work(gid, cache, stats):
        row = _row_for_gene(http_get_json(f"/gene/{gid}", cache=cache, stats=stats))
        if not row["primaryIdentifier"]:
            raise ValueError("response has no gene.primaryExternalId")
        return row

    with TsvWriter(
        out_dir / "alliance-genes.tsv",
        COLUMNS,
        release=release,
        source="/gene/{id}",
    ) as writer:
        report = run_parallel(
            gene_ids, work, lambda _gid, row: writer.write_row(**row),
            cache_name="genes", stats=stats, workers=args.workers, label="genes",
        )
        finish_run(report, stats, out_dir, args.max_failure_rate,
                   rows_written=writer.rows_written)


if __name__ == "__main__":
    main()
