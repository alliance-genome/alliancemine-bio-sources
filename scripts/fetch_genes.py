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
    python3 fetch_genes.py
    python3 fetch_genes.py --ids SGD:S000004103,SGD:S000002429
    python3 fetch_genes.py --limit 100   # smoke test first 100 yeast genes

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
    configure_logging,
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    http_get_json,
    open_cache,
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


def _extract_crossrefs(gene: dict) -> list[str]:
    """Pull compact 'DB:ID' tokens out of crossReferenceMap values.

    The API surface is nested: crossReferenceMap: {slot: {name: 'DB:ID', ...}}.
    We preserve the 'name' field of each entry that looks like 'TYPE:VALUE'."""
    cr_map = gene.get("crossReferenceMap") or {}
    tokens: list[str] = []
    for slot, entry in cr_map.items():
        if not isinstance(entry, dict):
            continue
        name = entry.get("name") or ""
        if ":" in name and not name.startswith("["):
            tokens.append(name)
    return tokens


def _row_for_gene(gene: dict) -> dict:
    primary_id = gene.get("id", "")
    secondary_ids = gene.get("secondaryIds") or []
    secondary = secondary_ids[0] if secondary_ids else ""

    synonyms_list = gene.get("synonyms") or []
    crossref_list = _extract_crossrefs(gene)

    loc = (gene.get("genomeLocations") or [{}])[0]
    so_name = get_nested(gene, "soTerm", "name", default="")

    species = gene.get("species") or {}
    species_taxon = species.get("taxonId", "")

    return {
        "primaryIdentifier": primary_id,
        "secondaryIdentifier": secondary,
        "synonyms": _list_literal(synonyms_list),
        "crossrefs": _list_literal(crossref_list),
        "name": gene.get("name", ""),
        "symbol": gene.get("symbol", ""),
        "modDescription": gene.get("geneSynopsis", ""),
        "automatedDescription": gene.get("automatedGeneSynopsis", ""),
        "species": species_taxon,
        "chromosome": loc.get("chromosome", ""),
        "start": loc.get("start", "") or "",
        "end": loc.get("end", "") or "",
        "strand": loc.get("strand", ""),
        "featureType": so_name,
        "dateProduced": gene.get("dateProduced", ""),
        "dataProvider": gene.get("dataProvider", ""),
        "modCrossRefCompleteUrl": gene.get("modCrossRefCompleteUrl", ""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance gene metadata into TSV")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV (default: ../data)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    if args.ids:
        gene_ids = [g.strip() for g in args.ids.split(",") if g.strip()]
    else:
        gene_ids = enumerate_yeast_genes()
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
    with open_cache("genes") as cache, \
         TsvWriter(
             out_dir / "alliance-genes.tsv",
             COLUMNS,
             release=release,
             source="/gene/{id}",
         ) as writer:
        for gid in iterator:
            try:
                gene = http_get_json(f"/gene/{gid}", cache=cache, stats=stats)
            except Exception as e:
                log.warning("Gene fetch failed for %s: %s", gid, e)
                continue
            writer.write_row(**_row_for_gene(gene))

    log.info("Done. API stats: %s", stats)


if __name__ == "__main__":
    main()
