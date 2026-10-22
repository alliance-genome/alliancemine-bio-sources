#!/usr/bin/env python3
"""
Build transgenic-alleles.tsv (see fetch_transgenic_alleles.py) from a restored
Alliance curation database.

A transgenic allele is a live allele (not obsolete, not internal) with a live
allele-construct association, as the fetcher's FMS path takes ALLELE records
with a "contains" construct relation. Columns:
    geneId                  first live is_allele_of gene
    geneSymbol              empty, as in the FMS path
    geneTaxon               allele taxon curie
    transgenicAlleleId      allele primary external id
    transgenicAlleleSymbol  symbol formatText (the FMS symbolText form)
    constructIds            construct primary external ids (else MOD internal ids),
                            pipe-joined by association id
    dataProvider            fetcher MOD name

Usage:
    python3 extract_transgenic_alleles.py
    python3 extract_transgenic_alleles.py --mods MGI,ZFIN
    python3 extract_transgenic_alleles.py --ids MGI:88276
    python3 extract_transgenic_alleles.py --limit 100

Output: data/transgenic-alleles.tsv
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
from extract_alleles_common import first_genes, mod_for, summary_alleles  # noqa: E402
from fetch_transgenic_alleles import COLUMNS, DEFAULT_MODS  # noqa: E402

log = logging.getLogger("alliance.extract.transgenic_alleles")


# AlleleDAO construct ids: live allele-construct associations, construct
# primary external id else MOD internal id.
CONSTRUCTS_SQL = """
    SELECT aca.alleleassociationsubject_id, coalesce(r.primaryexternalid, r.modinternalid)
    FROM alleleconstructassociation aca
    JOIN reagent r ON r.id = aca.alleleconstructassociationobject_id
    WHERE aca.obsolete = false AND aca.internal = false
      AND (r.primaryexternalid IS NOT NULL OR r.modinternalid IS NOT NULL)
    ORDER BY aca.alleleassociationsubject_id, aca.id
"""


def transgenic_row(pid: str, mod: str, taxon: str, symbol: str, construct_ids: list[str], gene_id: str) -> dict:
    return {
        "geneId": gene_id,
        "geneSymbol": "",
        "geneTaxon": normalize_taxon(taxon),
        "transgenicAlleleId": pid,
        "transgenicAlleleSymbol": symbol,
        "constructIds": "|".join(dict.fromkeys(construct_ids)),
        "dataProvider": mod,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance transgenic alleles from the curation DB")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs; keep alleles of these genes")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--limit", type=int, help="Emit at most N alleles")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    release = get_current_release()
    mods = selected_mods(DEFAULT_MODS, args.mods)
    genes = {g.strip() for g in args.ids.split(",") if g.strip()} if args.ids else None

    constructs: dict[str, list[str]] = collections.defaultdict(list)
    for allele_id, construct_id in copy_rows(CONSTRUCTS_SQL):
        constructs[allele_id].append(construct_id)
    first_gene = first_genes(copy_rows)
    per_mod: collections.Counter = collections.Counter()
    seen: set[str] = set()
    with TsvWriter(
        Path(args.out_dir) / "transgenic-alleles.tsv",
        COLUMNS,
        release=release,
        source=f"curation DB allele-construct associations ({','.join(mods)})",
    ) as writer:
        for db_id, pid, provider, taxon, _symbol, symbol, _xref in summary_alleles(copy_rows, mods):
            mod = mod_for(provider, taxon)
            construct_ids = constructs.get(db_id)
            gene_id = first_gene.get(db_id, ("", ""))[0]
            if (not construct_ids or pid in seen or mod not in mods
                    or (genes is not None and gene_id not in genes)):
                continue
            seen.add(pid)
            writer.write_row(**transgenic_row(pid, mod, taxon, symbol, construct_ids, gene_id))
            per_mod[mod] += 1
            if args.limit and writer.rows_written >= args.limit:
                break
    log.info("transgenic_alleles rows per MOD: %s", dict(sorted(per_mod.items())))


if __name__ == "__main__":
    main()
