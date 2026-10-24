#!/usr/bin/env python3
"""
Build variants.tsv (see fetch_variants.py) from a restored Alliance curation
database.

One row per distinct variant HGVS name linked (live association) to an allele
in scope, the first allele winning, MODs in the fetcher's taxon order. Values
use the VARIANT-ALLELE formats (extract_alleles_common.load_variants):
variantSymbol is chromosome plus the HGVS change, sequences carry the padded
base, consequenceNames has one entry per gene-level predicted consequence
(its VEP terms comma-joined, most severe first), geneId the genes of those
consequences, and dataProvider the prefix of the variant's first cross
reference.

Scope: alleles of the release's FMS VARIANT-ALLELE files that the curation DB
still serves; --all-alleles takes every summary allele.

Usage:
    python3 extract_variants.py
    python3 extract_variants.py --mods ZFIN,MGI
    python3 extract_variants.py --limit 1000

Output: data/variants.tsv
"""

from __future__ import annotations

import argparse
import collections
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import MOD_TAXON, TsvWriter, configure_logging, get_current_release, normalize_taxon, selected_mods  # noqa: E402
from curation_db import copy_rows  # noqa: E402
from extract_alleles_common import (  # noqa: E402
    VARIANT_ALLELE_MODS,
    allele_variant_links,
    load_variants,
    mod_for,
    seed_allele_ids,
    summary_alleles,
)
from fetch_variants import COLUMNS, DEFAULT_TAXONS  # noqa: E402

log = logging.getLogger("alliance.extract.variants")

# Fetcher MOD order (its taxon list order decides which allele keeps a shared variant).
_TAXON_ORDER = [t.replace("NCBITaxon", "NCBITaxon:") for t in DEFAULT_TAXONS]
MOD_ORDER = sorted(MOD_TAXON, key=lambda m: _TAXON_ORDER.index(MOD_TAXON[m]))


def variant_row(v: dict, allele_id: str, taxon: str) -> dict:
    return {
        "variantId": v["hgvs"],
        "variantSymbol": v["symbol"],
        "variantType": v["type_name"],
        "variantTypeId": v["type_id"],
        "hgvsName": v["hgvs"],
        "chromosome": v["chromosome"],
        "startPosition": v["start"],
        "endPosition": v["end"],
        "refSequence": v["ref"],
        "varSequence": v["var"],
        "consequenceNames": "|".join(v["consequences"]),
        "alleleId": allele_id,
        "geneId": "|".join(v["gene_ids"]),
        "taxon": normalize_taxon(taxon),
        "dataProvider": v["xrefs"][0].split(":")[0] if v["xrefs"] else "",
    }


def iter_rows(mods, by_mod, links, variants):
    """(MOD, row) per distinct variant HGVS name, first allele first."""
    seen: set[str] = set()
    for mod in [m for m in MOD_ORDER if m in mods]:
        for db_id, pid, taxon in by_mod.get(mod, []):
            for variant_id in links.get(db_id, []):
                v = variants.get(variant_id)
                if v is None or not v["hgvs"] or v["hgvs"] in seen:
                    continue
                seen.add(v["hgvs"])
                yield mod, variant_row(v, pid, taxon)


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance variants from the curation DB")
    parser.add_argument("--mods", help="Comma-separated MODs (default: all, env ALLIANCE_MODS)")
    parser.add_argument("--all-alleles", action="store_true",
                        help="Variants of every summary allele, not only the FMS VARIANT-ALLELE seed")
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
    variants = load_variants(copy_rows)
    links = allele_variant_links(copy_rows)

    by_mod: dict[str, list[tuple[str, str, str]]] = collections.defaultdict(list)
    for db_id, pid, provider, taxon, _symbol, _fmt, _xref in summary_alleles(copy_rows, mods):
        if db_id in links and (wanted is None or pid in wanted):
            by_mod[mod_for(provider, taxon)].append((db_id, pid, taxon))

    per_mod: collections.Counter = collections.Counter()
    with TsvWriter(
        Path(args.out_dir) / "variants.tsv",
        COLUMNS,
        release=release,
        source=f"curation DB allele-variant associations ({','.join(mods)})",
    ) as writer:
        for mod, row in iter_rows(mods, by_mod, links, variants):
            writer.write_row(**row)
            per_mod[mod] += 1
            if args.limit and writer.rows_written >= args.limit:
                break
    log.info("variants rows per MOD: %s", dict(sorted(per_mod.items())))


if __name__ == "__main__":
    main()
