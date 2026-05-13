#!/usr/bin/env python3
"""
Fetch molecular variants from the Alliance FMS VARIANT-ALLELE-JSON dumps.

The bulk file ships both allele-only rows (Category="allele") and rows
that carry a populated VariantId (Category="allele with N known variants").
Only the latter yield Variant items. AlleleId is preserved as the linking
field so the integration engine merges into existing Allele rows.

FMS subType is `NCBITaxon{N}` — single integer with no colon.

Usage:
    python3 fetch_variants.py
    python3 fetch_variants.py --taxons 7955,10090
    python3 fetch_variants.py --limit 1000

Output: data/variants.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    get_current_release,
    iter_fms_datatype,
    join_pipe,
    normalize_taxon,
)

log = logging.getLogger("alliance.fetch.variants")


# Subset of taxons that VARIANT-ALLELE-JSON shards. FB uses NCBITaxon7227.
DEFAULT_TAXONS = [
    "NCBITaxon7227",   # D. melanogaster
    "NCBITaxon7955",   # D. rerio
    "NCBITaxon6239",   # C. elegans
    "NCBITaxon10090",  # M. musculus
    "NCBITaxon10116",  # R. norvegicus
    "NCBITaxon559292", # S. cerevisiae
    "NCBITaxon9606",   # H. sapiens
    "NCBITaxon8355",   # X. laevis
    "NCBITaxon8364",   # X. tropicalis
]


COLUMNS = [
    "variantId",
    "variantSymbol",
    "variantType",       # SO term name (e.g. "insertion")
    "variantTypeId",     # SO curie (e.g. "SO:0000667")
    "hgvsName",
    "chromosome",
    "startPosition",
    "endPosition",
    "refSequence",
    "varSequence",
    "consequenceNames",  # pipe-joined SO consequence names
    "alleleId",
    "geneId",            # pipe-joined affected gene ids
    "taxon",
    "dataProvider",
]


def _row(record: dict) -> dict | None:
    variant_id = record.get("VariantId")
    if not variant_id:
        return None
    return {
        "variantId": variant_id,
        "variantSymbol": record.get("VariantSymbol") or "",
        "variantType": record.get("VariantsTypeName") or "",
        "variantTypeId": record.get("VariantsTypeId") or "",
        "hgvsName": record.get("VariantsHgvsNames") or "",
        "chromosome": record.get("Chromosome") or "",
        "startPosition": str(record.get("StartPosition") or ""),
        "endPosition": str(record.get("EndPosition") or ""),
        "refSequence": record.get("SequenceOfReference") or "",
        "varSequence": record.get("SequenceOfVariant") or "",
        "consequenceNames": join_pipe(record.get("MostSevereConsequenceName") or []),
        "alleleId": record.get("AlleleId") or "",
        "geneId": join_pipe(record.get("VariantAffectedGeneId") or []),
        "taxon": normalize_taxon(record.get("Taxon") or ""),
        "dataProvider": (record.get("VariantCrossReferences") or [""])[0].split(":")[0]
                       if record.get("VariantCrossReferences") else "",
    }


def _iter_taxon(taxon_sub: str) -> Iterator[dict]:
    for _sub, rec in iter_fms_datatype("VARIANT-ALLELE-JSON", sub_filter=taxon_sub):
        row = _row(rec)
        if row is not None:
            yield row


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance variants (FMS bulk)")
    parser.add_argument(
        "--taxons",
        default=",".join(DEFAULT_TAXONS),
        help="Comma-separated FMS subType taxon strings (e.g. NCBITaxon10090)",
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

    taxons = [t.strip() for t in args.taxons.split(",") if t.strip()]
    release = get_current_release()
    out_dir = Path(args.out_dir)
    seen: set = set()
    total = 0

    with TsvWriter(
        out_dir / "variants.tsv",
        COLUMNS,
        release=release,
        source=f"FMS VARIANT-ALLELE-JSON ({','.join(taxons)})",
    ) as writer:
        for taxon in taxons:
            log.info("Variants FMS: %s", taxon)
            for row in _iter_taxon(taxon):
                key = row["variantId"]
                if key in seen:
                    continue
                seen.add(key)
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    log.info("Hit limit; stopping")
                    log.info("Done. rows=%d", total)
                    return

    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
