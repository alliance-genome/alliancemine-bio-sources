#!/usr/bin/env python3
"""
Fetch alliance allele records as a TSV the existing
AllianceAllelesConverter can consume.

The converter was originally built around the legacy FMS
`VARIANT-ALLELE_COMBINED.tsv` file (frozen at release 4.0.0). The
current populated data lives in the `VARIANT-ALLELE-JSON` per-taxon
shards; this fetcher reads those JSON shards and re-renders them into
the 27-column TSV the converter expects.

TSV column order (must match AllianceAllelesConverter.process):
    0  Taxon                              (e.g. "NCBITaxon:7227")
    1  SpeciesName                        (e.g. "Drosophila melanogaster")
    2  AlleleId
    3  AlleleSymbol
    4  AlleleSynonyms                     (pipe-joined)
    5  VariantId
    6  VariantSymbol
    7  VariantSynonyms                    (pipe-joined)
    8  VariantCrossReferences             (pipe-joined)
    9  AlleleAssociatedGeneId             (first)
   10  AlleleAssociatedGeneSymbol         (first)
   11  VariantAffectedGeneId              (first)
   12  VariantAffectedGeneSymbol          (first)
   13  Category
   14  VariantsTypeId
   15  VariantsTypeName
   16  VariantsHgvsNames
   17  Assembly
   18  Chromosome
   19  StartPosition
   20  EndPosition
   21  SequenceOfReference
   22  SequenceOfVariant
   23  MostSevereConsequenceName          (pipe-joined)
   24  VariantInformationReference        (pipe-joined, first column form)
   25  HasDiseaseAnnotations
   26  HasPhenotypeAnnotations

Usage:
    python3 fetch_alleles.py
    python3 fetch_alleles.py --taxa NCBITaxon10090 --limit 1000

Output: data/alliance-alleles.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    get_current_release,
    iter_fms_datatype,
    join_pipe,
    normalize_taxon,
)

log = logging.getLogger("alliance.fetch.alleles")


DEFAULT_TAXA = [
    "NCBITaxon7227",
    "NCBITaxon7955",
    "NCBITaxon6239",
    "NCBITaxon10090",
    "NCBITaxon10116",
    "NCBITaxon559292",
    "NCBITaxon9606",
    "NCBITaxon8355",
    "NCBITaxon8364",
]


COLUMNS = [
    "Taxon",
    "SpeciesName",
    "AlleleId",
    "AlleleSymbol",
    "AlleleSynonyms",
    "VariantId",
    "VariantSymbol",
    "VariantSynonyms",
    "VariantCrossReferences",
    "AlleleAssociatedGeneId",
    "AlleleAssociatedGeneSymbol",
    "VariantAffectedGeneId",
    "VariantAffectedGeneSymbol",
    "Category",
    "VariantsTypeId",
    "VariantsTypeName",
    "VariantsHgvsNames",
    "Assembly",
    "Chromosome",
    "StartPosition",
    "EndPosition",
    "SequenceOfReference",
    "SequenceOfVariant",
    "MostSevereConsequenceName",
    "VariantInformationReference",
    "HasDiseaseAnnotations",
    "HasPhenotypeAnnotations",
]


def _row(rec: dict) -> dict:
    return {
        "Taxon": normalize_taxon(rec.get("Taxon", "")),
        "SpeciesName": rec.get("SpeciesName", "") or "",
        "AlleleId": rec.get("AlleleId", "") or "",
        "AlleleSymbol": rec.get("AlleleSymbol", "") or "",
        "AlleleSynonyms": join_pipe(rec.get("AlleleSynonyms") or []),
        "VariantId": rec.get("VariantId", "") or "",
        "VariantSymbol": rec.get("VariantSymbol", "") or "",
        "VariantSynonyms": join_pipe(rec.get("VariantSynonyms") or []),
        "VariantCrossReferences": join_pipe(rec.get("VariantCrossReferences") or []),
        "AlleleAssociatedGeneId": (rec.get("AlleleAssociatedGeneId") or [""])[0] if rec.get("AlleleAssociatedGeneId") else "",
        "AlleleAssociatedGeneSymbol": (rec.get("AlleleAssociatedGeneSymbol") or [""])[0] if rec.get("AlleleAssociatedGeneSymbol") else "",
        "VariantAffectedGeneId": (rec.get("VariantAffectedGeneId") or [""])[0] if rec.get("VariantAffectedGeneId") else "",
        "VariantAffectedGeneSymbol": (rec.get("VariantAffectedGeneSymbol") or [""])[0] if rec.get("VariantAffectedGeneSymbol") else "",
        "Category": rec.get("Category", "") or "",
        "VariantsTypeId": rec.get("VariantsTypeId", "") or "",
        "VariantsTypeName": rec.get("VariantsTypeName", "") or "",
        "VariantsHgvsNames": rec.get("VariantsHgvsNames", "") or "",
        "Assembly": rec.get("Assembly", "") or "",
        "Chromosome": rec.get("Chromosome", "") or "",
        "StartPosition": str(rec.get("StartPosition") or ""),
        "EndPosition": str(rec.get("EndPosition") or ""),
        "SequenceOfReference": rec.get("SequenceOfReference", "") or "",
        "SequenceOfVariant": rec.get("SequenceOfVariant", "") or "",
        "MostSevereConsequenceName": join_pipe(rec.get("MostSevereConsequenceName") or []),
        "VariantInformationReference": join_pipe(rec.get("VariantInformationReference") or []),
        "HasDiseaseAnnotations": rec.get("HasDiseaseAnnotations", "") or "",
        "HasPhenotypeAnnotations": rec.get("HasPhenotypeAnnotations", "") or "",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance alleles into legacy 27-col TSV (FMS bulk)")
    parser.add_argument(
        "--taxa",
        default=",".join(DEFAULT_TAXA),
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

    taxa = [t.strip() for t in args.taxa.split(",") if t.strip()]
    release = get_current_release()
    out_dir = Path(args.out_dir)
    seen: set = set()
    total = 0

    with TsvWriter(
        out_dir / "alliance-alleles.tsv",
        COLUMNS,
        release=release,
        source=f"FMS VARIANT-ALLELE-JSON → legacy 27-col TSV ({','.join(taxa)})",
    ) as writer:
        for taxon in taxa:
            log.info("Alleles FMS: %s", taxon)
            for _sub, rec in iter_fms_datatype("VARIANT-ALLELE-JSON", sub_filter=taxon):
                allele_id = rec.get("AlleleId")
                if not allele_id:
                    continue
                # Multiple rows can share an AlleleId (one row per variant);
                # the converter looks at the allele-level identity for the
                # Allele item but emits one Variant item per row, so we
                # keep all rows (don't dedup on AlleleId alone).
                row = _row(rec)
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    log.info("Done. rows=%d", total)
                    return

    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
