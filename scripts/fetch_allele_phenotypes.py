#!/usr/bin/env python3
"""
Fetch per-allele phenotype annotations from the Alliance FMS PHENOTYPE
bulk files. Filters to rows that carry a non-empty primaryGeneticEntityIDs
list — those are the allele-keyed annotations (vs gene-direct orthology
annotations which lack the genetic-entity reference and are already
covered by alliance-phenotypes).

Emits one row per (gene, allele, term) triple so the integration engine
can merge against existing PhenotypeAnnotation rows by alleleSubject +
ontologyTerm.

Usage:
    python3 fetch_allele_phenotypes.py
    python3 fetch_allele_phenotypes.py --mods MGI,ZFIN --limit 100

Output: data/allele-phenotypes.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    enumerate_allele_ids,
    get_current_release,
    iter_fms_datatype,
    mod_gene_taxa,
    selected_mods,
)

log = logging.getLogger("alliance.fetch.allele_phenotypes")


DEFAULT_MODS = ["MGI", "ZFIN", "FB", "WB", "RGD", "SGD", "XBXL", "XBXT"]


COLUMNS = [
    "geneId",
    "alleleId",
    "termId",
    "phenotypeStatement",
    "evidencePmid",
    "dataProvider",
]


def _iter_mod(mod: str, genes: set[str], alleles: set[str]) -> Iterator[dict]:
    """objectId is the annotated entity: a gene for FB/SGD/WB, a gene or an
    allele for MGI/ZFIN. primaryGeneticEntityIDs mixes alleles with genotypes
    and fish (AGMs). Only ids in the release's allele list become alleleId."""
    for _sub, rec in iter_fms_datatype("PHENOTYPE", sub_filter=mod):
        object_id = rec.get("objectId") or ""
        entities = rec.get("primaryGeneticEntityIDs") or []
        if not object_id or not entities:
            continue
        terms = rec.get("phenotypeTermIdentifiers") or []
        if not terms:
            continue
        gene_id = object_id if object_id in genes else ""
        allele_ids = [e for e in entities if e in alleles]
        if object_id in alleles and object_id not in allele_ids:
            allele_ids.insert(0, object_id)
        statement = rec.get("phenotypeStatement") or ""
        pmid = (rec.get("evidence") or {}).get("publicationId", "")
        for entity_id in allele_ids:
            for t in terms:
                term_id = t.get("termId") if isinstance(t, dict) else ""
                if not term_id:
                    continue
                yield {
                    "geneId": gene_id,
                    "alleleId": entity_id,
                    "termId": term_id,
                    "phenotypeStatement": statement,
                    "evidencePmid": pmid,
                    "dataProvider": mod,
                }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance allele phenotypes (FMS bulk)")
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

    mods = selected_mods(DEFAULT_MODS, args.mods)
    release = get_current_release()
    out_dir = Path(args.out_dir)
    genes = set(mod_gene_taxa(mods))
    alleles = set(enumerate_allele_ids(mods))
    seen: set = set()
    total = 0

    with TsvWriter(
        out_dir / "allele-phenotypes.tsv",
        COLUMNS,
        release=release,
        source=f"FMS PHENOTYPE ({','.join(mods)})",
    ) as writer:
        for mod in mods:
            log.info("Allele phenotypes FMS: %s", mod)
            for row in _iter_mod(mod, genes, alleles):
                key = (row["alleleId"], row["termId"])
                if key in seen:
                    continue
                seen.add(key)
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    log.info("Done. rows=%d", total)
                    return

    log.info("Done. rows=%d", total)


if __name__ == "__main__":
    main()
