#!/usr/bin/env python3
"""
Fetch chromosome lengths from NCBI Datasets sequence reports for the chromosomes
that alliance-genes.tsv places genes on.

alliance-genes creates Chromosome items with a name and no length. InterMine's
create-intergenic-region-features and create-gene-flanking-features skip a
chromosome whose length is null, so before this source only yeast (whose
lengths come from the sgd source) had intergenic and flanking regions.

Each taxon maps to the RefSeq assembly its MOD coordinates are on. A chromosome
the primary assembly lacks (GRCr8 and GRCz12tu have no MT) comes from the
fallback assembly. A chromosome whose longest gene ends past the NCBI length
is left out and logged: that means the assembly in ASSEMBLIES is stale.

Usage:
    python3 fetch_chromosome_lengths.py
    python3 fetch_chromosome_lengths.py --genes data/alliance-genes.tsv

Output: data/chromosome-lengths.tsv
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import TsvWriter, configure_logging, http_get_json  # noqa: E402

log = logging.getLogger("alliance.fetch.chromosome_lengths")

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"

# taxon -> (primary assembly, fallback assemblies). Checked 2026-10-26 against the
# gene coordinates in the rc3 alliance-genes.tsv: every gene ends inside its chromosome.
# Yeast is absent: the sgd source loads its chromosome sequences and lengths.
ASSEMBLIES: dict[str, tuple[str, tuple[str, ...]]] = {
    "9606": ("GCF_000001405.40", ()),                    # GRCh38.p14
    "10090": ("GCF_000001635.27", ()),                   # GRCm39
    "10116": ("GCF_036323735.1", ("GCF_015227675.2",)),  # GRCr8, MT from mRatBN7.2
    "7955": ("GCF_049306965.2", ("GCF_000002035.6",)),   # GRCz12tu, MT from GRCz11
    "7227": ("GCF_000001215.4", ()),                     # Release 6
    "6239": ("GCF_000002985.6", ()),                     # WBcel235
    "8364": ("GCF_000004195.4", ()),                     # X. tropicalis UCB_Xtro_10.0
    "8355": ("GCF_017654675.1", ()),                     # X. laevis Xenopus_laevis_v10.1
}

NAME_FIELDS = ("chr_name", "sequence_name", "ucsc_style_name", "refseq_accession", "genbank_accession")
MITO_NAMES = {"MT", "MtDNA", "mitochondrion_genome", "Mito", "chrmt"}
MITO_ALIASES = ("MT", "M", "chrM", "MtDNA", "mitochondrion", "Mito")

COLUMNS = ["taxonId", "chromosome", "length", "refseqAccession", "assembly"]


def sequence_reports(assembly: str) -> list[dict]:
    rows: list[dict] = []
    token = None
    while True:
        params = {"page_size": "1000"}
        if token:
            params["page_token"] = token
        data = http_get_json(f"{DATASETS}/{assembly}/sequence_reports?{urllib.parse.urlencode(params)}",
                             timeout=120)
        rows.extend(data.get("reports") or [])
        token = data.get("next_page_token")
        # NCBI allows 3 requests a second without an API key
        time.sleep(0.4)
        if not token:
            return rows


def index_reports(assembly: str, rows: list[dict]) -> dict[str, dict]:
    """Every name NCBI gives a sequence, as typed and lowercased, to its report row."""
    index: dict[str, dict] = {}
    for row in rows:
        row["assembly"] = assembly
        for field in NAME_FIELDS:
            name = row.get(field)
            if name and name != "na":
                index.setdefault(name, row)
                index.setdefault(name.lower(), row)
    return index


def lookup(index: dict[str, dict], name: str) -> dict | None:
    for cand in (name, name.lower(), "chr" + name, name.removeprefix("chr"), name.removeprefix("Chr")):
        if cand in index:
            return index[cand]
    if name in MITO_NAMES:
        for cand in MITO_ALIASES:
            if cand in index:
                return index[cand]
    return None


def gene_chromosome_ends(genes: Path) -> dict[str, dict[str, int]]:
    """taxon -> chromosome -> largest gene end, from alliance-genes.tsv."""
    ends: dict[str, dict[str, int]] = {}
    with open(genes, encoding="utf-8") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for rec in reader:
            chrom = (rec.get("chromosome") or "").strip()
            end = (rec.get("end") or "").strip()
            if chrom in ("", "null") or not end.isdigit():
                continue
            taxon = (rec.get("species") or "").removeprefix("NCBITaxon:")
            per_taxon = ends.setdefault(taxon, {})
            per_taxon[chrom] = max(per_taxon.get(chrom, 0), int(end))
    return ends


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch chromosome lengths from NCBI Datasets")
    default_data = Path(__file__).parent.parent / "data"
    parser.add_argument("--genes", help="alliance-genes.tsv whose chromosomes need lengths "
                        "(default: <out-dir>/alliance-genes.tsv)")
    parser.add_argument("--out-dir", default=str(default_data), help="Directory for the emitted TSV")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    ends = gene_chromosome_ends(Path(args.genes or Path(args.out_dir) / "alliance-genes.tsv"))
    written = missing = too_short = 0
    with TsvWriter(Path(args.out_dir) / "chromosome-lengths.tsv", COLUMNS,
                   source="NCBI Datasets sequence_reports") as writer:
        for taxon, (primary, fallbacks) in ASSEMBLIES.items():
            indexes = [index_reports(a, sequence_reports(a)) for a in (primary, *fallbacks)]
            for chrom, gene_end in sorted(ends.get(taxon, {}).items()):
                row = next((r for r in (lookup(i, chrom) for i in indexes) if r), None)
                if row is None:
                    log.warning("taxon %s chromosome %s: not in %s", taxon, chrom, primary)
                    missing += 1
                    continue
                if row["length"] < gene_end:
                    log.warning("taxon %s chromosome %s: gene ends at %d, %s length is %d; left out",
                                taxon, chrom, gene_end, row["assembly"], row["length"])
                    too_short += 1
                    continue
                writer.write_row(taxonId=taxon, chromosome=chrom, length=row["length"],
                                 refseqAccession=row.get("refseq_accession") or "",
                                 assembly=row["assembly"])
                written += 1
    log.info("Done. rows=%d  not_found=%d  shorter_than_genes=%d", written, missing, too_short)
    if too_short:
        sys.exit(f"{too_short} chromosomes are shorter than their genes; update ASSEMBLIES")


if __name__ == "__main__":
    main()
