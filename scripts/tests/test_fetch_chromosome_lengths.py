"""Unit tests for scripts/fetch_chromosome_lengths.py name matching (no network)."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fetch_chromosome_lengths import gene_chromosome_ends, index_reports, lookup  # noqa: E402

REPORTS = [
    {"chr_name": "1", "sequence_name": "chr1", "refseq_accession": "NC_086019.1", "length": 270518180},
    {"chr_name": "1L", "sequence_name": "Chr1L", "refseq_accession": "NC_054371.1", "length": 233740091},
    {"chr_name": "Un", "sequence_name": "211000022278279", "refseq_accession": "NW_007931104.1",
     "genbank_accession": "na", "length": 12714},
    {"chr_name": "MT", "sequence_name": "MT", "refseq_accession": "NC_001665.2", "length": 16313},
]


class LookupTest(unittest.TestCase):
    def setUp(self) -> None:
        self.index = index_reports("GCF_TEST.1", [dict(r) for r in REPORTS])

    def test_plain_and_prefixed_names(self) -> None:
        self.assertEqual(lookup(self.index, "1")["length"], 270518180)
        self.assertEqual(lookup(self.index, "Chr1L")["length"], 233740091)
        self.assertEqual(lookup(self.index, "chr1L")["length"], 233740091)

    def test_scaffold_by_sequence_name(self) -> None:
        self.assertEqual(lookup(self.index, "211000022278279")["refseq_accession"], "NW_007931104.1")

    def test_mitochondrion_aliases(self) -> None:
        for name in ("MtDNA", "mitochondrion_genome", "chrmt"):
            self.assertEqual(lookup(self.index, name)["length"], 16313, name)

    def test_na_is_not_a_name(self) -> None:
        self.assertIsNone(lookup(self.index, "na"))

    def test_assembly_recorded(self) -> None:
        self.assertEqual(lookup(self.index, "1")["assembly"], "GCF_TEST.1")


class GeneEndsTest(unittest.TestCase):
    def test_largest_end_per_chromosome(self) -> None:
        rows = [
            "# header comment",
            "primaryIdentifier\tspecies\tchromosome\tstart\tend",
            "RGD:1\tNCBITaxon:10116\t1\t10\t500",
            "RGD:2\tNCBITaxon:10116\t1\t20\t900",
            "RGD:3\tNCBITaxon:10116\tnull\tnull\tnull",
            "FB:1\tNCBITaxon:7227\t2L\t1\t40",
        ]
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("\n".join(rows) + "\n")
        try:
            self.assertEqual(gene_chromosome_ends(Path(fh.name)),
                             {"10116": {"1": 900}, "7227": {"2L": 40}})
        finally:
            Path(fh.name).unlink()


if __name__ == "__main__":
    unittest.main()
