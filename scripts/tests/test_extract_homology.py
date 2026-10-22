"""Unit tests for the pure row functions of extract_orthologs.py and extract_paralogs.py."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import extract_orthologs  # noqa: E402
import extract_paralogs  # noqa: E402
import fetch_orthologs  # noqa: E402
import fetch_paralogs  # noqa: E402

HUMAN = "NCBITaxon:9606"
MOUSE = "NCBITaxon:10090"


def ortholog(subject, obj, matched="PANTHER"):
    return extract_orthologs.to_row([subject, MOUSE, obj, HUMAN, matched, "OMA", "ZFIN"])


def paralog(subject, obj, rank="1", identity="50"):
    return extract_paralogs.to_row(
        [subject, "a", MOUSE, obj, "b", MOUSE, identity, "60", "300", rank, "PANTHER", "", "HGNC"])


class OrthologRowTest(unittest.TestCase):
    def test_row_has_fetcher_columns_and_stringent_label(self):
        row = ortholog("MGI:1", "HGNC:2", "Ensembl Compara|PANTHER")
        self.assertEqual(list(row), fetch_orthologs.COLUMNS)
        self.assertEqual(row["stringencyFilter"], "stringent")
        self.assertEqual(row["predictionMethodsMatched"], "Ensembl Compara|PANTHER")
        self.assertEqual(row["orthologGeneTaxon"], HUMAN)

    def test_keeps_seed_subjects_in_seed_order_once_per_pair(self):
        rows = [ortholog("MGI:2", "HGNC:9"), ortholog("HGNC:9", "MGI:2"),
                ortholog("MGI:1", "HGNC:8"), ortholog("MGI:1", "HGNC:8", "OMA"),
                ortholog("", "HGNC:7")]
        out = extract_orthologs.seeded_rows(rows, {"MGI:1": 0, "MGI:2": 1})
        self.assertEqual([(r["subjectGeneId"], r["orthologGeneId"]) for r in out],
                         [("MGI:1", "HGNC:8"), ("MGI:2", "HGNC:9")])
        self.assertEqual(out[0]["predictionMethodsMatched"], "PANTHER")

    def test_mirror_swaps_genes_and_keeps_methods(self):
        row = extract_orthologs.mirror(ortholog("MGI:1", "HGNC:2", "OMA"))
        self.assertEqual((row["subjectGeneId"], row["subjectGeneTaxon"]), ("HGNC:2", HUMAN))
        self.assertEqual((row["orthologGeneId"], row["orthologGeneTaxon"]), ("MGI:1", MOUSE))
        self.assertEqual(row["predictionMethodsMatched"], "OMA")

    def test_gap_seeds_are_known_seeds_of_uncovered_taxa(self):
        taxa = {"SGD:1": "NCBITaxon:559292", "MGI:1": MOUSE}
        self.assertEqual(extract_orthologs.gap_seeds(["SGD:1", "MGI:1", "X:1"], taxa, {MOUSE}), ["SGD:1"])


class ParalogRowTest(unittest.TestCase):
    def test_row_has_fetcher_columns(self):
        row = paralog("MGI:1", "MGI:2")
        self.assertEqual(list(row), fetch_paralogs.COLUMNS)

    def test_zero_metric_is_empty_like_the_fetcher(self):
        row = paralog("MGI:1", "MGI:2", identity="0")
        self.assertEqual(row["identity"], "")
        self.assertEqual(row["rank"], "1")

    def test_pair_keeps_direction_of_earliest_seed(self):
        rows = [paralog("MGI:2", "MGI:1", rank="3"), paralog("MGI:1", "MGI:2", rank="1"),
                paralog("MGI:2", "MGI:9", rank="5")]
        out = extract_paralogs.seeded_rows(rows, {"MGI:1": 0, "MGI:2": 1})
        self.assertEqual([(r["subjectGeneId"], r["paralogueGeneId"], r["rank"]) for r in out],
                         [("MGI:1", "MGI:2", "1"), ("MGI:2", "MGI:9", "5")])

    def test_non_seed_subject_is_dropped(self):
        out = extract_paralogs.seeded_rows([paralog("MGI:9", "MGI:2")], {"MGI:2": 0})
        self.assertEqual(out, [])


if __name__ == "__main__":
    unittest.main()
