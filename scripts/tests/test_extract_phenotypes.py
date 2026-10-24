"""Unit tests for the phenotype and AGM curation-DB extractors. Run: python3 -m unittest discover -s scripts/tests"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import extract_allele_phenotypes as eap  # noqa: E402
import extract_phenotypes as ep  # noqa: E402


class RelationTest(unittest.TestCase):
    def test_gene_level_annotation_gives_has_phenotype(self) -> None:
        self.assertEqual(ep.relation_for(True), "has_phenotype")

    def test_allele_or_model_only_gives_is_implicated_in(self) -> None:
        self.assertEqual(ep.relation_for(False), "is_implicated_in")


class PublicationIdsTest(unittest.TestCase):
    REFS = {"AGRKB:1": ["PMID:2", "MGI:5"], "AGRKB:2": ["PMID:1"], "AGRKB:3": ["MGI:7"],
            "AGRKB:4": ["FB:FBrf1"]}

    def test_pmid_per_reference_distinct_and_sorted(self) -> None:
        self.assertEqual(ep.publication_ids(["AGRKB:1", "AGRKB:2", "AGRKB:1"], self.REFS, "MGI:1"),
                         ["PMID:1", "PMID:2"])

    def test_reference_without_pmid_gives_its_mod_id(self) -> None:
        self.assertEqual(ep.publication_ids(["AGRKB:3", "AGRKB:2"], self.REFS, "MGI:1"), ["MGI:7", "PMID:1"])

    def test_other_mod_ids_and_unknown_refs_are_dropped(self) -> None:
        self.assertEqual(ep.publication_ids(["AGRKB:4", "AGRKB:9"], self.REFS, "MGI:1"), [])

    def test_non_agrkb_evidence_is_its_own_publication(self) -> None:
        self.assertEqual(ep.publication_ids(["ORPHA:552", "MIM:1", "AGRKB:2"], self.REFS, "HGNC:5"),
                         ["MIM:1", "ORPHA:552", "PMID:1"])


class ToRowTest(unittest.TestCase):
    def test_row_matches_fetcher_columns(self) -> None:
        group = {"gene": "SGD:S1", "taxon": "NCBITaxon:559292", "statement": "inviable",
                 "direct": False, "refs": {"AGRKB:2", "AGRKB:1"}}
        row = ep.to_row(group, {"AGRKB:1": ["PMID:7"]})
        self.assertEqual(list(row), ep.COLUMNS)
        self.assertEqual(row["uniqueId"], "SGD:S1|inviable")
        self.assertEqual(row["referenceIds"], "AGRKB:1|AGRKB:2")
        self.assertEqual(row["pubmedIds"], "PMID:7")
        self.assertEqual(row["relation"], "is_implicated_in")


class GeneFilterTest(unittest.TestCase):
    def test_small_sets_become_a_values_join(self) -> None:
        sql = ep.gene_filter_sql(["MGI:1", "ZFIN:x'y"])
        self.assertIn("('MGI:1')", sql)
        self.assertIn("('ZFIN:x''y')", sql)

    def test_no_or_large_sets_need_no_join(self) -> None:
        self.assertEqual(ep.gene_filter_sql(None), "")
        self.assertEqual(ep.gene_filter_sql([f"G:{i}" for i in range(ep.SQL_ID_FILTER_MAX + 1)]), "")


class ModForTest(unittest.TestCase):
    def test_xenbase_splits_by_taxon(self) -> None:
        self.assertEqual(eap.mod_for("XB", "NCBITaxon:8355"), "XBXL")
        self.assertEqual(eap.mod_for("XB", "NCBITaxon:8364"), "XBXT")
        self.assertEqual(eap.mod_for("XB", "NCBITaxon:1"), "")

    def test_other_providers_are_mod_names(self) -> None:
        self.assertEqual(eap.mod_for("MGI", "NCBITaxon:10090"), "MGI")

    def test_providers_sql_keeps_mod_order_and_merges_xenbase(self) -> None:
        self.assertEqual(eap.providers_sql(["ZFIN", "XBXL", "XBXT", "MGI"]), "'ZFIN','XB','MGI'")


class SplitModsTest(unittest.TestCase):
    def test_termless_mods_stay_on_fms(self) -> None:
        db, fms = eap.split_mods(["MGI", "ZFIN", "FB", "WB", "RGD", "SGD"])
        self.assertEqual(db, ["MGI", "RGD", "SGD"])
        self.assertEqual(fms, ["ZFIN", "FB", "WB"])


class EvidencePublicationTest(unittest.TestCase):
    def test_pmid_wins(self) -> None:
        self.assertEqual(eap.evidence_publication(["MGI:9", "PMID:3"], "MGI"), "PMID:3")

    def test_mod_id_when_no_pmid(self) -> None:
        self.assertEqual(eap.evidence_publication(["FB:FBrf2", "MGI:9"], "FB"), "FB:FBrf2")
        self.assertEqual(eap.evidence_publication(["Xenbase:XB-ART-1"], "XBXT"), "Xenbase:XB-ART-1")

    def test_nothing_known(self) -> None:
        self.assertEqual(eap.evidence_publication([], "WB"), "")


if __name__ == "__main__":
    unittest.main()
