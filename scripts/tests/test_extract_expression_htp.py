"""Unit tests for the row mapping in extract_expression.py and extract_htp.py.
Run: python3 -m unittest discover -s scripts/tests"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import extract_expression as ee  # noqa: E402
import extract_htp as eh  # noqa: E402
from fetch_htp import DATASET_COLUMNS, SAMPLE_COLUMNS  # noqa: E402


class PickReferenceTest(unittest.TestCase):
    def test_pmid_wins(self) -> None:
        self.assertEqual(ee.pick_reference(["FB:FBrf0219073", "PMID:7600985"], "FB"), "PMID:7600985")

    def test_own_mod_id_without_pmid(self) -> None:
        self.assertEqual(ee.pick_reference(["MGI:1", "ZFIN:ZDB-PUB-010810-1"], "ZFIN"), "ZFIN:ZDB-PUB-010810-1")
        self.assertEqual(ee.pick_reference(["Xenbase:XB-ART-1"], "XBXT"), "Xenbase:XB-ART-1")

    def test_other_mod_id_then_empty(self) -> None:
        self.assertEqual(ee.pick_reference(["WB:WBPaper00068304"], "MGI"), "WB:WBPaper00068304")
        self.assertEqual(ee.pick_reference([], "MGI"), "")


class ExpressionRowTest(unittest.TestCase):
    def test_build_row(self) -> None:
        middle = ["caudal fin proximal region", "Adult", "MMO:0000669", "western blot assay",
                  "", "", "", "", "", "", "", "",
                  "ZFA:0001058", "caudal fin", "BSPO:0000077", "proximal region"]
        rec = ["NCBITaxon:7955", "ZFIN:ZDB-GENE-100519-4", "aldh1l1", *middle,
               "https://zfin.org/ZDB-FIG-170119-12|https://zfin.org/ZDB-FIG-170119-11", "11|10|12"]
        refs = {"10": ["PMID:28096348", "ZFIN:ZDB-PUB-1"], "11": ["PMID:28096348"], "12": ["ZFIN:ZDB-PUB-2"]}
        row = ee.build_row(rec, "ZFIN", refs)
        self.assertEqual(len(row), len(ee.COLUMNS))
        self.assertEqual(row[:4], ["Danio rerio", "NCBITaxon:7955", "ZFIN:ZDB-GENE-100519-4", "aldh1l1"])
        self.assertEqual(row[4:20], middle)
        self.assertEqual(row[20], "https://zfin.org/ZDB-FIG-170119-11,https://zfin.org/ZDB-FIG-170119-12")
        self.assertEqual(row[21], "ZFIN")
        self.assertEqual(row[22], "PMID:28096348,ZFIN:ZDB-PUB-2")

    def test_no_xrefs_or_refs(self) -> None:
        rec = ["NCBITaxon:10116", "RGD:2525", "Dspp", *([""] * 16), "", ""]
        row = ee.build_row(rec, "RGD", {})
        self.assertEqual(row[20:], ["", "RGD", ""])

    def test_header_block_is_fourteen_lines(self) -> None:
        # AllianceExpressionConverter skips the first 15 lines: 14 comments plus the header.
        lines = ee.header_lines(["MGI", "XBXL"], "9.2.0")
        self.assertEqual(len(lines), 14)
        self.assertTrue(all(line.startswith("#") for line in lines))
        self.assertIn("NCBITaxon:10090, NCBITaxon:8355", lines[8])


class HtpRowTest(unittest.TestCase):
    def test_sex_label(self) -> None:
        self.assertEqual(eh.sex_label("unknown sex"), "unknown")
        self.assertEqual(eh.sex_label("pooled sexes"), "pooled")
        self.assertEqual(eh.sex_label("hermaphrodite"), "hermaphrodite")
        self.assertEqual(eh.sex_label(""), "")

    def test_dataset_row(self) -> None:
        rec = ["ZFIN:ZDB-HTPDSET-200812-4", "RNA-seq", "Summary", "baseline|developmental time course",
               "GEO:GSE45706", "7|5|9", "", "ZFIN"]
        refs = {"5": ["PMID:2"], "7": ["ZFIN:ZDB-PUB-1"], "9": ["PMID:2"]}
        row = eh.dataset_row(rec, refs, "ZFIN")
        self.assertEqual(list(row), DATASET_COLUMNS)
        self.assertEqual(row["publicationPmids"], "PMID:2|ZFIN:ZDB-PUB-1")
        self.assertEqual(row["geoXref"], "GEO:GSE45706")
        self.assertEqual(row["dataProvider"], "ZFIN")

    def test_sample_row(self) -> None:
        rec = ["GEO:GSM1112155", "64-cell, 2 hpf", "~2000 individuals", "OBI:0000869", "ZFS:0000007",
               "Cleavage:64-cell", "ZFA:0001094", "whole organism", "ZFIN:ZDB-FISH-150901-27842",
               "unknown sex", "MMO:0000659", "Zv9", "ZFIN:ZDB-HTPDSET-200812-4", "NCBITaxon:7955"]
        row = eh.sample_row(rec, "ZFIN")
        self.assertEqual(list(row), SAMPLE_COLUMNS)
        self.assertEqual(row["sex"], "unknown")
        self.assertEqual(row["anatomyStatement"], "whole organism")
        self.assertEqual(row["dataProvider"], "ZFIN")

    def test_sample_row_legacy_yeast_taxon(self) -> None:
        rec = ["GEO:GSM1", *([""] * 12), "NCBITaxon:4932"]
        self.assertEqual(eh.sample_row(rec, "SGD")["taxon"], "NCBITaxon:559292")


if __name__ == "__main__":
    unittest.main()
