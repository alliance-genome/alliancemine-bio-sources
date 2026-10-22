"""Unit tests for scripts/extract_gene_crossrefs.py row mapping (no database)."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from extract_gene_crossrefs import row_from_db, xref_sql, xref_type, xref_url  # noqa: E402
from fetch_gene_crossrefs import COLUMNS  # noqa: E402


class XrefTypeTest(unittest.TestCase):
    def test_fms_spellings(self) -> None:
        self.assertEqual(xref_type("default", "NCBI_Gene"), "generic_cross_reference")
        self.assertEqual(xref_type("expression_atlas", "ENSEMBL"), "gene/expressionAtlas")
        self.assertEqual(xref_type("biogrid/orcs", "NCBI_Gene"), "gene/biogrid_orcs")
        self.assertEqual(xref_type("gene/references", "MGI"), "gene/references")
        self.assertEqual(xref_type("SGD", "PANTHER"), "gene/panther")
        self.assertEqual(xref_type("default", "PANTHER"), "generic_cross_reference")
        self.assertEqual(xref_type("", ""), "")


class XrefUrlTest(unittest.TestCase):
    def test_local_part_fills_template(self) -> None:
        self.assertEqual(
            xref_url("MGI:101757", "http://www.informatics.jax.org/marker/MGI:[%s]", "MGI", "MGI:101757"),
            "http://www.informatics.jax.org/marker/MGI:101757")

    def test_panther_tree_page(self) -> None:
        template = "http://pantherdb.org/treeViewer/treeViewer.jsp?book=PAN_BOOK&species=agr&seq=SGD=[%s]"
        self.assertEqual(
            xref_url("PANTHER:PTHR15180:SF1", template, "PANTHER", "SGD:S000000001"),
            "http://pantherdb.org/treeViewer/treeViewer.jsp?book=PTHR15180&species=agr&seq=SGD=S000000001")

    def test_template_without_placeholder_or_empty(self) -> None:
        self.assertEqual(xref_url("SGD:S1", "https://www.yeastgenome.org/", "SGD", "SGD:S1"),
                         "https://www.yeastgenome.org/")
        self.assertEqual(xref_url("X:1", "", "", "SGD:S1"), "")


class RowTest(unittest.TestCase):
    def test_row_has_fetcher_columns(self) -> None:
        row = row_from_db("MGI:101757", "NCBITaxon:10090", "UniProtKB:P18760",
                          "default", "http://www.uniprot.org/uniprot/[%s]", "UniProtKB")
        self.assertEqual(list(row), COLUMNS)
        self.assertEqual(row["xrefUrl"], "http://www.uniprot.org/uniprot/P18760")
        self.assertEqual(row["xrefType"], "generic_cross_reference")

    def test_panther_family_gets_mod_tree_page(self) -> None:
        pages = {"MGI": "http://pantherdb.org/treeViewer/treeViewer.jsp?book=PAN_BOOK&species=agr&seq=MGI=MGI=[%s]"}
        row = row_from_db("MGI:101757", "NCBITaxon:10090", "PANTHER:PTHR11913", "default",
                          "http://pantherdb.org/panther/family.do?clsAccession=[%s]", "PANTHER", pages)
        self.assertEqual(row["xrefType"], "gene/panther")
        self.assertEqual(
            row["xrefUrl"],
            "http://pantherdb.org/treeViewer/treeViewer.jsp?book=PTHR11913&species=agr&seq=MGI=MGI=101757")
        xenbase = row_from_db("Xenbase:XB-GENE-1", "NCBITaxon:8355", "PANTHER:PTHR1", "default",
                              "http://pantherdb.org/panther/family.do?clsAccession=[%s]", "PANTHER", pages)
        self.assertEqual(xenbase["xrefType"], "generic_cross_reference")

    def test_own_id_on_default_page_is_gene(self) -> None:
        row = row_from_db("HGNC:100", "NCBITaxon:9606", "HGNC:100", "default",
                          "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:[%s]", "HGNC")
        self.assertEqual(row["xrefType"], "gene")

    def test_sql_includes_gcrp_and_page_rank(self) -> None:
        sql = xref_sql(["NCBITaxon:559292"])
        self.assertIn("gcrpcrossreference_id", sql)
        self.assertIn("DISTINCT ON (g.primaryexternalid, x.referencedcurie)", sql)


if __name__ == "__main__":
    unittest.main()
