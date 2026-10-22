"""Unit tests for scripts/extract_genes.py row mapping (no database)."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import extract_genes  # noqa: E402
from extract_genes import SEP, row_from_db  # noqa: E402
from fetch_genes import COLUMNS, _row_for_gene  # noqa: E402


def _db_rec(**overrides) -> dict:
    rec = {field: "" for field in extract_genes.DB_FIELDS}
    rec.update({
        "primaryExternalId": "SGD:S000000001",
        "taxon": "NCBITaxon:559292",
        "dataProvider": "SGD",
        "geneType": "protein_coding_gene",
        "dateUpdated": "2000-05-19T00:05:00.000Z",
        "symbol": "TFC3",
        "fullName": "TFC3",
        "systematicName": "YAL001C",
        "synonyms": SEP.join(["FUN24", "TSV115", "tau 138"]),
        "secondaryIds": SEP.join(["SGD:L000000641", "SGD:L000002287"]),
        "crossrefs": SEP.join(["SGD:YAL001C", "SGD:S000000001", "SGD:S000000001", "DIP:6739"]),
        "modDescription": "Subunit of RNA polymerase III",
        "chromosome": "chrI",
        "start": "147594",
        "end": "151166",
        "strand": "-",
        "modXrefCurie": "SGD:S000000001",
        "modXrefTemplate": "https://www.yeastgenome.org/locus/[%s]",
    })
    rec.update(overrides)
    return rec


def _api_payload() -> dict:
    """The /gene/{id} payload matching _db_rec(), as the API serves it."""
    return {"gene": {
        "primaryExternalId": "SGD:S000000001",
        "taxon": {"curie": "NCBITaxon:559292"},
        "dataProvider": {"abbreviation": "SGD"},
        "geneType": {"name": "protein_coding_gene"},
        "geneSymbol": {"displayText": "TFC3"},
        "geneFullName": {"displayText": "TFC3"},
        "geneSystematicName": {"displayText": "YAL001C"},
        "geneSynonyms": [{"displayText": s} for s in ("FUN24", "TSV115", "tau 138")],
        "geneSecondaryIds": [{"secondaryId": "SGD:L000000641"}, {"secondaryId": "SGD:L000002287"}],
        "crossReferences": [{"referencedCurie": c}
                            for c in ("SGD:YAL001C", "SGD:S000000001", "SGD:S000000001", "DIP:6739")],
        "relatedNotes": [
            {"freeText": "Subunit of RNA polymerase III",
             "noteType": {"name": "MOD_provided_gene_description"}},
            {"freeText": "Contributes to DNA binding", "noteType": {"name": "automated_gene_description"}},
        ],
        "geneGenomicLocationAssociations": [{
            "start": 147594, "end": 151166, "strand": "-",
            "geneGenomicLocationAssociationObject": {"name": "chrI"},
        }],
        "dataProviderCrossReference": {
            "referencedCurie": "SGD:S000000001",
            "resourceDescriptorPage": {"urlTemplate": "https://www.yeastgenome.org/locus/[%s]"},
        },
    }}


class RowFromDbTest(unittest.TestCase):
    def test_matches_api_mapping_except_date(self) -> None:
        auto = {"SGD:S000000001": "Contributes to DNA binding"}
        db_row = row_from_db(_db_rec(), auto)
        api_row = _row_for_gene(_api_payload())
        self.assertEqual(set(db_row), set(COLUMNS))
        for col in COLUMNS:
            if col != "dateProduced":
                self.assertEqual(db_row[col], api_row[col], col)
        self.assertEqual(db_row["dateProduced"], "2000-05-19T00:05:00.000Z")

    def test_list_columns(self) -> None:
        row = row_from_db(_db_rec(), {})
        self.assertEqual(row["secondaryIdentifier"], "SGD:L000000641")
        self.assertEqual(row["synonyms"], "[FUN24, TSV115, tau 138, YAL001C, SGD:L000002287]")
        self.assertEqual(row["crossrefs"], "[SGD:YAL001C, SGD:S000000001, DIP:6739]")
        self.assertEqual(row["modCrossRefCompleteUrl"], "https://www.yeastgenome.org/locus/S000000001")

    def test_missing_location_and_names(self) -> None:
        row = row_from_db(_db_rec(chromosome="", start="", end="", strand="", synonyms="",
                                  secondaryIds="", systematicName="", crossrefs="",
                                  dateUpdated="", dateCreated="2006-05-15T01:38:14.985Z"), {})
        self.assertEqual([row[c] for c in ("chromosome", "start", "end", "strand")], ["null"] * 4)
        self.assertEqual(row["synonyms"], "[]")
        self.assertEqual(row["crossrefs"], "[]")
        self.assertEqual(row["secondaryIdentifier"], "")
        self.assertEqual(row["automatedDescription"], "")
        self.assertEqual(row["dateProduced"], "2006-05-15T01:38:14.985Z")

    def test_mod_url_needs_placeholder(self) -> None:
        row = row_from_db(_db_rec(modXrefTemplate="https://www.yeastgenome.org/"), {})
        self.assertEqual(row["modCrossRefCompleteUrl"], "")


class SqlTest(unittest.TestCase):
    def test_ids_are_quoted(self) -> None:
        sql = extract_genes.gene_sql(["NCBITaxon:9606"], ["HGNC:1", "X:o'k"])
        self.assertIn("'NCBITaxon:9606'", sql)
        self.assertIn("'X:o''k'", sql)
        self.assertIn("AND b.primaryexternalid IN", sql)

    def test_mods_for_ids(self) -> None:
        self.assertEqual(extract_genes._mods_for_ids(["Xenbase:XB-GENE-1", "HGNC:5"]),
                         ["XBXL", "XBXT", "HUMAN"])


if __name__ == "__main__":
    unittest.main()
