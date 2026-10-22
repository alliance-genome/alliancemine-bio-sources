import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import extract_interactions as ei  # noqa: E402


def record(subject, obj):
    return {
        "subjectId": subject, "subjectSymbol": subject.lower(), "subjectTaxon": "NCBITaxon:1",
        "objectId": obj, "objectSymbol": obj.lower(), "objectTaxon": "NCBITaxon:2",
        "interactorAType": "MI:0326", "interactorATypeName": "protein",
        "interactorBType": "MI:0319", "interactorBTypeName": "deoxyribonucleic acid",
        "interactorARole": "MI:0496", "interactorARoleName": "bait",
        "interactorBRole": "MI:0498", "interactorBRoleName": "prey",
        "uniqueId": "u1", "referenceId": "AGRKB:1",
        "pubmedId": "PMID:1",
    }


class QueriedGeneTest(unittest.TestCase):
    def test_smaller_seeded_id_wins(self):
        self.assertEqual(ei.queried_gene("B", "A", {"A", "B"}), "A")

    def test_larger_id_when_smaller_not_seeded(self):
        self.assertEqual(ei.queried_gene("A", "B", {"B"}), "B")

    def test_unseeded_pair_is_dropped(self):
        self.assertIsNone(ei.queried_gene("A", "B", {"C"}))

    def test_self_interaction(self):
        self.assertEqual(ei.queried_gene("A", "A", {"A"}), "A")


class TsvRowTest(unittest.TestCase):
    def test_object_side_swaps_genes_and_interactor_columns(self):
        row = ei.tsv_row(record("B", "A"), "A", "X (2000) J 1(2):3")
        self.assertEqual((row["geneId"], row["geneSymbol"], row["geneTaxon"]), ("A", "a", "NCBITaxon:2"))
        self.assertEqual((row["partnerGeneId"], row["partnerSymbol"]), ("B", "b"))
        self.assertEqual((row["interactorARole"], row["interactorARoleName"]), ("MI:0498", "prey"))
        self.assertEqual((row["interactorBRole"], row["interactorBRoleName"]), ("MI:0496", "bait"))
        self.assertEqual((row["interactorAType"], row["interactorBType"]), ("MI:0319", "MI:0326"))
        self.assertEqual((row["pubmedId"], row["shortCitation"]), ("PMID:1", "X (2000) J 1(2):3"))
        self.assertNotIn("subjectId", row)

    def test_subject_side_keeps_order(self):
        row = ei.tsv_row(record("A", "B"), "A", "")
        self.assertEqual((row["geneId"], row["partnerGeneId"]), ("A", "B"))
        self.assertEqual((row["interactorARole"], row["interactorBRole"]), ("MI:0496", "MI:0498"))


class SqlRowTest(unittest.TestCase):
    def test_genetic_row_maps_terms_and_phenotypes(self):
        values = [f"v{i}" for i in range(len(ei._BASE_FIELDS))]
        values += [f"t{i}" for i in range(12)] + ["colony size|growth"]
        rec = ei.sql_row_to_record("genetic", values)
        self.assertEqual((rec["uniqueId"], rec["pubmedId"]), ("v8", "v10"))
        self.assertEqual((rec["interactionType"], rec["interactionTypeName"]), ("t0", "t1"))
        self.assertEqual((rec["interactionSource"], rec["interactionSourceName"]), ("t10", "t11"))
        self.assertEqual(rec["phenotypesOrTraits"], "colony size|growth")

    def test_molecular_row_maps_detection_and_aggregation(self):
        values = [""] * len(ei._BASE_FIELDS) + [f"t{i}" for i in range(16)]
        rec = ei.sql_row_to_record("molecular", values)
        self.assertEqual((rec["detectionMethod"], rec["detectionMethodName"]), ("t12", "t13"))
        self.assertEqual((rec["aggregationDatabase"], rec["aggregationDatabaseName"]), ("t14", "t15"))


class ShortCitationTest(unittest.TestCase):
    def test_short_citation_drops_et_al_and_keeps_empty_parts(self):
        summary = {"short_citation": "Costanzo M et al. (2016) Science 353(6306)",
                   "resource_title": "Science", "volume": "353", "issue_name": "6306"}
        self.assertEqual(ei.short_citation(summary), "Costanzo M (2016) Science 353(6306):")

    def test_short_citation_with_pages(self):
        summary = {"short_citation": "Caudy AA et al. (2003) Nature 425(6956):411-4",
                   "resource_title": "Nature", "volume": "425", "issue_name": "6956",
                   "page_range": "411-4"}
        self.assertEqual(ei.short_citation(summary), "Caudy AA (2003) Nature 425(6956):411-4")

    def test_short_citation_page_text_repeats_volume(self):
        summary = {"short_citation": "Whitworth AJ et al. (2008) Dis Model Mech 1(2-3):168-74; discussion 173",
                   "resource_title": "Dis Model Mech", "volume": "1", "issue_name": "2-3",
                   "page_range": "168-74; discussion 173"}
        self.assertEqual(ei.short_citation(summary),
                         "Whitworth AJ (2008) Dis Model Mech 1(2-3):168-74; discussion 173")

    def test_short_citation_pages_without_volume(self):
        summary = {"short_citation": "van den Heuvel S (2005) WormBook 1-16",
                   "resource_title": "WormBook", "page_range": "1-16"}
        self.assertEqual(ei.short_citation(summary), "van den Heuvel S (2005) WormBook ():1-16")

    def test_short_citation_without_resource(self):
        summary = {"short_citation": "Lindsley D.L. et al. (1992) The Genome of Drosophila melanogaster."}
        self.assertEqual(ei.short_citation(summary), "Lindsley D.L. (1992)   ():")


if __name__ == "__main__":
    unittest.main()
