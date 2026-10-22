"""Unit tests for scripts/disease_db.py and the extract_disease*.py row builders.
Run: python3 -m unittest discover -s scripts/tests"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import disease_db as ddb  # noqa: E402
import extract_disease as ed  # noqa: E402
import extract_disease_alleles as eda  # noqa: E402
import extract_disease_annotations as edn  # noqa: E402
import extract_disease_models as edm  # noqa: E402
import extract_experimental_disease as eed  # noqa: E402

GENE = ddb.Entity("MGI:1", "NCBITaxon:10090", "Abc", "gene")
GENE2 = ddb.Entity("MGI:2", "NCBITaxon:10090", "Def", "gene")
ALLELE = ddb.Entity("MGI:10", "NCBITaxon:10090", "Abc<tm1>", "allele")
AGM = ddb.Entity("MGI:100", "NCBITaxon:10090", "Abc<tm1>/Abc<tm1> [background:] B6", "agm")
ECO = {"ECO:0000033": ("author statement supported by traceable reference", "TAS"),
       "ECO:0000315": ("mutant phenotype evidence used in manual assertion", "IMP")}


def _ann(**over) -> ddb.Annotation:
    fields = dict(
        id="1", kind="agm", subject=AGM, disease_id="DOID:9", disease_name="disease nine",
        relation="is_model_of", negated=False, evidence_codes=["ECO:0000033"], reference="AGRKB:101",
        provider="MGI", secondary_provider="", date="20200101", qualifiers=[], with_genes=[],
        conditions=[], curie="AGRKB:100000000000002",
        links=[("inferred_gene", GENE), ("inferred_allele", ALLELE)],
    )
    fields.update(over)
    return ddb.Annotation(**fields)


class RelationTest(unittest.TestCase):
    def test_model_relations_become_implicated_on_genes(self) -> None:
        for rel in ("is_model_of", "is_ameliorated_model_of", "is_exacerbated_model_of"):
            self.assertEqual(ddb.gene_relation(rel, False), "is_implicated_in")
        self.assertEqual(ddb.gene_relation("is_model_of", True), "is_not_implicated_in")

    def test_negated_marker(self) -> None:
        self.assertEqual(ddb.gene_relation("is_marker_for", True), "is_not_marker_for")
        self.assertEqual(ddb.gene_relation("is_marker_for", False), "is_marker_for")

    def test_yyyymmdd(self) -> None:
        self.assertEqual(ddb.yyyymmdd("2012-01-06 18:16:47+00"), "20120106")
        self.assertEqual(ddb.yyyymmdd(""), "")

    def test_condition_label(self) -> None:
        self.assertEqual(ddb.condition_label("has_condition"), "Has Condition")


class AnnotationLinksTest(unittest.TestCase):
    def test_gene_annotation_is_its_own_gene(self) -> None:
        ann = _ann(kind="gene", subject=GENE, links=[])
        self.assertEqual(ann.genes(), [(ddb.SUBJECT, GENE)])
        self.assertEqual(ann.alleles(), [])

    def test_inferred_and_asserted_gene_listed_once(self) -> None:
        ann = _ann(links=[("asserted_gene", GENE), ("inferred_gene", GENE), ("asserted_gene", GENE2)])
        self.assertEqual([g.id for _, g in ann.genes()], ["MGI:1", "MGI:2"])
        self.assertEqual(ann.genes()[0][0], "inferred_gene")


class RefMapTest(unittest.TestCase):
    def test_pmid_first(self) -> None:
        self.assertEqual(ddb.pmid_first(["RGD:5", "PMID:7", "MGI:1"]), ["PMID:7", "MGI:1", "RGD:5"])

    def test_preferred_is_first_curie(self) -> None:
        refs = ddb.RefMap({"AGRKB:1": ["PMID:7", "MGI:5"]})
        self.assertEqual(refs.preferred("AGRKB:1"), "PMID:7")

    def test_unmapped_falls_back_to_agrkb(self) -> None:
        refs = ddb.RefMap()
        self.assertEqual(refs.preferred("AGRKB:9"), "AGRKB:9")
        self.assertEqual(refs.misses, {"AGRKB:9"})
        self.assertEqual(refs.curies(""), [])


class DetailRowsTest(unittest.TestCase):
    def test_groups_by_gene_disease_relation_qualifiers(self) -> None:
        anns = [
            _ann(id="1", curie="AGRKB:100000000000009", evidence_codes=["ECO:0000315"], reference="AGRKB:1"),
            _ann(id="2", curie="AGRKB:100000000000003", reference="AGRKB:2"),
            _ann(id="3", kind="gene", subject=GENE, links=[], relation="is_implicated_in",
                 qualifiers=["severity_of"], reference="AGRKB:3"),
        ]
        refs = ddb.RefMap({"AGRKB:1": ["PMID:1", "MGI:1"]})
        rows = edn.direct_rows(anns, {"DOID:9": ["DOID:4"]}, refs)
        self.assertEqual(len(rows), 2)
        plain = next(r for r in rows if not r["diseaseQualifiers"])
        self.assertEqual(plain["generatedRelationString"], "is_implicated_in")
        self.assertEqual(plain["evidenceCodes"], "ECO:0000033|ECO:0000315")
        self.assertEqual(plain["referenceIds"], "AGRKB:1|AGRKB:2")
        self.assertEqual(plain["pubmedIds"], "PMID:1|MGI:1|AGRKB:2")
        self.assertEqual(plain["parentSlimIds"], "DOID:9|DOID:4")
        self.assertEqual(plain["viaOrthologyOrder"], "0")
        self.assertEqual(plain["uniqueId"], "AGRKB:100000000000003")

    def test_taxon_filter(self) -> None:
        rows = edn.direct_rows([_ann()], {}, ddb.RefMap(), {"NCBITaxon:9606"})
        self.assertEqual(rows, [])

    def test_orthology_row(self) -> None:
        self.assertTrue(edn.is_orthology_row({"viaOrthologyOrder": "1"}))
        self.assertFalse(edn.is_orthology_row({"viaOrthologyOrder": "0"}))


class DiseaseAllianceRowsTest(unittest.TestCase):
    def test_agm_annotation_rows(self) -> None:
        ann = _ann(evidence_codes=["ECO:0000033", "ECO:0000315"],
                   conditions=[("has_condition", "standard conditions"), ("ameliorated_by", "diet:x")])
        rows = list(ed.rows_for(ann, ddb.RefMap({"AGRKB:101": ["PMID:7", "MGI:7"]}), ECO,
                                {"NCBITaxon:10090": "Mus musculus"}))
        self.assertEqual(len(rows), 6)  # (agm, gene, allele) x 2 ECO codes
        agm, gene, allele = rows[0], rows[2], rows[4]
        self.assertEqual(agm["DBobjectType"], "affected_genomic_model")
        self.assertEqual(agm["AssociationType"], "is_model_of")
        self.assertEqual(agm["InferredFromID"] if "InferredFromID" in agm else "", "")
        self.assertEqual(gene["DBobjectType"], "gene")
        self.assertEqual(gene["AssociationType"], "is_implicated_in")
        self.assertEqual(gene["InferredFromID"], "MGI:100")
        self.assertEqual(allele["DBobjectType"], "allele")
        self.assertEqual(agm["Reference"], "PMID:7")
        self.assertEqual(agm["ExperimentalCondition"], "Has Condition: standard conditions")
        self.assertEqual(agm["Modifier"], "Ameliorated By: diet:x")
        self.assertEqual(gene["ExperimentalCondition"], "")
        self.assertEqual(agm["SpeciesName"], "Mus musculus")
        self.assertEqual(rows[1]["EvidenceCodeName"], ECO["ECO:0000315"][0])

    def test_source_via_secondary_provider(self) -> None:
        self.assertEqual(ed.source_name(_ann(provider="OMIM", secondary_provider="RGD")), "OMIM Via RGD")

    def test_negated_subject_relation(self) -> None:
        ann = _ann(kind="gene", subject=GENE, links=[], relation="is_implicated_in", negated=True)
        rows = list(ed.rows_for(ann, ddb.RefMap(), ECO, {}))
        self.assertEqual(rows[0]["AssociationType"], "is_not_implicated_in")


class AlleleAndExperimentalRowsTest(unittest.TestCase):
    def test_allele_rows_merge_evidence(self) -> None:
        anns = [_ann(id="1", reference="AGRKB:1"),
                _ann(id="2", reference="AGRKB:2", evidence_codes=["ECO:0000315"])]
        rows = eda.allele_rows(anns, ddb.RefMap({"AGRKB:1": ["PMID:1"]}))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["alleleId"], "MGI:10")
        self.assertEqual(row["relationName"], "is_implicated_in")
        self.assertEqual(row["evidenceCodes"], "ECO:0000033|ECO:0000315")
        self.assertEqual(row["evidencePmids"], "PMID:1|AGRKB:2")
        self.assertEqual(row["evidenceCurie"], "AGRKB:1")
        self.assertEqual(row["negated"], "false")

    def test_experimental_rows_skip_direct_gene_annotations(self) -> None:
        anns = [_ann(), _ann(id="2", kind="gene", subject=GENE2, links=[], relation="is_implicated_in")]
        rows = eed.experimental_rows(anns, ddb.RefMap(), ECO)
        self.assertEqual([r["geneId"] for r in rows], ["MGI:1"])
        self.assertEqual(rows[0]["evidenceAbbrs"], "TAS")
        self.assertEqual(rows[0]["geneTaxon"], "NCBITaxon:10090")


class DiseaseModelsTest(unittest.TestCase):
    MODEL = {"id": "MGI:100", "name": "Abc<tm1>/Abc<tm1>", "subtype": "genotype", "provider": "MGI",
             "diseases": edm.disease_models("true=DOID:2=beta|false=DOID:1=Alpha"),
             "phenotypes": edm.sorted_phenotypes("small|Big"), "has_phenotypes": True}

    def test_disease_models(self) -> None:
        self.assertEqual(self.MODEL["diseases"], [
            {"diseaseId": "DOID:1", "diseaseName": "Alpha", "associationType": "IS_MODEL_OF"},
            {"diseaseId": "DOID:2", "diseaseName": "does not model beta", "associationType": "IS_NOT_MODEL_OF"},
        ])

    def test_rows_per_disease(self) -> None:
        rows = edm.model_rows("MGI:1", "NCBITaxon:10090", self.MODEL)
        self.assertEqual([r["diseaseId"] for r in rows], ["DOID:1", "DOID:2"])
        self.assertEqual(rows[0]["associatedPhenotypes"], "Big|small")
        self.assertEqual(rows[0]["hasDiseaseAnnotations"], "true")
        self.assertEqual(rows[0]["modifierRelationshipTypes"], "ameliorated_by|exacerbated_by")

    def test_model_without_disease_gives_one_row(self) -> None:
        model = dict(self.MODEL, diseases=[], phenotypes=[], has_phenotypes=False)
        rows = edm.model_rows("MGI:1", "NCBITaxon:10090", model)
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["diseaseId"], rows[0]["hasDiseaseAnnotations"],
                          rows[0]["hasPhenotypeAnnotations"]), ("", "false", "false"))


if __name__ == "__main__":
    unittest.main()
