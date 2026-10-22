"""Unit tests for the pure mapping functions of the allele extractors."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import extract_alleles  # noqa: E402
import extract_alleles_common as eac  # noqa: E402
import extract_allele_detail  # noqa: E402
import extract_transgenic_alleles  # noqa: E402
import extract_variants  # noqa: E402


class CountLabels(unittest.TestCase):
    def test_alteration_type_buckets(self):
        self.assertEqual(eac.alteration_type(0), "allele")
        self.assertEqual(eac.alteration_type(1), "allele with one variant")
        self.assertEqual(eac.alteration_type(5), "allele with multiple variants")

    def test_fms_category_counts(self):
        self.assertEqual(eac.fms_category(0), "allele")
        self.assertEqual(eac.fms_category(1), "allele with 1 known variant")
        self.assertEqual(eac.fms_category(3), "allele with 3 known variants")

    def test_yes_dash(self):
        self.assertEqual(eac.yes_dash(True), "yes")
        self.assertEqual(eac.yes_dash(False), "-")


class VariantFormats(unittest.TestCase):
    def test_variant_symbol_swaps_accession_for_chromosome(self):
        self.assertEqual(eac.variant_symbol("2L", "NT_033779.5:g.18844476_18844477ins"),
                         "2L:18844476_18844477ins")
        self.assertEqual(eac.variant_symbol("23", "NC_007134.7:g.28979215G>A"), "23:28979215G>A")

    def test_variant_symbol_without_chromosome_keeps_hgvs(self):
        self.assertEqual(eac.variant_symbol("", "NC_1.1:g.5A>T"), "NC_1.1:g.5A>T")

    def test_padded_sequences(self):
        self.assertEqual(eac.padded_sequences("deletion", "TC", "", "C"), ("CTC", "C"))
        self.assertEqual(eac.padded_sequences("insertion", "", "", "T"), ("T", "T."))
        self.assertEqual(eac.padded_sequences("insertion", "", "GT", "A"), ("A", "AGT"))
        self.assertEqual(eac.padded_sequences("point_mutation", "G", "A", ""), ("G", "A"))

    def test_severity_joined(self):
        self.assertEqual(eac.severity_joined([("5_prime_UTR_variant", 250), ("splice_region_variant", 160)]),
                         "splice_region_variant,5_prime_UTR_variant")
        self.assertEqual(eac.severity_joined([("x", None), ("y", 9), ("y", 9)]), "y,x")
        self.assertEqual(eac.severity_joined([("", None)]), "")

    def test_consequences_grouped_per_location(self):
        rows = [
            ["L1", "P1", "intron_variant", "20", "FB:FBgn1", "a"],
            ["L1", "P1", "splice_region_variant", "12", "FB:FBgn1", "a"],
            ["L1", "P2", "intron_variant", "20", "FB:FBgn2", "b"],
            ["L2", "P3", "", "", "", ""],
        ]
        out = eac._consequences_by_location(rows)
        self.assertEqual(out["L1"]["consequences"], ["splice_region_variant,intron_variant", "intron_variant"])
        self.assertEqual(out["L1"]["genes"], [("FB:FBgn1", "a"), ("FB:FBgn2", "b")])
        self.assertEqual(out["L2"], {"consequences": [], "genes": []})


class References(unittest.TestCase):
    def test_pmid_preferred_then_mod_then_any_paper(self):
        self.assertEqual(eac.reference_curie(["FB:FBrf1", "PMID:9", "DOI:x"], "FB"), "PMID:9")
        self.assertEqual(eac.reference_curie(["MGI:2", "FB:FBrf1"], "FB"), "FB:FBrf1")
        self.assertEqual(eac.reference_curie(["MGI:2", "DOI:x"], "FB"), "MGI:2")
        self.assertEqual(eac.reference_curie(["DOI:x", "PMCID:PMC1"], "FB"), "")

    def test_information_references_dedup_in_order(self):
        refs = [["PMID:1"], ["FB:FBrf2"], ["PMID:1", "FB:FBrf3"], ["DOI:x"]]
        self.assertEqual(eac.information_references(refs, "FB"), ["PMID:1", "FB:FBrf2"])

    def test_references_by_variant(self):
        rows = [["V1", "R1", "PMID:1"], ["V1", "R1", "FB:a"], ["V1", "R2", "FB:b"]]
        self.assertEqual(eac._references_by_variant(rows), {"V1": [["PMID:1", "FB:a"], ["FB:b"]]})


class Scope(unittest.TestCase):
    def test_mod_for_provider_and_xenopus_taxon(self):
        self.assertEqual(eac.mod_for("FB", "NCBITaxon:7227"), "FB")
        self.assertEqual(eac.mod_for("XB", "NCBITaxon:8364"), "XBXT")
        self.assertEqual(eac.mod_for("Alliance", "NCBITaxon:9606"), "HUMAN")

    def test_provider_filter(self):
        self.assertEqual(eac.provider_filter(["XBXL", "XBXT", "FB", "HUMAN"]), "'FB', 'XB'")
        self.assertEqual(eac.provider_filter(["HUMAN"]), "''")

    def test_split_agg(self):
        self.assertEqual(eac.split_agg(""), [])
        self.assertEqual(eac.split_agg("a\x1fb"), ["a", "b"])


VARIANT = {
    "hgvs": "NT_033779.5:g.18844476_18844477ins", "symbol": "2L:18844476_18844477ins",
    "type_id": "SO:0000667", "type_name": "insertion", "assembly": "R6", "chromosome": "2L",
    "start": "18844476", "end": "18844477", "ref": "T", "var": "T.",
    "consequences": ["intron_variant", "intron_variant"], "gene_ids": ["FB:FBgn1", "FB:FBgn2"],
    "gene_symbols": ["a", "b"], "xrefs": ["FB:FBal0218094"], "references": [["PMID:1"]],
}


class Rows(unittest.TestCase):
    def test_species_names_fall_back_to_provider(self):
        by_taxon, by_provider = extract_alleles.species_names(
            [["NCBITaxon:559292", "SGD", "Saccharomyces cerevisiae"], ["NCBITaxon:8355", "XB", "Xenopus laevis"],
             ["NCBITaxon:8364", "XB", "Xenopus tropicalis"]])
        self.assertEqual(by_taxon["NCBITaxon:8364"], "Xenopus tropicalis")
        self.assertEqual(by_provider["XB"], "Xenopus laevis")

    def test_detail_row(self):
        self.assertEqual(extract_allele_detail.detail_row("SGD:S1", 2, "SGD:S1"), {
            "alleleId": "SGD:S1", "alterationType": "allele with multiple variants",
            "apiCategory": "allele_summary", "apiCrossReference": "SGD:S1"})

    def test_allele_fields(self):
        row = extract_alleles.allele_fields("FB:FBal1", "NCBITaxon:7227", "Drosophila melanogaster",
                                            "w<sup>1</sup>", ["b", "a"], "FB:FBgn1", "w", False, True, 1)
        self.assertEqual(row["AlleleSynonyms"], "a|b")
        self.assertEqual(row["Category"], "allele with 1 known variant")
        self.assertEqual((row["HasDiseaseAnnotations"], row["HasPhenotypeAnnotations"]), ("-", "yes"))
        self.assertTrue(set(row) <= set(extract_alleles.COLUMNS))

    def test_allele_variant_fields(self):
        row = extract_alleles.variant_fields(VARIANT, "FB")
        self.assertEqual(row["VariantAffectedGeneId"], "FB:FBgn1")
        self.assertEqual(row["MostSevereConsequenceName"], "intron_variant|intron_variant")
        self.assertEqual(row["VariantInformationReference"], "PMID:1")
        self.assertTrue(set(row) <= set(extract_alleles.COLUMNS))

    def test_variant_row(self):
        row = extract_variants.variant_row(VARIANT, "FB:FBal0218094", "NCBITaxon:7227")
        self.assertEqual(row["geneId"], "FB:FBgn1|FB:FBgn2")
        self.assertEqual(row["dataProvider"], "FB")
        self.assertEqual(list(row), extract_variants.COLUMNS)

    def test_variant_order_follows_fetcher_taxa(self):
        self.assertEqual(extract_variants.MOD_ORDER[:3], ["FB", "ZFIN", "WB"])

    def test_transgenic_row(self):
        row = extract_transgenic_alleles.transgenic_row(
            "MGI:1", "MGI", "NCBITaxon:10090", "Tg(x)1", ["MGI:1_con", "MGI:1_con", "MGI:2_con"], "")
        self.assertEqual(row["constructIds"], "MGI:1_con|MGI:2_con")
        self.assertEqual(list(row), extract_transgenic_alleles.COLUMNS)


if __name__ == "__main__":
    unittest.main()
