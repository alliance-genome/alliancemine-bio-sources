"""Unit tests for scripts/fetch_disease_models.py. Run: python3 -m unittest discover -s scripts/tests"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import fetch_disease_models as fdm  # noqa: E402

OBO = """[Term]
id: DOID:1
name: alpha disease
synonym: "beta disease" EXACT []

[Term]
id: DOID:2
name: beta disease
synonym: "gamma disease" EXACT []
synonym: "delta disease" RELATED []
"""


def _row(**over) -> dict:
    row = {"modelId": "M1", "diseaseId": "DOID:1", "diseaseName": "alpha disease",
           "subjectGeneId": "G1", "associationType": "is_model_of"}
    row.update(over)
    return row


class LoadDoNamesTest(unittest.TestCase):
    def _load(self) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "do.obo"
            path.write_text(OBO, encoding="utf-8")
            return fdm.load_do_names(path)

    def test_primary_name_beats_earlier_synonym(self) -> None:
        self.assertEqual(self._load()["beta disease"], "DOID:2")

    def test_exact_synonym_resolves_and_related_does_not(self) -> None:
        names = self._load()
        self.assertEqual(names["gamma disease"], "DOID:2")
        self.assertNotIn("delta disease", names)


class DedupeKeyTest(unittest.TestCase):
    def test_association_type_separates_rows(self) -> None:
        self.assertNotEqual(fdm.dedupe_key(_row()), fdm.dedupe_key(_row(associationType="is_not_model_of")))

    def test_gene_separates_rows(self) -> None:
        self.assertNotEqual(fdm.dedupe_key(_row()), fdm.dedupe_key(_row(subjectGeneId="G2")))

    def test_unmapped_diseases_keep_their_names(self) -> None:
        a = _row(diseaseId="", diseaseName="x")
        b = _row(diseaseId="", diseaseName="y")
        self.assertNotEqual(fdm.dedupe_key(a), fdm.dedupe_key(b))

    def test_identical_rows_collapse(self) -> None:
        self.assertEqual(fdm.dedupe_key(_row()), fdm.dedupe_key(_row()))


class UnmappedRateTest(unittest.TestCase):
    def test_threshold_is_one_percent(self) -> None:
        self.assertFalse(fdm.unmapped_rate_exceeded(10, 1000))
        self.assertTrue(fdm.unmapped_rate_exceeded(11, 1000))

    def test_no_named_rows_passes(self) -> None:
        self.assertFalse(fdm.unmapped_rate_exceeded(0, 0))


if __name__ == "__main__":
    unittest.main()
