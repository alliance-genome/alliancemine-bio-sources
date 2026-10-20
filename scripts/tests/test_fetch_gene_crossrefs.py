"""Unit tests for scripts/fetch_gene_crossrefs.py."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import fetch_gene_crossrefs  # noqa: E402


class IterFilesTest(unittest.TestCase):
    def test_xenbase_gene_ids_are_kept(self) -> None:
        records = [
            ("Xenbase", {"GeneID": "Xenbase:XB-GENE-1", "GlobalCrossReferenceID": "NCBI_Gene:1"}),
            ("Other", {"GeneID": "UNKNOWN:1", "GlobalCrossReferenceID": "NCBI_Gene:2"}),
        ]
        with mock.patch.object(fetch_gene_crossrefs, "iter_fms_datatype", return_value=records):
            rows = list(fetch_gene_crossrefs._iter_files())
        self.assertEqual([r["geneId"] for r in rows], ["Xenbase:XB-GENE-1"])


if __name__ == "__main__":
    unittest.main()
