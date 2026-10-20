"""Unit tests for the gene-crossrefs filter in scripts/fetch_chembl.py."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest import mock

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common  # noqa: E402
import fetch_chembl  # noqa: E402

CROSSREFS = (
    "# header\n"
    "geneId\tgeneTaxon\txrefId\txrefUrl\txrefType\n"
    "SGD:S000000001\t559292\tUniProtKB:P11111\thttp://x\tgene\n"
    "SGD:S000000001\t559292\tUniProtKB:P11111\thttp://x\tgene\n"
    "SGD:S000000002\t559292\tNCBI_Gene:5\thttp://x\tgene\n"
    "short\trow\n"
)
MAPPING = (
    "# comment\n"
    "P11111\tCHEMBL1\tTarget one\tSINGLE PROTEIN\n"
    "Q99999\tCHEMBL2\tTarget two\tSINGLE PROTEIN\n"
)


class _Resp:
    def read(self) -> bytes:
        return MAPPING.encode()


def test_rows_from_crossrefs_keeps_only_uniprot_rows(tmp_path):
    path = tmp_path / "gene-crossrefs.tsv"
    path.write_text(CROSSREFS)
    rows = list(common.rows_from_crossrefs(path))
    assert [r["uniProtAccession"] for r in rows] == ["P11111", "P11111"]
    assert rows[0]["geneTaxon"] == "559292"


def _run(tmp_path, *extra):
    argv = ["fetch_chembl.py", "--out-dir", str(tmp_path), *extra]
    with mock.patch.object(sys, "argv", argv), \
            mock.patch.object(fetch_chembl, "get_current_release", return_value="9.2.0"), \
            mock.patch.object(fetch_chembl.urllib.request, "urlopen", return_value=_Resp()):
        fetch_chembl.main()


def test_main_drops_accessions_missing_from_crossrefs(tmp_path):
    (tmp_path / "gene-crossrefs.tsv").write_text(CROSSREFS)
    _run(tmp_path)
    rows = [l.split("\t") for l in (tmp_path / "chembl-targets.tsv").read_text().splitlines() if not l.startswith("#")]
    assert rows[0][0] == "uniProtAccession"
    assert [r[0] for r in rows[1:]] == ["P11111"]
    assert rows[1][-1] == "559292"


def test_main_exits_when_crossrefs_missing(tmp_path):
    with pytest.raises(SystemExit) as exc:
        _run(tmp_path)
    assert exc.value.code == 1
    assert not (tmp_path / "chembl-targets.tsv").exists()
