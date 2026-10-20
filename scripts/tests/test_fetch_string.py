"""Unit tests for scripts/fetch_string.py pure functions."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import fetch_string  # noqa: E402


def _gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def test_pick_accession_prefers_known_over_first():
    assert fetch_string.pick_accession(["D6VPL5", "P34111"], {"P34111"}) == "P34111"


def test_pick_accession_empty_when_none_known():
    assert fetch_string.pick_accession(["D6VPL5"], {"P34111"}) == ""


def test_pick_accession_first_without_known_set():
    assert fetch_string.pick_accession(["D6VPL5", "P34111"], None) == "D6VPL5"
    assert fetch_string.pick_accession([], None) == ""


def test_alias_map_uses_known_accession(tmp_path):
    aliases = _gz(
        tmp_path / "a.gz",
        "#string_protein_id\talias\tsource\n"
        "4932.YAL001C\tD6VPL5\tUniProt_AC\n"
        "4932.YAL001C\tP34111\tUniProt_AC\n"
        "4932.YAL001C\tTFC3\tUniProt_ID\n"
        "4932.YAL002W\tQ00000\tUniProt_AC\n",
    )
    with mock.patch.object(fetch_string, "_download", return_value=aliases):
        assert fetch_string._alias_map("4932", {"P34111"}) == {"4932.YAL001C": "P34111"}
        assert fetch_string._alias_map("4932") == {"4932.YAL001C": "D6VPL5", "4932.YAL002W": "Q00000"}


def test_iter_links_emits_each_pair_once(tmp_path):
    links = _gz(
        tmp_path / "l.gz",
        "protein1 protein2 combined_score\n"
        "4932.A 4932.B 999\n"
        "4932.B 4932.A 999\n"
        "4932.A 4932.C 100\n"
        "4932.C 4932.D 800\n",
    )
    with mock.patch.object(fetch_string, "_download", return_value=links):
        rows = list(fetch_string._iter_links("4932", 700))
    assert rows == [("4932.A", "4932.B", 999), ("4932.C", "4932.D", 800)]
