"""Failure-path tests for scripts/common.py. Run: python3 -m pytest scripts/tests"""

from __future__ import annotations

import io
import sqlite3
import sys
import time
import urllib.error
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common  # noqa: E402


@pytest.fixture(autouse=True)
def pinned_release(monkeypatch, tmp_path):
    monkeypatch.setenv("ALLIANCE_RELEASE", "9.2.0")
    common.get_current_release.cache_clear()
    monkeypatch.setattr(common, "CACHE_DIR", tmp_path / "cache")
    monkeypatch.setattr(common, "FMS_CACHE_DIR", tmp_path / "fms")
    yield
    common.get_current_release.cache_clear()


def _report(failures: int, total: int = 100) -> common.RunReport:
    return common.RunReport(
        label="genes", total=total, done=total,
        failures=[(f"G{i}", "x") for i in range(failures)],
    )


def test_threshold_exit_inside_writer_leaves_no_file(tmp_path):
    target = tmp_path / "out.tsv"
    with pytest.raises(SystemExit):
        with common.TsvWriter(target, ["a"]) as writer:
            writer.write_row(a="1")
            common.finish_run(_report(5), common.FetchStats(), tmp_path, 0.005,
                              rows_written=writer.rows_written)
    assert not target.exists()
    assert not target.with_suffix(".tsv.tmp").exists()


def test_zero_rows_exits_and_leaves_no_file(tmp_path):
    target = tmp_path / "out.tsv"
    with pytest.raises(SystemExit):
        with common.TsvWriter(target, ["a"]) as writer:
            common.finish_run(_report(0), common.FetchStats(), tmp_path, 0.005,
                              rows_written=writer.rows_written)
    assert not target.exists()


def test_success_inside_writer_renames_file(tmp_path):
    target = tmp_path / "out.tsv"
    with common.TsvWriter(target, ["a"]) as writer:
        writer.write_row(a="1")
        common.finish_run(_report(0), common.FetchStats(), tmp_path, 0.005,
                          rows_written=writer.rows_written)
    assert target.exists()
    assert not target.with_suffix(".tsv.tmp").exists()


@pytest.mark.parametrize("code", [400, 404])
def test_http_400_and_404_count_as_failures(code):
    def work(eid, cache, stats):
        raise urllib.error.HTTPError(eid, code, "nope", {}, None)

    report = common.run_parallel(
        ["G1", "G2"], work, lambda *_: None,
        cache_name="t", stats=common.FetchStats(), workers=1,
    )
    assert report.failures == [("G1", f"HTTP {code}"), ("G2", f"HTTP {code}")]
    assert report.failure_rate == 1.0


def test_fms_download_failure_leaves_no_target(monkeypatch):
    class Boom(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self, n=-1):
            raise ConnectionError("cut")

    monkeypatch.setattr(common.urllib.request, "urlopen", lambda *a, **k: Boom())
    with pytest.raises(ConnectionError):
        common.fms_download_gz("https://x/f.json.gz")
    assert list(common.FMS_CACHE_DIR.iterdir()) == []


def test_missing_fms_mod_raises(monkeypatch):
    monkeypatch.setattr(common, "fms_list_files", lambda *a, **k: [])
    with pytest.raises(RuntimeError, match="No FMS files"):
        list(common.iter_fms_datatype("AGM", sub_filter="WB"))
