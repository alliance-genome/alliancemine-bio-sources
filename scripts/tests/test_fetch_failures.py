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


def test_cache_key_includes_release_and_api_base(monkeypatch):
    url = "https://x/gene/1"
    k920 = common._url_key(url)
    monkeypatch.setenv("ALLIANCE_RELEASE", "9.1.0")
    common.get_current_release.cache_clear()
    assert common._url_key(url) != k920
    monkeypatch.setenv("ALLIANCE_RELEASE", "9.2.0")
    common.get_current_release.cache_clear()
    assert common._url_key(url) == k920
    monkeypatch.setattr(common, "API_BASE", "https://stage/api")
    assert common._url_key(url) != k920


def test_cached_response_from_other_release_is_not_served(monkeypatch):
    calls = []

    class Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        calls.append(req.full_url)
        return Resp(b'{"release": "%d"}' % len(calls))

    monkeypatch.setattr(common.urllib.request, "urlopen", fake_urlopen)
    conn = common.connect_cache("t")
    assert common.http_get_json("/gene/1", cache=conn) == {"release": "1"}
    assert common.http_get_json("/gene/1", cache=conn) == {"release": "1"}
    monkeypatch.setenv("ALLIANCE_RELEASE", "9.1.0")
    common.get_current_release.cache_clear()
    assert common.http_get_json("/gene/1", cache=conn) == {"release": "2"}
    assert len(calls) == 2


def test_pathquery_cache_expires_by_ts(monkeypatch):
    calls = []

    class Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        calls.append(1)
        return Resp(b'{"results": [["a"]], "wasSuccessful": true}')

    monkeypatch.setattr(common.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(common, "PATHQUERY_MAX_AGE", 100)
    conn = common.connect_cache("pq")
    run = lambda: list(common.intermine_paginate("http://m", "<query/>", cache=conn))
    assert run() == [["a"]]
    assert run() == [["a"]]
    assert len(calls) == 1
    conn.execute("UPDATE cache SET ts = ?", (int(time.time()) - 1000,))
    assert run() == [["a"]]
    assert len(calls) == 2


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


@pytest.mark.parametrize("module, datatypes", [
    ("fetch_agms", ["AGM"]),
    ("fetch_htp", ["HTPDATASET", "HTPDATASAMPLE"]),
    ("fetch_disease_alleles", ["DAF"]),
    ("fetch_experimental_disease", ["DAF"]),
])
def test_default_mods_all_have_fms_files(monkeypatch, module, datatypes):
    # 9.1.0 manifest: every MOD but XBXL and XBXT has a file for these datatypes.
    manifest_mods = {"MGI", "ZFIN", "FB", "WB", "RGD", "SGD", "HUMAN"}
    rows = [
        {"dataType": {"name": dt}, "dataSubType": {"name": m},
         "s3Url": f"x/{dt}_{m}.json.gz", "uploadDate": "2026-01-01"}
        for dt in datatypes for m in manifest_mods
    ]
    monkeypatch.setattr(common, "_fms_manifest", lambda release: rows)
    monkeypatch.setattr(common, "fms_download_gz", lambda url: Path(url))
    monkeypatch.setattr(common, "fms_iter_json_records", lambda path: iter([{}]))
    monkeypatch.delenv("ALLIANCE_MODS", raising=False)
    mod = __import__(module)
    mods = common.selected_mods(mod.FMS_DEFAULT_MODS, None)
    assert mods
    for dt in datatypes:
        for m in mods:
            assert list(common.iter_fms_datatype(dt, sub_filter=m)) == [(m, {})]
