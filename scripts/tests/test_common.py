"""Unit tests for scripts/common.py. Run: python3 -m unittest discover -s scripts/tests"""

from __future__ import annotations

import os
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common  # noqa: E402


class ReleaseTest(unittest.TestCase):
    def setUp(self) -> None:
        common.get_current_release.cache_clear()
        common.get_api_release.cache_clear()

    tearDown = setUp

    def test_env_pin_wins_without_network(self) -> None:
        with mock.patch.dict(os.environ, {"ALLIANCE_RELEASE": "9.2.0"}), \
             mock.patch.object(common, "http_get_json", side_effect=AssertionError("no HTTP")):
            self.assertEqual(common.get_current_release(), "9.2.0")

    def test_unpinned_follows_api_not_newest_fms(self) -> None:
        def fake(path, **_kw):
            if path == "/releaseInfo":
                return {"releaseVersion": "9.1.0"}
            raise AssertionError(path)
        env = {k: v for k, v in os.environ.items() if k != "ALLIANCE_RELEASE"}
        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch.object(common, "http_get_json", side_effect=fake):
            self.assertEqual(common.get_current_release(), "9.1.0")

    def test_falls_back_to_fms_current(self) -> None:
        def fake(path, **_kw):
            if path == "/releaseInfo":
                raise OSError("down")
            if path.endswith("/releaseversion/current"):
                return {"releaseVersion": "9.0.0"}
            raise AssertionError(path)
        env = {k: v for k, v in os.environ.items() if k != "ALLIANCE_RELEASE"}
        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch.object(common, "http_get_json", side_effect=fake):
            self.assertEqual(common.get_current_release(), "9.0.0")


class ModsTest(unittest.TestCase):
    def test_cli_wins(self) -> None:
        with mock.patch.dict(os.environ, {"ALLIANCE_MODS": "SGD"}):
            self.assertEqual(common.selected_mods(["MGI"], "mgi,FB"), ["MGI", "FB"])

    def test_env_narrows_defaults(self) -> None:
        with mock.patch.dict(os.environ, {"ALLIANCE_MODS": "SGD,MGI"}):
            self.assertEqual(common.selected_mods(["MGI", "RGD"]), ["MGI"])
            self.assertEqual(common.selected_mods(), ["SGD", "MGI"])

    def test_default_is_all_nine(self) -> None:
        env = {k: v for k, v in os.environ.items() if k != "ALLIANCE_MODS"}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertEqual(len(common.selected_mods()), 9)

    def test_unknown_mod_rejected(self) -> None:
        with self.assertRaises(ValueError):
            common.parse_mods("SGD,XB")

    def test_ids_and_limit(self) -> None:
        self.assertEqual(common.resolve_gene_ids("A, B,C", None, 2), ["A", "B"])


class FmsListFilesTest(unittest.TestCase):
    def setUp(self) -> None:
        common._fms_manifest.cache_clear()

    tearDown = setUp

    def test_keeps_latest_upload_per_subtype(self) -> None:
        def row(sub, path, date, dt="BGI"):
            return {"dataType": {"name": dt}, "dataSubType": {"name": sub},
                    "s3Url": f"https://x/{path}", "uploadDate": date}
        manifest = [
            row("SGD", "9.0.0/sgd", "2026-03-18"),
            row("SGD", "9.2.0/sgd", "2026-09-24"),
            row("ZFIN", "9.1.0/zfin3", "2026-08-04"),
            row("ZFIN", "9.1.0/zfin0", "2026-07-07"),
            row("SGD", "9.2.0/daf", "2026-09-24", dt="DAF"),
        ]
        with mock.patch.object(common, "http_get_json", return_value=manifest):
            files = common.fms_list_files("BGI", release="9.2.0")
        self.assertEqual(sorted(f["url"] for f in files), ["https://x/9.1.0/zfin3", "https://x/9.2.0/sgd"])


class RunParallelTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.patch = mock.patch.object(common, "CACHE_DIR", Path(self.tmp.name))
        self.patch.start()

    def tearDown(self) -> None:
        self.patch.stop()
        self.tmp.cleanup()

    def test_order_failures_and_thread_local_cache(self) -> None:
        ids = [f"G{i}" for i in range(50)]
        seen_conns: dict[int, int] = {}
        lock = threading.Lock()

        def work(gid, cache, stats):
            with lock:
                seen_conns.setdefault(threading.get_ident(), id(cache))
                self.assertEqual(seen_conns[threading.get_ident()], id(cache))
            cache.execute("INSERT OR REPLACE INTO cache VALUES (?, ?, ?, ?)", (gid, gid, gid, 0))
            time.sleep(0.001 * (hash(gid) % 5))
            if gid in ("G7", "G33"):
                raise RuntimeError("boom")
            return [gid.lower()]

        consumed: list[str] = []
        stats = common.FetchStats()
        report = common.run_parallel(
            ids, work, lambda gid, rows: consumed.extend(rows),
            cache_name="t", stats=stats, workers=4, label="genes",
        )
        expected = [g.lower() for g in ids if g not in ("G7", "G33")]
        self.assertEqual(consumed, expected)
        self.assertEqual([f[0] for f in report.failures], ["G7", "G33"])
        self.assertEqual(report.done, 50)
        with common.open_cache("t") as conn:
            self.assertEqual(conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0], 50)

    def test_finish_run_threshold(self) -> None:
        stats = common.FetchStats()
        ok = common.RunReport(label="genes", total=1000, done=1000, failures=[("G1", "x")])
        common.finish_run(ok, stats, Path(self.tmp.name), 0.005)
        self.assertTrue((Path(self.tmp.name) / f"{Path(sys.argv[0]).stem}.failures.tsv").exists())
        bad = common.RunReport(label="genes", total=100, done=100, failures=[("G1", "x")])
        with self.assertRaises(SystemExit):
            common.finish_run(bad, stats, Path(self.tmp.name), 0.005)


class StatsTest(unittest.TestCase):
    def test_incr_is_thread_safe(self) -> None:
        stats = common.FetchStats()
        threads = [threading.Thread(target=lambda: [stats.incr("requests") for _ in range(10000)])
                   for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertEqual(stats.requests, 80000)


if __name__ == "__main__":
    unittest.main()
