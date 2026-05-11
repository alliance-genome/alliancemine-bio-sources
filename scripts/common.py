"""
Shared helpers for the Alliance API fetchers.

Each fetcher in this directory pulls data from the public Alliance API at
https://www.alliancegenome.org/api and emits a TSV on disk that the matching
Java BioFileConverter (under alliance-*/src/main/java/...) reads during the
InterMine build. This module keeps the HTTP plumbing, pagination, retry,
on-disk caching, and TSV writing in one place so the fetchers themselves stay
focused on the per-endpoint mapping from JSON to columns.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator

API_BASE = os.environ.get("ALLIANCE_API_BASE", "https://www.alliancegenome.org/api")
CACHE_DIR = Path(os.environ.get("ALLIANCE_FETCH_CACHE", Path(__file__).parent / ".cache"))

log = logging.getLogger("alliance.fetch")


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-5s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# HTTP GET with retry + cache


@dataclass
class FetchStats:
    requests: int = 0
    cache_hits: int = 0
    retries: int = 0
    errors: int = 0


def _cache_path(cache_name: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"{cache_name}.sqlite"


@contextmanager
def open_cache(cache_name: str) -> Iterator[sqlite3.Connection]:
    """Open the SQLite cache for a fetcher. Keyed by URL hash.
    WAL + busy_timeout so multiple ThreadPool workers can write concurrently
    without hitting "database is locked" / "readonly database" errors."""
    conn = sqlite3.connect(_cache_path(cache_name), timeout=30.0)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=30000")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute(
        "CREATE TABLE IF NOT EXISTS cache ("
        " key TEXT PRIMARY KEY,"
        " url TEXT NOT NULL,"
        " body TEXT NOT NULL,"
        " ts INTEGER NOT NULL"
        ")"
    )
    try:
        yield conn
    finally:
        conn.commit()
        conn.close()


def _url_key(url: str) -> str:
    return hashlib.sha256(url.encode("utf-8")).hexdigest()


def http_get_json(
    path: str,
    *,
    cache: sqlite3.Connection | None = None,
    stats: FetchStats | None = None,
    max_retries: int = 5,
    timeout: float = 30.0,
) -> dict[str, Any]:
    """GET path (prepended with API_BASE) as JSON. Retries with exponential
    backoff on 429/5xx and caches successful responses."""
    url = path if path.startswith("http") else f"{API_BASE}{path}"
    key = _url_key(url)

    if cache is not None:
        row = cache.execute("SELECT body FROM cache WHERE key = ?", (key,)).fetchone()
        if row is not None:
            if stats is not None:
                stats.cache_hits += 1
            return json.loads(row[0])

    delay = 1.0
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "AllianceMine-fetcher/1.0 (+https://github.com/alliance-genome/alliancemine-bio-sources)",
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read().decode("utf-8", errors="replace")
            if stats is not None:
                stats.requests += 1
            if cache is not None:
                cache.execute(
                    "INSERT OR REPLACE INTO cache (key, url, body, ts) VALUES (?, ?, ?, ?)",
                    (key, url, body, int(time.time())),
                )
            return json.loads(body)
        except urllib.error.HTTPError as e:
            # 405 + HTML body = Cloudflare/WAF bot-protection page when concurrency
            # exceeds origin tolerance. Treat as retryable with longer backoff.
            if e.code in (405, 429, 500, 502, 503, 504) and attempt + 1 < max_retries:
                if stats is not None:
                    stats.retries += 1
                log.warning("HTTP %s on %s, retry %d after %.1fs", e.code, path, attempt + 1, delay)
                time.sleep(delay)
                delay = min(delay * 2, 60.0)
                continue
            if stats is not None:
                stats.errors += 1
            raise
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt + 1 < max_retries:
                if stats is not None:
                    stats.retries += 1
                log.warning("URLError on %s (%s), retry %d after %.1fs", path, e, attempt + 1, delay)
                time.sleep(delay)
                delay = min(delay * 2, 30.0)
                continue
            if stats is not None:
                stats.errors += 1
            raise
    raise RuntimeError(f"exhausted retries for {url}")


# ---------------------------------------------------------------------------
# Pagination helper — Alliance API convention: ?limit=N&page=P with
# {results:[...], total, returnedRecords}


def paginate(
    path: str,
    *,
    limit: int = 1000,
    extra_params: dict[str, str] | None = None,
    cache: sqlite3.Connection | None = None,
    stats: FetchStats | None = None,
) -> Iterator[dict[str, Any]]:
    page = 1
    seen = 0
    while True:
        params = {"limit": str(limit), "page": str(page)}
        if extra_params:
            params.update(extra_params)
        query = urllib.parse.urlencode(params)
        sep = "&" if "?" in path else "?"
        data = http_get_json(f"{path}{sep}{query}", cache=cache, stats=stats)
        results = data.get("results") or []
        total = data.get("total", 0)
        for r in results:
            yield r
        seen += len(results)
        if not results or seen >= total:
            break
        page += 1


# ---------------------------------------------------------------------------
# TSV writer


class TsvWriter:
    """Small wrapper that atomically writes a TSV with a header comment block
    and named columns. Writes to a .tmp file then renames on close, so readers
    never see a partial file."""

    def __init__(self, path: Path, columns: list[str], *, release: str | None = None, source: str | None = None) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.columns = columns
        self.tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        self.rows_written = 0
        self._fh = open(self.tmp, "w", encoding="utf-8")
        self._write_header(release=release, source=source)

    def _write_header(self, *, release: str | None, source: str | None) -> None:
        self._fh.write("#" * 70 + "\n")
        self._fh.write(f"# Generated by {Path(sys.argv[0]).name} on {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
        if source:
            self._fh.write(f"# Source: {source}\n")
        if release:
            self._fh.write(f"# Alliance release: {release}\n")
        self._fh.write("#" * 70 + "\n")
        self._fh.write("\t".join(self.columns) + "\n")

    def write_row(self, **fields: Any) -> None:
        row = []
        for col in self.columns:
            val = fields.get(col, "")
            if val is None:
                val = ""
            # Flatten pipes / tabs / newlines that would break TSV parsing
            s = str(val).replace("\t", " ").replace("\n", " ").replace("\r", " ")
            row.append(s)
        self._fh.write("\t".join(row) + "\n")
        self.rows_written += 1

    def close(self) -> None:
        self._fh.close()
        os.replace(self.tmp, self.path)
        log.info("Wrote %s (%d rows)", self.path, self.rows_written)

    def __enter__(self) -> "TsvWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is None:
            self.close()
        else:
            self._fh.close()
            if self.tmp.exists():
                self.tmp.unlink()


# ---------------------------------------------------------------------------
# Helpers for digging into the Alliance JSON shape


def get_nested(obj: Any, *path: str, default: Any = None) -> Any:
    """obj['a']['b']['c'] with graceful None fallback."""
    cur = obj
    for p in path:
        if cur is None or not isinstance(cur, dict):
            return default
        cur = cur.get(p)
    return cur if cur is not None else default


def join_pipe(values: Iterable[Any]) -> str:
    """Pipe-join non-empty stringifications."""
    return "|".join(str(v) for v in values if v not in (None, ""))


# ---------------------------------------------------------------------------
# Track D4: collapse legacy/ambiguous taxon IDs to AllianceMine canonical IDs.
# Without this, sources tagging yeast to 4932 produce ghost Organism rows
# alongside the canonical S288C (559292) row, splitting yeast genes/alleles
# across two organisms.

LEGACY_TAXON_REMAP: dict[str, str] = {
    "4932": "559292",                      # S. cerevisiae generic -> S288C reference
    "NCBITaxon:4932": "NCBITaxon:559292",
    "taxon:4932": "taxon:559292",
}


def normalize_taxon(taxon_id: str | None) -> str:
    """Map legacy/ambiguous taxon IDs to canonical AllianceMine IDs.
    Idempotent — returns input unchanged if no remap rule applies, or "" if
    the input is None."""
    if taxon_id is None:
        return ""
    return LEGACY_TAXON_REMAP.get(taxon_id, taxon_id)


# ---------------------------------------------------------------------------
# Release discovery


def get_current_release() -> str:
    """Fetch the current Alliance release version from the legacy FMS endpoint.
    Used for provenance in the TSV header. Falls back to 'unknown' if FMS
    is fully decommissioned by the time this runs."""
    try:
        data = http_get_json("https://fms.alliancegenome.org/api/releaseversion/all")
        if isinstance(data, list) and data:
            latest = sorted(data, key=lambda x: x.get("releaseDate", ""), reverse=True)[0]
            return latest.get("releaseVersion", "unknown")
    except Exception as e:
        log.debug("Could not fetch release version: %s", e)
    return "unknown"


# ---------------------------------------------------------------------------
# Gene-ID enumeration


# ---------------------------------------------------------------------------
# InterMine PathQuery REST client (for cross-mine federation in Phase 6d)


def intermine_paginate(
    base_url: str,
    query_xml: str,
    *,
    page_size: int = 1000,
    cache: sqlite3.Connection | None = None,
    stats: FetchStats | None = None,
) -> Iterator[list[Any]]:
    """Yield each row from a paginated InterMine PathQuery REST call.

    Each row is the JSON-array form: [val1, val2, ...] in the order the query's
    view declares. The caller maps positional values to TSV columns.

    base_url should NOT include /service - we append /service/query/results.
    """
    import urllib.parse
    import urllib.error
    base = base_url.rstrip("/")
    if not base.endswith("/service"):
        base = base + "/service"

    start = 0
    while True:
        body = urllib.parse.urlencode({
            "query": query_xml,
            "format": "json",
            "size": str(page_size),
            "start": str(start),
        }).encode()
        url = f"{base}/query/results"

        cache_key = _url_key(f"{url}|{start}|{page_size}|{query_xml}")
        if cache is not None:
            row = cache.execute("SELECT body FROM cache WHERE key = ?", (cache_key,)).fetchone()
            if row is not None:
                if stats is not None:
                    stats.cache_hits += 1
                data = json.loads(row[0])
                results = data.get("results") or []
                for r in results:
                    yield r
                if len(results) < page_size:
                    return
                start += page_size
                continue

        delay = 1.0
        for attempt in range(5):
            try:
                req = urllib.request.Request(url, data=body)
                with urllib.request.urlopen(req, timeout=60) as resp:
                    text = resp.read().decode("utf-8", errors="replace")
                if stats is not None:
                    stats.requests += 1
                if cache is not None:
                    cache.execute(
                        "INSERT OR REPLACE INTO cache (key, url, body, ts) VALUES (?, ?, ?, ?)",
                        (cache_key, url, text, int(time.time())),
                    )
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503, 504) and attempt + 1 < 5:
                    if stats is not None:
                        stats.retries += 1
                    log.warning("HTTP %s on PathQuery, retry %d after %.1fs", e.code, attempt + 1, delay)
                    time.sleep(delay)
                    delay = min(delay * 2, 30.0)
                    continue
                if stats is not None:
                    stats.errors += 1
                raise
        else:
            raise RuntimeError(f"exhausted retries for {url}")

        data = json.loads(text)
        if not data.get("wasSuccessful", True):
            raise RuntimeError(f"PathQuery error: {data.get('error', 'unknown')}")
        results = data.get("results") or []
        for r in results:
            yield r
        if len(results) < page_size:
            return
        start += page_size


def enumerate_yeast_genes() -> list[str]:
    """Return SGD gene CURIEs for Saccharomyces cerevisiae S288C.

    First tries the FMS BGI SGD JSON file (while FMS is up). If that fails we
    fall back to parsing SGD's public gene list. For a multi-MOD build extend
    this with per-MOD seeders."""
    return enumerate_mod_genes(["SGD"])


def _resolve_bgi_urls(mods: list[str]) -> dict[str, str]:
    """Use the FMS snapshot API to find the current BGI URL for each MOD.

    Different MODs may have fallen back to older releases in any given snapshot
    (e.g. WB/MGI pinning to 8.3.0 even when the active release is 9.0.0), so we
    ask the snapshot API to tell us which version is current rather than
    hardcoding a release-specific URL."""
    release = get_current_release()
    url = f"https://fms.alliancegenome.org/api/snapshot/release/{release}"
    data = json.loads(urllib.request.urlopen(url, timeout=30).read())
    bgi_by_mod: dict[str, str] = {}
    for f in data.get("snapShot", {}).get("dataFiles", []):
        if f.get("dataType", {}).get("name") != "BGI":
            continue
        if not f.get("s3Path", "").endswith(".json.gz"):
            continue
        sub = f.get("dataSubType", {}).get("name", "")
        if sub in mods:
            bgi_by_mod[sub] = f["s3Url"]
    return bgi_by_mod


# MOD abbreviation -> curie prefix that its gene CURIEs carry.
_MOD_PREFIX = {
    "SGD": "SGD:",
    "MGI": "MGI:",
    "RGD": "RGD:",
    "ZFIN": "ZFIN:",
    "FB": "FB:",
    "WB": "WB:",
    "XBXL": "Xenbase:",
    "XBXT": "Xenbase:",
    "HUMAN": "HGNC:",
}


def enumerate_mod_genes(mods: list[str]) -> list[str]:
    """Return gene CURIEs across one or more MODs by parsing each MOD's BGI
    JSON from FMS. Mods are abbreviations like "SGD", "MGI", "RGD", "ZFIN",
    "FB", "WB", "XBXL", "XBXT", "HUMAN"."""
    import gzip
    try:
        urls = _resolve_bgi_urls(mods)
    except Exception as e:
        log.error("Could not resolve BGI URLs via snapshot API: %s", e)
        raise
    missing = [m for m in mods if m not in urls]
    if missing:
        log.warning("No BGI file found in snapshot for MOD(s): %s", missing)

    combined: list[str] = []
    for mod in mods:
        if mod not in urls:
            continue
        url = urls[mod]
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                raw = resp.read()
            text = gzip.decompress(raw).decode("utf-8")
            d = json.loads(text)
        except Exception as e:
            log.error("Failed to fetch BGI for %s: %s", mod, e)
            continue
        prefix = _MOD_PREFIX.get(mod)
        added = 0
        for gene in d.get("data", []):
            xrefs = gene.get("basicGeneticEntity", {}).get("crossReferences", [])
            for x in xrefs:
                xid = x.get("id", "")
                if prefix and xid.startswith(prefix):
                    combined.append(xid)
                    added += 1
                    break
        log.info("Enumerated %d %s genes from FMS BGI", added, mod)
    log.info("Total across MODs: %d genes", len(combined))
    return combined
