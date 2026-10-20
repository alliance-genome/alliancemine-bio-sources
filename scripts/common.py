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

import functools
import hashlib
import json
import logging
import os
import sqlite3
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator

API_BASE = os.environ.get("ALLIANCE_API_BASE", "https://www.alliancegenome.org/api")
FMS_BASE = os.environ.get("ALLIANCE_FMS_BASE", "https://fms.alliancegenome.org/api")
# The site WAF answers HTTP 405 to bulk clients with a non-browser User-Agent
# (stage.alliancegenome.org, 2026-10-14); ALLIANCE_USER_AGENT overrides.
USER_AGENT = os.environ.get(
    "ALLIANCE_USER_AGENT",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 14_7) AppleWebKit/605.1.15 "
    "(KHTML, like Gecko) Version/17.6 Safari/605.1.15",
)
CACHE_DIR = Path(os.environ.get("ALLIANCE_FETCH_CACHE", Path(__file__).parent / ".cache"))
FMS_CACHE_DIR = Path(os.environ.get("ALLIANCE_FMS_CACHE", Path(__file__).parent / ".fms_cache"))

# Seconds a cached PathQuery page stays valid; other mines publish on their own schedule.
PATHQUERY_MAX_AGE = int(os.environ.get("ALLIANCE_PATHQUERY_MAX_AGE", str(7 * 86400)))

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
    throttled: int = 0
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False, compare=False)

    def incr(self, name: str, n: int = 1) -> None:
        with self._lock:
            setattr(self, name, getattr(self, name) + n)


def _bump(stats: FetchStats | None, name: str) -> None:
    if stats is not None:
        stats.incr(name)


def _cache_path(cache_name: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"{cache_name}.sqlite"


def connect_cache(cache_name: str, *, initialise: bool = True) -> sqlite3.Connection:
    """Open one SQLite cache connection in autocommit mode, so every cached
    response is durable immediately (reruns resume from it) and no thread
    holds the write lock between requests. WAL + busy_timeout let several
    connections (one per worker thread) share the file. Switching to WAL
    needs an exclusive lock, so worker connections skip initialisation and
    rely on it having been done once up front."""
    conn = sqlite3.connect(
        _cache_path(cache_name), timeout=30.0, isolation_level=None, check_same_thread=False
    )
    conn.execute("PRAGMA busy_timeout=30000")
    conn.execute("PRAGMA synchronous=NORMAL")
    if initialise:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute(
            "CREATE TABLE IF NOT EXISTS cache ("
            " key TEXT PRIMARY KEY,"
            " url TEXT NOT NULL,"
            " body TEXT NOT NULL,"
            " ts INTEGER NOT NULL"
            ")"
        )
    return conn


@contextmanager
def open_cache(cache_name: str) -> Iterator[sqlite3.Connection]:
    """Open the SQLite cache for a fetcher. Keyed by URL hash."""
    conn = connect_cache(cache_name)
    try:
        yield conn
    finally:
        conn.close()


def _url_key(url: str) -> str:
    """Cache key: sha256 of release, API base and URL, so a payload cached
    for one release or API host is never served for another."""
    scope = f"{get_current_release()}|{API_BASE}|{url}"
    return hashlib.sha256(scope.encode("utf-8")).hexdigest()


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
    key = _url_key(url) if cache is not None else ""

    if cache is not None:
        row = cache.execute("SELECT body FROM cache WHERE key = ?", (key,)).fetchone()
        if row is not None:
            _bump(stats, "cache_hits")
            return json.loads(row[0])

    delay = 1.0
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": USER_AGENT,
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read().decode("utf-8", errors="replace")
            _bump(stats, "requests")
            payload = json.loads(body)
            if cache is not None:
                cache.execute(
                    "INSERT OR REPLACE INTO cache (key, url, body, ts) VALUES (?, ?, ?, ?)",
                    (key, url, body, int(time.time())),
                )
            return payload
        except urllib.error.HTTPError as e:
            # 405 + HTML body = Cloudflare/WAF bot-protection page when concurrency
            # exceeds origin tolerance. Treat as retryable with longer backoff.
            if e.code in (405, 429):
                _bump(stats, "throttled")
            if e.code in (405, 429, 500, 502, 503, 504) and attempt + 1 < max_retries:
                _bump(stats, "retries")
                log.warning("HTTP %s on %s, retry %d after %.1fs", e.code, path, attempt + 1, delay)
                time.sleep(delay)
                delay = min(delay * 2, 60.0)
                continue
            _bump(stats, "errors")
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if attempt + 1 < max_retries:
                _bump(stats, "retries")
                log.warning("URLError on %s (%s), retry %d after %.1fs", path, e, attempt + 1, delay)
                time.sleep(delay)
                delay = min(delay * 2, 30.0)
                continue
            _bump(stats, "errors")
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


def publication_ids(result: dict) -> list[str]:
    """PubMed + MOD publication curies of an annotation. The API replaced the
    flat pubmedPubModIDs list with pubmedPublications / pubModPublications
    cross-reference lists; read whichever is present."""
    flat = result.get("pubmedPubModIDs")
    if flat:
        return list(flat)
    ids: list[str] = []
    for key in ("pubmedPublications", "pubModPublications"):
        for pub in result.get(key) or []:
            curie = pub.get("referencedCurie", "") if isinstance(pub, dict) else ""
            if curie and curie not in ids:
                ids.append(curie)
    return ids


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


@functools.lru_cache(maxsize=None)
def get_api_release() -> str | None:
    """Release the live Alliance API is serving (GET /releaseInfo)."""
    try:
        data = http_get_json("/releaseInfo", max_retries=3)
    except Exception as e:
        log.warning("Could not read API releaseInfo: %s", e)
        return None
    version = data.get("releaseVersion") if isinstance(data, dict) else None
    return version or None


@functools.lru_cache(maxsize=None)
def get_current_release() -> str:
    """The Alliance release every fetcher works against.

    ALLIANCE_RELEASE pins it (the build sets this so FMS file selection, gene
    seeding and TSV provenance all agree). Unpinned runs follow what the API
    serves, then FMS's own "current" release, never the newest release by
    date: FMS lists the next release weeks before it goes live."""
    pinned = os.environ.get("ALLIANCE_RELEASE", "").strip()
    if pinned:
        return pinned
    api_release = get_api_release()
    if api_release:
        log.warning("ALLIANCE_RELEASE not set; using API release %s", api_release)
        return api_release
    try:
        data = http_get_json(f"{FMS_BASE}/releaseversion/current", max_retries=3)
        version = data.get("releaseVersion") if isinstance(data, dict) else None
        if version:
            log.warning("ALLIANCE_RELEASE not set; using FMS current release %s", version)
            return version
    except Exception as e:
        log.warning("Could not read FMS current release: %s", e)
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
            row = cache.execute(
                "SELECT body, ts FROM cache WHERE key = ?", (cache_key,)
            ).fetchone()
            if row is not None and time.time() - row[1] <= PATHQUERY_MAX_AGE:
                _bump(stats, "cache_hits")
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
                _bump(stats, "requests")
                if cache is not None:
                    cache.execute(
                        "INSERT OR REPLACE INTO cache (key, url, body, ts) VALUES (?, ?, ?, ?)",
                        (cache_key, url, text, int(time.time())),
                    )
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503, 504) and attempt + 1 < 5:
                    _bump(stats, "retries")
                    log.warning("HTTP %s on PathQuery, retry %d after %.1fs", e.code, attempt + 1, delay)
                    time.sleep(delay)
                    delay = min(delay * 2, 30.0)
                    continue
                _bump(stats, "errors")
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
    """Return SGD gene CURIEs for Saccharomyces cerevisiae S288C."""
    return enumerate_mod_genes(["SGD"])


# MOD abbreviation (FMS BGI dataSubType) -> curie prefix of its gene CURIEs.
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

ALL_MODS = ["SGD", "MGI", "RGD", "ZFIN", "FB", "WB", "XBXL", "XBXT", "HUMAN"]

MOD_TAXON = {
    "SGD": "NCBITaxon:559292",
    "MGI": "NCBITaxon:10090",
    "RGD": "NCBITaxon:10116",
    "ZFIN": "NCBITaxon:7955",
    "FB": "NCBITaxon:7227",
    "WB": "NCBITaxon:6239",
    "XBXL": "NCBITaxon:8355",
    "XBXT": "NCBITaxon:8364",
    "HUMAN": "NCBITaxon:9606",
}


def parse_mods(value: str | None) -> list[str]:
    mods = [m.strip().upper() for m in (value or "").split(",") if m.strip()]
    unknown = [m for m in mods if m not in _MOD_PREFIX]
    if unknown:
        raise ValueError(f"Unknown MOD(s) {unknown}; expected a subset of {ALL_MODS}")
    return mods


def selected_mods(defaults: list[str] | None = None, cli_value: str | None = None) -> list[str]:
    """MODs a fetcher should cover. An explicit --mods wins; otherwise the
    fetcher's own defaults, narrowed by ALLIANCE_MODS when that is set (so a
    single env var scopes a whole fetch_all run without widening a fetcher
    past the MODs its endpoint supports)."""
    if cli_value:
        return parse_mods(cli_value)
    base = list(defaults) if defaults is not None else list(ALL_MODS)
    env_mods = parse_mods(os.environ.get("ALLIANCE_MODS"))
    if env_mods:
        return [m for m in base if m in env_mods]
    return base


def resolve_gene_ids(
    ids: str | None,
    mods: str | None,
    limit: int | None,
    *,
    default_mods: list[str] | None = None,
) -> list[str]:
    """Seed list for the per-gene fetchers: --ids if given, else every gene
    of the selected MODs from the pinned release's BGI files."""
    if ids:
        gene_ids = [g.strip() for g in ids.split(",") if g.strip()]
    else:
        gene_ids = enumerate_mod_genes(selected_mods(default_mods, mods))
    if limit:
        gene_ids = gene_ids[:limit]
    return gene_ids


# BGI files stay at each MOD's last upload (WB 8.3.0, FB 9.0.0 in the 9.2.0
# manifest). Alliance regenerates GENE-DESCRIPTION-JSON from the release
# database, so it lists genes added since: 760 WB, 1,875 FB and 22,985 RGD
# genes in 9.2.0 that the stage API serves and BGI lacks.
_GENE_LIST_DATATYPES = ("GENE-DESCRIPTION-JSON", "BGI")


def _resolve_gene_list_urls(mods: list[str]) -> dict[str, list[str]]:
    """Gene-list file URLs per MOD for the pinned release, newest upload of
    each datatype in _GENE_LIST_DATATYPES order."""
    urls: dict[str, list[str]] = {}
    for datatype in _GENE_LIST_DATATYPES:
        for f in fms_list_files(datatype):
            if f["subType"] in mods and f["url"].endswith(".json.gz"):
                urls.setdefault(f["subType"], []).append(f["url"])
    return urls


def _bgi_gene_ids(path: Path, prefix: str) -> list[str]:
    ids: list[str] = []
    for gene in fms_iter_json_records(path):
        entity = gene.get("basicGeneticEntity") or {}
        primary = entity.get("primaryId", "")
        if primary.startswith(prefix):
            ids.append(primary)
            continue
        for x in entity.get("crossReferences") or []:
            xid = x.get("id", "")
            if xid.startswith(prefix):
                ids.append(xid)
                break
    return ids


def _file_gene_ids(path: Path, prefix: str) -> list[str]:
    """Gene CURIEs of a GENE-DESCRIPTION-JSON (gene_id) or BGI file."""
    ids = [r["gene_id"] for r in fms_iter_json_records(path)
           if str(r.get("gene_id", "")).startswith(prefix)]
    return ids or _bgi_gene_ids(path, prefix)


def _gene_list_path(release: str, mod: str) -> Path:
    return CACHE_DIR / "gene-ids" / release / f"{mod}.descriptions+bgi.txt"


def _mod_gene_ids(mod: str, release: str, urls: dict[str, list[str]]) -> list[str]:
    """Gene CURIEs of one MOD: the union of its GENE-DESCRIPTION-JSON and BGI
    files, cached under $ALLIANCE_FETCH_CACHE/gene-ids/<release>/ when the
    release is known."""
    path = _gene_list_path(release, mod)
    cacheable = release != "unknown"
    if cacheable and path.exists():
        ids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        log.info("Enumerated %d %s genes from %s", len(ids), mod, path)
        return ids
    if not urls:
        urls.update(_resolve_gene_list_urls(ALL_MODS))
    if mod not in urls:
        raise RuntimeError(f"No gene-list file for {mod} in FMS release {release}")
    ids: list[str] = []
    seen: set[str] = set()
    for url in urls[mod]:
        file_ids = _file_gene_ids(fms_download_gz(url), _MOD_PREFIX[mod])
        if not file_ids:
            raise RuntimeError(f"{url} yielded no {mod} gene IDs")
        added = [g for g in file_ids if g not in seen]
        seen.update(added)
        ids.extend(added)
        log.info("Enumerated %d %s genes (%d new) from FMS %s", len(file_ids), mod, len(added), url)
    if cacheable:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text("\n".join(ids) + "\n")
        os.replace(tmp, path)
    return ids


def enumerate_mod_genes(mods: list[str]) -> list[str]:
    """Return gene CURIEs across one or more MODs from each MOD's gene-list
    files for the pinned release. The per-MOD lists are cached so the
    fetchers in one run enumerate once instead of re-downloading the files."""
    release = get_current_release()
    urls: dict[str, list[str]] = {}
    combined: list[str] = []
    seen: set[str] = set()
    for mod in mods:
        for gid in _mod_gene_ids(mod, release, urls):
            if gid not in seen:
                seen.add(gid)
                combined.append(gid)
    log.info("Total across MODs %s: %d genes", ",".join(mods), len(combined))
    return combined


def mod_gene_taxa(mods: list[str] | None = None) -> dict[str, str]:
    """Gene CURIE -> NCBITaxon CURIE from the per-MOD gene lists. The
    9.2.0 API dropped gene.taxon from several payloads (/gene/{id}/models),
    and the curie prefix alone cannot tell X. laevis from X. tropicalis."""
    release = get_current_release()
    urls: dict[str, list[str]] = {}
    taxa: dict[str, str] = {}
    for mod in mods or ALL_MODS:
        for gid in _mod_gene_ids(mod, release, urls):
            taxa.setdefault(gid, MOD_TAXON[mod])
    return taxa


# ---------------------------------------------------------------------------
# FMS (file management system) bulk extraction
#
# AGR ships per-MOD bulk JSON dumps on a "release" cadence via FMS. Some
# datatypes are populated in FMS but empty on the live REST API (notably
# /gene/{id}/diseases-by-experiment as of 2026-05-11), so fetchers default to
# the FMS path when available and fall back to the API only when explicitly
# requested. FMS is itself slated for sunset; the long-term replacement is
# direct AGR Postgres access — keep fetcher code source-agnostic.


def _fms_file_list_url(release: str) -> str:
    return f"{FMS_BASE}/datafile/by/release/{release}?onlyCurrent=true"


def fms_list_files(
    datatype: str,
    *,
    release: str | None = None,
    sub_filter: str | None = None,
) -> list[dict]:
    """Return FMS file manifests for the given datatype (e.g. 'DAF',
    'CONSTRUCT', 'ALLELE', 'PHENOTYPE'). sub_filter narrows by dataSubType
    name (a MOD prefix like 'MGI' or a NCBITaxon ID like 'NCBITaxon10090').

    A release manifest lists every upload made during the release cycle, so
    one subtype can appear several times (e.g. five ZFIN BGI files in 9.1.0).
    Only the most recent upload per subtype is returned.
    """
    if release is None:
        release = get_current_release()
    if release == "unknown":
        raise RuntimeError("Cannot list FMS files: release is unknown")
    latest: dict[str, dict] = {}
    for row in _fms_manifest(release):
        dt = (row.get("dataType") or {}).get("name", "")
        if dt != datatype:
            continue
        sub = (row.get("dataSubType") or {}).get("name", "")
        if sub_filter and sub != sub_filter:
            continue
        if sub in latest and latest[sub].get("uploadDate", "") >= row.get("uploadDate", ""):
            continue
        latest[sub] = row
    out: list[dict] = []
    for sub, row in sorted(latest.items()):
        url = row.get("s3Url") or row.get("s3Path") or ""
        if not url.startswith("http"):
            url = "https://" + url.lstrip("/")
        out.append({
            "url": url,
            "datatype": datatype,
            "subType": sub,
            "release": release,
            "uploadDate": row.get("uploadDate", ""),
        })
    return out


@functools.lru_cache(maxsize=None)
def _fms_manifest(release: str) -> list[dict]:
    return http_get_json(_fms_file_list_url(release)) or []


def fms_tsv_column(datatype: str, column: str, *, subtypes: list[str] | None = None) -> list[str]:
    """Distinct values of one column across the release's FMS TSV files of
    a datatype (optionally only the given dataSubTypes), in file order."""
    import gzip
    seen: set[str] = set()
    out: list[str] = []
    files = fms_list_files(datatype)
    if subtypes is not None:
        by_sub = {fi["subType"]: fi for fi in files}
        files = [by_sub[sub] for sub in subtypes if sub in by_sub]
    for fi in files:
        local = fms_download_gz(fi["url"])
        idx = None
        with gzip.open(local, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if idx is None:
                    if column not in cols:
                        raise RuntimeError(f"{fi['url']} has no {column} column")
                    idx = cols.index(column)
                    continue
                if len(cols) > idx:
                    value = cols[idx].strip()
                    if value and value not in seen:
                        seen.add(value)
                        out.append(value)
        log.info("%s %s: %d distinct %s so far", datatype, fi["subType"], len(out), column)
    return out


def enumerate_disease_ids() -> list[str]:
    """Distinct DOIDs of the pinned release's DISEASE-ALLIANCE COMBINED file."""
    return fms_tsv_column("DISEASE-ALLIANCE", "DOID", subtypes=["COMBINED"])


def enumerate_allele_ids(mods: list[str]) -> list[str]:
    """Distinct allele CURIEs of the pinned release's per-taxon VARIANT-ALLELE
    TSVs. The COMBINED file is frozen at 4.0.0 and is not used."""
    subtypes = [MOD_TAXON[m].replace(":", "") for m in mods]
    return fms_tsv_column("VARIANT-ALLELE", "AlleleId", subtypes=list(dict.fromkeys(subtypes)))


def fms_download_gz(url: str) -> Path:
    """Download a gzipped FMS payload to FMS_CACHE_DIR keyed on the URL.
    Returns the local cached path. Idempotent — uses HEAD ETag matching
    when possible, otherwise size-stat. Caller deals with gunzip + JSON
    parsing (use fms_iter_json_records)."""
    import urllib.request as _u
    FMS_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(url.encode("utf-8")).hexdigest()[:24]
    suffix = ".json.gz" if url.endswith(".gz") else ".bin"
    target = FMS_CACHE_DIR / f"{key}{suffix}"
    if target.exists() and target.stat().st_size > 0:
        log.debug("FMS cache HIT %s -> %s", url, target.name)
        return target
    log.info("FMS download %s", url)
    req = _u.Request(url, headers={"User-Agent": USER_AGENT})
    part = target.with_suffix(target.suffix + ".part")
    try:
        with _u.urlopen(req, timeout=300) as resp, part.open("wb") as fh:
            while True:
                chunk = resp.read(1 << 16)
                if not chunk:
                    break
                fh.write(chunk)
        os.replace(part, target)
    finally:
        part.unlink(missing_ok=True)
    return target


def fms_iter_json_records(
    path: Path,
    *,
    key: str = "data",
) -> Iterator[dict]:
    """Stream records from a downloaded FMS JSON.gz file. Most AGR FMS files
    are wrapped {metaData, data: [...]} with a single top-level array under
    the 'data' key. Loaded in-memory because gzip is small enough (<200MB
    decompressed) for the volumes involved — switch to ijson if a single
    file grows past memory budget."""
    import gzip
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        try:
            doc = json.load(fh)
        except json.JSONDecodeError as e:
            log.error("FMS file %s parse error: %s", path, e)
            return
    records = doc.get(key) if isinstance(doc, dict) else doc
    if not isinstance(records, list):
        log.warning("FMS file %s has no '%s' array", path, key)
        return
    for r in records:
        if isinstance(r, dict):
            yield r


def iter_fms_datatype(
    datatype: str,
    *,
    release: str | None = None,
    sub_filter: str | None = None,
) -> Iterator[tuple[str, dict]]:
    """High-level: enumerate every record across every shard of an FMS
    datatype. Yields (subType, record) tuples so the caller can tag rows
    with the MOD without re-parsing the URL."""
    files = fms_list_files(datatype, release=release, sub_filter=sub_filter)
    if not files:
        raise RuntimeError(
            f"No FMS files for datatype={datatype} sub_filter={sub_filter} in release "
            f"{release or get_current_release()}; narrow --mods or ALLIANCE_MODS to skip it"
        )
    for fi in files:
        local = fms_download_gz(fi["url"])
        for rec in fms_iter_json_records(local):
            yield fi["subType"], rec


# ---------------------------------------------------------------------------
# Bounded-parallel per-entity runner for the API fan-out fetchers

DEFAULT_WORKERS = int(os.environ.get("ALLIANCE_FETCH_WORKERS", "6"))
DEFAULT_MAX_FAILURE_RATE = float(os.environ.get("ALLIANCE_MAX_FAILURE_RATE", "0.005"))


def add_parallel_args(parser) -> None:
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Concurrent API requests (default {DEFAULT_WORKERS}, env ALLIANCE_FETCH_WORKERS); "
             "the AGR WAF answers 405 well before 20",
    )
    parser.add_argument(
        "--max-failure-rate",
        type=float,
        default=DEFAULT_MAX_FAILURE_RATE,
        help="Exit non-zero when more than this fraction of entities fail "
             f"(default {DEFAULT_MAX_FAILURE_RATE}, env ALLIANCE_MAX_FAILURE_RATE)",
    )


@dataclass
class RunReport:
    label: str
    total: int
    done: int = 0
    failures: list[tuple[str, str]] = field(default_factory=list)
    elapsed: float = 0.0

    @property
    def failure_rate(self) -> float:
        return len(self.failures) / self.total if self.total else 0.0


def run_parallel(
    ids: list[str],
    work: Callable[[str, sqlite3.Connection, FetchStats], Any],
    consume: Callable[[str, Any], None],
    *,
    cache_name: str,
    stats: FetchStats,
    workers: int = DEFAULT_WORKERS,
    label: str = "genes",
    progress_every: int = 1000,
    progress_seconds: float = 60.0,
) -> RunReport:
    """Run work(id, cache, stats) for every id on a bounded thread pool.

    consume(id, result) runs on the calling thread, in input order, so TSV
    writing and de-duplication stay single-threaded and the output is
    deterministic. Each worker thread owns one autocommit cache connection.
    An exception in work() marks that id failed and the run carries on;
    the caller decides from the report whether the run as a whole failed."""
    report = RunReport(label=label, total=len(ids))
    local = threading.local()
    conns: list[sqlite3.Connection] = []
    conns_lock = threading.Lock()

    def thread_cache() -> sqlite3.Connection:
        conn = getattr(local, "conn", None)
        if conn is None:
            conn = connect_cache(cache_name, initialise=False)
            local.conn = conn
            with conns_lock:
                conns.append(conn)
        return conn

    def task(entity_id: str) -> tuple[str, Any]:
        try:
            return "ok", work(entity_id, thread_cache(), stats)
        except urllib.error.HTTPError as e:
            return "failed", f"HTTP {e.code}"
        except Exception as e:
            return "failed", f"{type(e).__name__}: {e}"

    connect_cache(cache_name).close()
    start = time.time()
    last_log = start
    workers = max(1, workers)
    log.info("%s: %d to fetch with %d workers", label, len(ids), workers)
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pending: deque = deque()
            source = iter(ids)

            def submit_next() -> None:
                entity_id = next(source, None)
                if entity_id is not None:
                    pending.append((entity_id, pool.submit(task, entity_id)))

            for _ in range(workers * 4):
                submit_next()
            while pending:
                entity_id, fut = pending.popleft()
                status, result = fut.result()
                submit_next()
                report.done += 1
                if status == "ok":
                    consume(entity_id, result)
                else:
                    report.failures.append((entity_id, result))
                    log.warning("%s: %s failed: %s", label, entity_id, result)
                now = time.time()
                if (report.done % progress_every == 0 or now - last_log >= progress_seconds
                        or report.done == report.total):
                    _log_progress(report, stats, now - start)
                    last_log = now
    finally:
        for conn in conns:
            conn.close()
    report.elapsed = time.time() - start
    return report


def _log_progress(report: RunReport, stats: FetchStats, elapsed: float) -> None:
    rate = report.done / elapsed if elapsed > 0 else 0.0
    remaining = (report.total - report.done) / rate if rate > 0 else 0.0
    log.info(
        "%s: %d/%d (%.1f%%)  %.1f %s/s  %.1f req/s  cache_hits=%d  failed=%d  "
        "throttled=%d  eta %.1fh",
        report.label, report.done, report.total,
        100.0 * report.done / report.total if report.total else 100.0,
        rate, report.label, stats.requests / elapsed if elapsed > 0 else 0.0,
        stats.cache_hits, len(report.failures), stats.throttled,
        remaining / 3600,
    )


def finish_run(
    report: RunReport,
    stats: FetchStats,
    out_dir: Path,
    max_failure_rate: float,
    *,
    rows_written: int,
) -> None:
    """Log the run summary, write failed ids next to the TSVs, and exit
    non-zero when the failure rate is above max_failure_rate or no row was
    written. Call it inside the TsvWriter block: the SystemExit then makes the
    writer drop its temp file instead of renaming a partial TSV into place.
    Failed ids are not cached, so a rerun retries exactly those (everything
    else is served from the cache)."""
    log.info(
        "%s: done %d/%d in %.0fs  failed=%d (%.3f%%)  rows=%d  API stats: %s",
        report.label, report.done, report.total, report.elapsed,
        len(report.failures), 100.0 * report.failure_rate, rows_written, stats,
    )
    if stats.cache_hits:
        log.info("%s: %d responses served from the cache", report.label, stats.cache_hits)
    failures_path = Path(out_dir) / f"{Path(sys.argv[0]).stem}.failures.tsv"
    if report.failures:
        failures_path.parent.mkdir(parents=True, exist_ok=True)
        with failures_path.open("w", encoding="utf-8") as fh:
            for entity_id, err in report.failures:
                fh.write(f"{entity_id}\tfailed\t{err}\n")
        for entity_id, err in report.failures[:20]:
            log.warning("  failed %s: %s", entity_id, err)
        log.warning("Failed ids written to %s", failures_path)
    elif failures_path.exists():
        failures_path.unlink()
    if report.failure_rate > max_failure_rate:
        log.error(
            "%s: failure rate %.3f%% exceeds %.3f%%",
            report.label, 100.0 * report.failure_rate, 100.0 * max_failure_rate,
        )
        sys.exit(1)
    if rows_written == 0:
        log.error("%s: no rows written", report.label)
        sys.exit(1)


def rows_from_crossrefs(path: Path) -> Iterator[dict]:
    """Yield {uniProtAccession, geneId, geneTaxon} per UniProtKB xref in gene-crossrefs.tsv."""
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0] == "geneId":
                continue
            if len(cols) < 5:
                continue
            gene_id, taxon, xref_id = cols[:3]
            if not xref_id.startswith("UniProtKB:"):
                continue
            yield {
                "uniProtAccession": xref_id.split(":", 1)[1],
                "geneId": gene_id,
                "geneTaxon": taxon,
            }
