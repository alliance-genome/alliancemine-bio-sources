"""
Read access to a restored Alliance curation database dump.

The extract_*.py scripts build the same TSVs as the fetch_*.py API fetchers
from SQL against a pg_restore of the curation database. Rows stream through
psql's COPY ... TO STDOUT in CSV form, so the scripts need only a psql client
and no Python driver.

Environment variables:
    CURATION_PSQL    - psql command line, split with shlex; default
                       'psql -X -q' plus CURATION_DB_URL. Example for a local
                       container: 'docker exec -i agr-curation-920 psql -X -q -U postgres -d curation'
    CURATION_DB_URL  - libpq connection string used by the default command
"""

from __future__ import annotations

import csv
import io
import logging
import os
import shlex
import subprocess
import sys
from typing import Iterator

log = logging.getLogger("alliance.curation_db")

csv.field_size_limit(sys.maxsize)


def psql_command() -> list[str]:
    explicit = os.environ.get("CURATION_PSQL", "").strip()
    if explicit:
        return shlex.split(explicit)
    url = os.environ.get("CURATION_DB_URL", "").strip()
    if not url:
        raise SystemExit("Set CURATION_PSQL or CURATION_DB_URL to reach the curation database")
    return ["psql", "-X", "-q", url]


def copy_rows(sql: str) -> Iterator[list[str]]:
    """Yield each result row of sql as a list of strings; NULL becomes ''."""
    copy = f"COPY ({sql.strip().rstrip(';')}) TO STDOUT WITH (FORMAT csv)"
    cmd = psql_command() + ["-v", "ON_ERROR_STOP=1", "-c", copy]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None and proc.stderr is not None
    # Popen takes no newline= before Python 3.10; csv needs newline="" so
    # line breaks inside quoted values survive.
    stdout = io.TextIOWrapper(proc.stdout, encoding="utf-8", newline="")
    try:
        yield from csv.reader(stdout)
    except GeneratorExit:
        # Caller stopped early (--limit): drop the rest of the result.
        proc.kill()
        proc.wait()
        raise
    err = proc.stderr.read().decode("utf-8", "replace")
    if proc.wait() != 0:
        raise RuntimeError(f"psql failed (rc={proc.returncode}): {err.strip()}")


def scalar(sql: str) -> str:
    rows = list(copy_rows(sql))
    return rows[0][0] if rows and rows[0] else ""
