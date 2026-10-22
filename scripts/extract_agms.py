#!/usr/bin/env python3
"""
Build agms.tsv (same columns as fetch_agms.py) from a restored Alliance
curation database instead of the FMS AGM files.

modelName is the AGM full-name slot annotation (display text, the FMS `name`),
modelSubtype the subtype vocabulary term name, taxon the NCBITaxon curie and
dataProvider the curation data provider (Xenbase split by taxon).

Usage:
    python3 extract_agms.py
    python3 extract_agms.py --mods MGI,RGD --limit 50

Output: data/agms.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    get_current_release,
    normalize_taxon,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402
from extract_allele_phenotypes import mod_for, providers_sql  # noqa: E402
from fetch_agms import COLUMNS, DEFAULT_MODS  # noqa: E402

log = logging.getLogger("alliance.extract.agms")

AGMS_SQL = """
SELECT org.abbreviation, be.primaryexternalid, coalesce(fn.displaytext, ''),
       coalesce(st.name, ''), coalesce(tx.curie, '')
FROM affectedgenomicmodel a
JOIN biologicalentity be ON be.id = a.id
JOIN organization org ON org.id = be.dataprovider_id
LEFT JOIN vocabularyterm st ON st.id = a.subtype_id
LEFT JOIN ontologyterm tx ON tx.id = be.taxon_id
LEFT JOIN slotannotation fn ON fn.singleagm_id = a.id
     AND fn.slotannotationtype = 'AgmFullNameSlotAnnotation' AND NOT fn.obsolete
WHERE NOT be.obsolete AND NOT be.internal
  AND coalesce(be.primaryexternalid, '') <> ''
  AND org.abbreviation = ANY(ARRAY[{providers}])
ORDER BY array_position(ARRAY[{providers}], org.abbreviation::text), tx.curie, be.primaryexternalid
"""


def iter_rows(mods: list[str]) -> Iterator[dict]:
    for provider, model_id, name, subtype, taxon in copy_rows(AGMS_SQL.format(providers=providers_sql(mods))):
        mod = mod_for(provider, taxon)
        if mod not in mods:
            continue
        yield {
            "modelId": model_id,
            "modelName": name,
            "modelSubtype": subtype,
            "taxon": normalize_taxon(taxon),
            "dataProvider": mod,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance AGMs from the curation DB")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    mods = selected_mods(DEFAULT_MODS, args.mods)
    release = get_current_release()
    seen: set = set()

    with TsvWriter(
        Path(args.out_dir) / "agms.tsv",
        COLUMNS,
        release=release,
        source=f"curation DB affectedgenomicmodel ({','.join(mods)})",
    ) as writer:
        for row in iter_rows(mods):
            if row["modelId"] in seen:
                continue
            seen.add(row["modelId"])
            writer.write_row(**row)
            if args.limit is not None and writer.rows_written >= args.limit:
                break
        if writer.rows_written == 0:
            raise SystemExit("agms: no rows extracted")

    log.info("Done. rows=%d", writer.rows_written)


if __name__ == "__main__":
    main()
