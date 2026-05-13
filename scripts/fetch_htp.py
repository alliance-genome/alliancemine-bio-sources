#!/usr/bin/env python3
"""
Fetch high-throughput experiment metadata from the Alliance FMS
HTPDATASET + HTPDATASAMPLE bulk files. Maps to existing DataSet +
Sample + SampleCharacteristic classes in the genomic model.

HTPDATASET rows describe the parent study (title, summary, GEO/SRA
cross-reference, publications). HTPDATASAMPLE rows describe individual
sequencing/array samples within a study (anatomy, stage, sex, assay
type, biosample id).

Two TSVs:
    data/htp-datasets.tsv
    data/htp-samples.tsv

Usage:
    python3 fetch_htp.py
    python3 fetch_htp.py --mods MGI,ZFIN
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
    get_nested,
    iter_fms_datatype,
    join_pipe,
    normalize_taxon,
)

log = logging.getLogger("alliance.fetch.htp")


DEFAULT_MODS = ["MGI", "ZFIN", "FB", "WB", "RGD", "SGD", "XBXL", "XBXT"]


DATASET_COLUMNS = [
    "datasetId",
    "title",
    "summary",
    "categoryTags",
    "geoXref",
    "publicationPmids",
    "dateAssigned",
    "dataProvider",
]

SAMPLE_COLUMNS = [
    "sampleId",
    "sampleTitle",
    "abundance",
    "sampleType",          # OBI term curie
    "stageId",             # ZFS/FBdv/MmusDv/etc term
    "stageName",
    "anatomyId",
    "anatomyStatement",
    "biosampleId",
    "sex",
    "assayType",           # MMO term curie
    "assemblyVersion",
    "datasetIds",          # pipe-joined parent dataset ids
    "taxon",
    "dataProvider",
]


def _ds_row(rec: dict, mod: str) -> dict | None:
    ds = rec.get("datasetId") or {}
    primary = ds.get("primaryId")
    if not primary:
        return None
    geo = ""
    for x in (ds.get("alternateIds") or []):
        if isinstance(x, str) and x.startswith("GEO:"):
            geo = x
            break
    pmids = join_pipe(
        p.get("publicationId", "") for p in (rec.get("publications") or []) if p.get("publicationId")
    )
    return {
        "datasetId": primary,
        "title": rec.get("title") or "",
        "summary": rec.get("summary") or "",
        "categoryTags": join_pipe(rec.get("categoryTags") or []),
        "geoXref": geo,
        "publicationPmids": pmids,
        "dateAssigned": rec.get("dateAssigned") or "",
        "dataProvider": mod,
    }


def _sample_row(rec: dict, mod: str) -> dict | None:
    sid = (rec.get("sampleId") or {}).get("primaryId")
    if not sid:
        return None
    stage = (rec.get("sampleAge") or {}).get("stage") or {}
    locations = rec.get("sampleLocations") or []
    first_loc = locations[0] if locations else {}
    return {
        "sampleId": sid,
        "sampleTitle": rec.get("sampleTitle") or "",
        "abundance": rec.get("abundance") or "",
        "sampleType": rec.get("sampleType") or "",
        "stageId": stage.get("stageTermId") or "",
        "stageName": stage.get("stageName") or "",
        "anatomyId": first_loc.get("anatomicalStructureTermId") or "",
        "anatomyStatement": first_loc.get("whereExpressedStatement") or "",
        "biosampleId": get_nested(rec, "genomicInformation", "biosampleId", default=""),
        "sex": rec.get("sex") or "",
        "assayType": rec.get("assayType") or "",
        "assemblyVersion": (rec.get("assemblyVersions") or [""])[0],
        "datasetIds": join_pipe(rec.get("datasetIds") or []),
        "taxon": normalize_taxon(rec.get("taxonId") or ""),
        "dataProvider": mod,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance HTP datasets/samples (FMS bulk)")
    parser.add_argument(
        "--mods",
        default=",".join(DEFAULT_MODS),
        help="Comma-separated MOD prefixes",
    )
    parser.add_argument("--limit", type=int, help="Cap rows emitted (per TSV)")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for emitted TSVs",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    mods = [m.strip() for m in args.mods.split(",") if m.strip()]
    release = get_current_release()
    out_dir = Path(args.out_dir)

    ds_seen: set = set()
    ds_total = 0
    with TsvWriter(
        out_dir / "htp-datasets.tsv",
        DATASET_COLUMNS,
        release=release,
        source=f"FMS HTPDATASET ({','.join(mods)})",
    ) as writer:
        for mod in mods:
            log.info("HTPDATASET FMS: %s", mod)
            for _sub, rec in iter_fms_datatype("HTPDATASET", sub_filter=mod):
                row = _ds_row(rec, mod)
                if row is None or row["datasetId"] in ds_seen:
                    continue
                ds_seen.add(row["datasetId"])
                writer.write_row(**row)
                ds_total += 1
                if args.limit is not None and ds_total >= args.limit:
                    break
            if args.limit is not None and ds_total >= args.limit:
                break
    log.info("Datasets done. rows=%d", ds_total)

    s_seen: set = set()
    s_total = 0
    with TsvWriter(
        out_dir / "htp-samples.tsv",
        SAMPLE_COLUMNS,
        release=release,
        source=f"FMS HTPDATASAMPLE ({','.join(mods)})",
    ) as writer:
        for mod in mods:
            log.info("HTPDATASAMPLE FMS: %s", mod)
            for _sub, rec in iter_fms_datatype("HTPDATASAMPLE", sub_filter=mod):
                row = _sample_row(rec, mod)
                if row is None or row["sampleId"] in s_seen:
                    continue
                s_seen.add(row["sampleId"])
                writer.write_row(**row)
                s_total += 1
                if args.limit is not None and s_total >= args.limit:
                    break
            if args.limit is not None and s_total >= args.limit:
                break
    log.info("Samples done. rows=%d", s_total)


if __name__ == "__main__":
    main()
