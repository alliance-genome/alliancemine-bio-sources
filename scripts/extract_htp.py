#!/usr/bin/env python3
"""
Build data/htp-datasets.tsv and data/htp-samples.tsv from a restored
Alliance curation database, with the columns fetch_htp.py writes from the
FMS HTPDATASET and HTPDATASAMPLE files.

Datasets come from htpexpressiondatasetannotation, samples from
htpexpressiondatasetsampleannotation; both skip obsolete and internal rows.
A dataset or sample id seen under an earlier MOD is not written again,
as in the fetcher.

Differences from the FMS files, all in columns AllianceHtpConverter reads
only for SampleCharacteristic values:
  - The curation DB keeps no sample whereExpressedStatement, so
    anatomyStatement carries the name of the anatomy term.
  - dateAssigned is the annotation's datecreated (ISO 8601, UTC), blank
    where the DB has none.
Samples without a sample id (all MGI samples in the 9.2.0 dump) are
skipped, as the fetcher skips FMS records without sampleId.primaryId.

Publication ids resolve through reference_crossreference.

Usage:
    python3 extract_htp.py
    python3 extract_htp.py --mods MGI,ZFIN
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    MOD_TAXON,
    TsvWriter,
    configure_logging,
    get_current_release,
    normalize_taxon,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402
import extract_expression  # noqa: E402
from extract_expression import pick_reference  # noqa: E402
from fetch_htp import DATASET_COLUMNS, DEFAULT_MODS, SAMPLE_COLUMNS  # noqa: E402

log = logging.getLogger("alliance.extract.htp")

# organization.abbreviation of each MOD's data provider.
PROVIDER = {"SGD": "SGD", "MGI": "MGI", "RGD": "RGD", "ZFIN": "ZFIN", "FB": "FB",
            "WB": "WB", "XBXL": "XB", "XBXT": "XB"}

# Genetic sex vocabulary names that the FMS files spell differently.
SEX_LABEL = {"unknown sex": "unknown", "pooled sexes": "pooled"}

LIVE = "NOT coalesce({t}.obsolete, false) AND NOT coalesce({t}.internal, false)"


def sex_label(name: str) -> str:
    return SEX_LABEL.get(name, name)


def sql_list(values: list[str]) -> str:
    return ", ".join(f"'{v}'" for v in values)


def dataset_filter_sql(providers: list[str]) -> str:
    return f"""
        SELECT d.* FROM htpexpressiondatasetannotation d
        JOIN organization o ON o.id = d.dataprovider_id
        WHERE {LIVE.format(t='d')} AND o.abbreviation IN ({sql_list(providers)})"""


def dataset_sql(providers: list[str]) -> str:
    return f"""
        WITH ds AS ({dataset_filter_sql(providers)}),
        tags AS (
            SELECT c.htpexpressiondatasetannotation_id AS id, string_agg(v.name, '|' ORDER BY v.name) AS tags
            FROM htpexpressiondatasetannotation_categorytags c
            JOIN vocabularyterm v ON v.id = c.categorytags_id
            GROUP BY 1
        ),
        geo AS (
            SELECT externaldatabaseentity_id AS id, min(secondaryidentifiers) AS geo
            FROM externaldatabaseentity_secondaryidentifiers
            WHERE secondaryidentifiers LIKE 'GEO:%'
            GROUP BY 1
        ),
        refs AS (
            SELECT htpexpressiondatasetannotation_id AS id, string_agg(references_id::text, '|') AS refs
            FROM htpexpressiondatasetannotation_reference
            GROUP BY 1
        )
        SELECT e.curie, coalesce(d.name, ''), coalesce(n.freetext, ''), coalesce(tags.tags, ''),
               coalesce(geo.geo, ''), coalesce(refs.refs, ''),
               coalesce(to_char(d.datecreated AT TIME ZONE 'UTC', 'YYYY-MM-DD"T"HH24:MI:SS"+00:00"'), ''),
               o.abbreviation
        FROM ds d
        JOIN externaldatabaseentity e ON e.id = d.htpexpressiondataset_id
        JOIN organization o ON o.id = d.dataprovider_id
        LEFT JOIN note n ON n.id = d.relatednote_id AND NOT coalesce(n.obsolete, false)
        LEFT JOIN tags ON tags.id = d.id
        LEFT JOIN geo ON geo.id = e.id
        LEFT JOIN refs ON refs.id = d.id
        WHERE coalesce(e.curie, '') <> ''
        ORDER BY array_position(ARRAY[{sql_list(providers)}]::text[], o.abbreviation::text), e.curie"""


def sample_sql(taxa: list[str]) -> str:
    return f"""
        WITH loc AS (
            SELECT DISTINCT ON (l.htpexpressiondatasetsampleannotation_id)
                   l.htpexpressiondatasetsampleannotation_id AS id, t.curie, t.name
            FROM htpexpressiondatasetsampleannotation_anatomicalsite l
            JOIN anatomicalsite a ON a.id = l.htpexpressionsamplelocations_id
            LEFT JOIN ontologyterm t ON t.id = a.anatomicalstructure_id
            ORDER BY l.htpexpressiondatasetsampleannotation_id, l.htpexpressionsamplelocations_id
        ),
        asm AS (
            SELECT htpexpressiondatasetsampleannotation_id AS id, min(assemblyversions) AS assembly
            FROM htpexpressiondatasetsampleannotation_assemblyversions
            GROUP BY 1
        ),
        dsl AS (
            SELECT x.htpexpressiondatasetsampleannotation_id AS id, string_agg(e.curie, '|' ORDER BY e.curie) AS ids
            FROM htpexpressiondatasetsampleannotation_externaldatabaseentity x
            JOIN externaldatabaseentity e ON e.id = x.datasetids_id
            GROUP BY 1
        )
        SELECT e.curie, coalesce(s.htpexpressionsampletitle, ''), coalesce(s.abundance, ''),
               coalesce(st.curie, ''), coalesce(stg.curie, ''), coalesce(ba.whenexpressedstagename, ''),
               coalesce(loc.curie, ''), coalesce(loc.name, ''),
               coalesce(agm.primaryexternalid, al.primaryexternalid, ''), coalesce(sex.name, ''),
               coalesce(mmo.curie, ''), coalesce(asm.assembly, ''), coalesce(dsl.ids, ''), tx.curie
        FROM htpexpressiondatasetsampleannotation s
        JOIN externaldatabaseentity e ON e.id = s.htpexpressionsample_id
        JOIN ontologyterm tx ON tx.id = s.taxon_id
        LEFT JOIN ontologyterm st ON st.id = s.htpexpressionsampletype_id
        LEFT JOIN ontologyterm mmo ON mmo.id = s.expressionassayused_id
        LEFT JOIN vocabularyterm sex ON sex.id = s.geneticsex_id
        LEFT JOIN biosampleage ba ON ba.id = s.htpexpressionsampleage_id
        LEFT JOIN temporalcontext tc ON tc.id = ba.stage_id
        LEFT JOIN ontologyterm stg ON stg.id = tc.developmentalstagestart_id
        LEFT JOIN biosamplegenomicinformation gi ON gi.id = s.genomicinformation_id
        LEFT JOIN biologicalentity agm ON agm.id = gi.biosampleagm_id
        LEFT JOIN biologicalentity al ON al.id = gi.biosampleallele_id
        LEFT JOIN loc ON loc.id = s.id
        LEFT JOIN asm ON asm.id = s.id
        LEFT JOIN dsl ON dsl.id = s.id
        WHERE {LIVE.format(t='s')} AND tx.curie IN ({sql_list(taxa)}) AND coalesce(e.curie, '') <> ''
        ORDER BY array_position(ARRAY[{sql_list(taxa)}]::text[], tx.curie::text), e.curie"""


def dataset_row(rec: list[str], refs: dict[str, list[str]], mod: str) -> dict:
    ds_id, title, summary, tags, geo, ref_ids, date, provider = rec
    pubs: list[str] = []
    for ref_id in sorted(r for r in ref_ids.split("|") if r):
        pub = pick_reference(refs.get(ref_id, []), mod)
        if pub and pub not in pubs:
            pubs.append(pub)
    return {
        "datasetId": ds_id,
        "title": title,
        "summary": summary,
        "categoryTags": tags,
        "geoXref": geo,
        "publicationPmids": "|".join(pubs),
        "dateAssigned": date,
        "dataProvider": provider,
    }


def sample_row(rec: list[str], mod: str) -> dict:
    (sid, title, abundance, sample_type, stage_id, stage_name, anatomy_id, anatomy_name,
     biosample, sex, assay, assembly, dataset_ids, taxon) = rec
    return {
        "sampleId": sid,
        "sampleTitle": title,
        "abundance": abundance,
        "sampleType": sample_type,
        "stageId": stage_id,
        "stageName": stage_name,
        "anatomyId": anatomy_id,
        "anatomyStatement": anatomy_name,
        "biosampleId": biosample,
        "sex": sex_label(sex),
        "assayType": assay,
        "assemblyVersion": assembly,
        "datasetIds": dataset_ids,
        "taxon": normalize_taxon(taxon),
        "dataProvider": mod,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance HTP datasets/samples from the curation DB")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
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

    mods = [m for m in selected_mods(DEFAULT_MODS, args.mods) if m in PROVIDER]
    release = get_current_release()
    out_dir = Path(args.out_dir)
    providers = list(dict.fromkeys(PROVIDER[m] for m in mods))
    provider_mod = {PROVIDER[m]: m for m in reversed(mods)}
    started = time.monotonic()

    refs = extract_expression.reference_curies(
        "SELECT r.references_id FROM htpexpressiondatasetannotation_reference r "
        f"JOIN ({dataset_filter_sql(providers)}) d ON d.id = r.htpexpressiondatasetannotation_id")
    ds_seen: set = set()
    with TsvWriter(out_dir / "htp-datasets.tsv", DATASET_COLUMNS, release=release,
                   source=f"curation DB htpexpressiondatasetannotation ({','.join(mods)})") as writer:
        for rec in copy_rows(dataset_sql(providers)):
            if rec[0] in ds_seen:
                continue
            ds_seen.add(rec[0])
            writer.write_row(**dataset_row(rec, refs, provider_mod[rec[-1]]))
            if args.limit is not None and len(ds_seen) >= args.limit:
                break
    log.info("Datasets done. rows=%d (%.0fs)", len(ds_seen), time.monotonic() - started)

    by_taxon = {MOD_TAXON[m]: m for m in reversed(mods)}
    s_seen: set = set()
    with TsvWriter(out_dir / "htp-samples.tsv", SAMPLE_COLUMNS, release=release,
                   source=f"curation DB htpexpressiondatasetsampleannotation ({','.join(mods)})") as writer:
        for rec in copy_rows(sample_sql(list(by_taxon))):
            if rec[0] in s_seen:
                continue
            s_seen.add(rec[0])
            writer.write_row(**sample_row(rec, by_taxon[rec[-1]]))
            if args.limit is not None and len(s_seen) >= args.limit:
                break
    log.info("Samples done. rows=%d (%.0fs)", len(s_seen), time.monotonic() - started)


if __name__ == "__main__":
    main()
