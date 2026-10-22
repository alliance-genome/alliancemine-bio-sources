#!/usr/bin/env python3
"""
Build orthologs.tsv from a restored Alliance curation database.

Same file, columns and value formats as fetch_orthologs.py. The API's
/gene/{id}/orthologs endpoint returns the generated orthology rows whose
subject is the seed gene and whose strictFilter is true, so this script
selects those rows and keeps the ones whose subject is a seed gene. Rows
with an obsolete gene on either side are dropped, as the API drops them. Method
lists keep the join-table row order, which is the order the API returns.

API fallback: the 9.2.0 curation dump holds no orthology rows for SGD
genes. Seeds whose taxon is the subject of no database row are fetched from
the API as fetch_orthologs.py does, and each fetched row is mirrored so the
other species' genes also list their SGD orthologs (the database stores
every pair in both directions with the same methods). --no-api-fallback
turns this off.

Usage:
    python3 extract_orthologs.py                   # every gene of the 9 MODs
    python3 extract_orthologs.py --mods MGI,HUMAN
    python3 extract_orthologs.py --ids MGI:101758
    python3 extract_orthologs.py --limit 200

Output: data/orthologs.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Callable

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    mod_gene_taxa,
    paginate,
    resolve_gene_ids,
    run_parallel,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402
import fetch_orthologs  # noqa: E402

log = logging.getLogger("alliance.extract.orthologs")

# The API labels every row it returns by default "stringent": the rows with
# strictFilter = true.
STRINGENCY = "stringent"

_GENES = """
    select b.id, b.primaryexternalid, t.curie
    from gene g
    join biologicalentity b on b.id = g.id
    left join ontologyterm t on t.id = b.taxon_id
    where not coalesce(b.obsolete, false)
"""


def _methods_cte(name: str, kind: str) -> str:
    table = f"genetogeneorthologygenerated_predictionmethods{kind}"
    return f"""
    {name} as (
        select j.genetogeneorthologygenerated_id id,
               string_agg(v.name, '|' order by j.ctid) names
        from {table} j
        join vocabularyterm v on v.id = j.predictionmethods{kind}_id
        where j.genetogeneorthologygenerated_id in (select id from strict)
        group by 1
    )"""


SQL = f"""
    with strict as (
        select o.id, o.subjectgene_id, o.objectgene_id
        from genetogeneorthology o
        join genetogeneorthologygenerated g on g.id = o.id
        where g.strictfilter and not coalesce(o.obsolete, false)
    ),
    genes as ({_GENES}),
    {_methods_cte("matched", "matched")},
    {_methods_cte("not_matched", "notmatched")},
    {_methods_cte("not_called", "notcalled")}
    select s.primaryexternalid, s.curie, b.primaryexternalid, b.curie,
           m.names, nm.names, nc.names
    from strict o
    join genes s on s.id = o.subjectgene_id
    join genes b on b.id = o.objectgene_id
    left join matched m on m.id = o.id
    left join not_matched nm on nm.id = o.id
    left join not_called nc on nc.id = o.id
"""


def to_row(values: list[str]) -> dict:
    subject, subject_taxon, ortholog, ortholog_taxon, matched, not_matched, not_called = values
    return {
        "subjectGeneId": subject,
        "subjectGeneTaxon": subject_taxon,
        "orthologGeneId": ortholog,
        "orthologGeneTaxon": ortholog_taxon,
        "stringencyFilter": STRINGENCY,
        "predictionMethodsMatched": matched,
        "predictionMethodsNotMatched": not_matched,
        "predictionMethodsNotCalled": not_called,
    }


def mirror(row: dict) -> dict:
    """The same orthology pair seen from the ortholog's side."""
    return {
        **row,
        "subjectGeneId": row["orthologGeneId"],
        "subjectGeneTaxon": row["orthologGeneTaxon"],
        "orthologGeneId": row["subjectGeneId"],
        "orthologGeneTaxon": row["subjectGeneTaxon"],
    }


def seeded_rows(rows, seed_rank: dict[str, int]) -> list[dict]:
    """Rows whose subject is a seed gene, one per ordered gene pair, in seed
    order (the order fetch_orthologs.py writes them)."""
    kept: dict[tuple[str, str], dict] = {}
    for row in rows:
        pair = (row["subjectGeneId"], row["orthologGeneId"])
        if all(pair) and pair[0] in seed_rank and pair not in kept:
            kept[pair] = row
    return sorted(kept.values(), key=lambda r: (seed_rank[r["subjectGeneId"]], r["orthologGeneId"]))


def gap_seeds(gene_ids: list[str], seed_taxa: dict[str, str], covered_taxa: set[str]) -> list[str]:
    """Seeds whose taxon is the subject of no database row."""
    return [g for g in gene_ids if g in seed_taxa and seed_taxa[g] not in covered_taxa]


def api_rows(seeds: list[str], endpoint: str, to_api_row: Callable[[str, dict], dict],
             args: argparse.Namespace) -> list[dict]:
    """Rows for seeds from the API endpoint /gene/{id}/<endpoint>, through the
    fetcher's own row mapping and response cache."""
    stats = FetchStats()
    rows: list[dict] = []

    def work(gid, cache, stats):
        return [to_api_row(gid, r) for r in paginate(f"/gene/{gid}/{endpoint}", cache=cache, stats=stats)]

    def consume(_gid, result):
        rows.extend(result)

    report = run_parallel(seeds, work, consume, cache_name=endpoint, stats=stats,
                          workers=args.workers, label=f"{endpoint} api-fallback genes")
    finish_run(report, stats, Path(args.out_dir), args.max_failure_rate, rows_written=len(rows) or 1)
    return rows


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--mods", help="Comma-separated MODs to seed from (default: all 9, env ALLIANCE_MODS)")
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--no-api-fallback", action="store_true",
                        help="Skip the API for seeds whose species has no rows in the database")
    add_parallel_args(parser)
    parser.add_argument("--verbose", action="store_true")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance ortholog per-algorithm detail from the curation DB")
    add_args(parser)
    args = parser.parse_args()

    configure_logging(args.verbose)
    started = time.monotonic()

    release = get_current_release()
    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit)
    seed_rank = {gid: i for i, gid in enumerate(gene_ids)}
    rows = [to_row(values) for values in copy_rows(SQL)]
    log.info("orthologs: %d strict rows in the database", len(rows))

    covered = {r["subjectGeneTaxon"] for r in rows}
    gap = gap_seeds(gene_ids, mod_gene_taxa(selected_mods(None, args.mods)), covered)
    if gap and not args.no_api_fallback:
        log.info("orthologs: %d seeds of taxa with no database rows go to the API", len(gap))
        fetched = api_rows(gap, "orthologs", fetch_orthologs._row, args)
        rows.extend(fetched)
        rows.extend(mirror(r) for r in fetched if r["orthologGeneTaxon"] in covered)
        log.info("orthologs: %d API rows, mirrored for database taxa", len(fetched))
    elif gap:
        log.warning("orthologs: %d seeds of taxa with no database rows skipped", len(gap))

    rows = seeded_rows(rows, seed_rank)
    with TsvWriter(
        Path(args.out_dir) / "orthologs.tsv",
        fetch_orthologs.COLUMNS,
        release=release,
        source="curation DB genetogeneorthologygenerated (strictFilter)",
    ) as writer:
        for row in rows:
            writer.write_row(**row)

    seeds_with_rows = len({r["subjectGeneId"] for r in rows})
    log.info("orthologs: rows=%d seeds=%d seeds_with_rows=%d elapsed=%.0fs",
             len(rows), len(gene_ids), seeds_with_rows, time.monotonic() - started)


if __name__ == "__main__":
    main()
