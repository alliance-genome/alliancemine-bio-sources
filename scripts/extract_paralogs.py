#!/usr/bin/env python3
"""
Build paralogs.tsv from a restored Alliance curation database.

Same file, columns and value formats as fetch_paralogs.py. The database
stores each paralog pair once per direction, and identity, similarity,
length and rank belong to the subject gene. fetch_paralogs.py keeps the
first direction it meets while walking the seed genes in order, so this
script keeps, per unordered pair, the row whose subject comes first in the
seed list. Method lists keep the join-table row order, which is the order
the API returns. Pairs with an obsolete gene stay: the paralogs endpoint
returns them (5,127 pairs in 9.2.0), unlike the orthologs endpoint.

API fallback: the 9.2.0 curation dump holds no paralogy rows for SGD genes.
Seeds whose taxon is the subject of no database row are fetched from the API
as fetch_paralogs.py does (paralogs pair genes of one species, so these rows
cover the whole gap). --no-api-fallback turns this off.

Usage:
    python3 extract_paralogs.py                    # every gene of the 9 MODs
    python3 extract_paralogs.py --mods FB,WB
    python3 extract_paralogs.py --ids SGD:S000004103
    python3 extract_paralogs.py --limit 200

Output: data/paralogs.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    TsvWriter,
    configure_logging,
    get_current_release,
    mod_gene_taxa,
    resolve_gene_ids,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402
from extract_orthologs import add_args, api_rows, gap_seeds  # noqa: E402
import fetch_paralogs  # noqa: E402

log = logging.getLogger("alliance.extract.paralogs")

_GENES = """
    select b.id, b.primaryexternalid, t.curie, sym.displaytext symbol
    from gene g
    join biologicalentity b on b.id = g.id
    left join ontologyterm t on t.id = b.taxon_id
    left join slotannotation sym
           on sym.singlegene_id = g.id
          and sym.slotannotationtype = 'GeneSymbolSlotAnnotation'
          and not coalesce(sym.obsolete, false)
"""


def _methods_cte(name: str, kind: str) -> str:
    return f"""
    {name} as (
        select j.genetogeneparalogy_id id,
               string_agg(v.name, '|' order by j.ctid) names
        from genetogeneparalogy_predictionmethods{kind} j
        join vocabularyterm v on v.id = j.predictionmethods{kind}_id
        group by 1
    )"""


SQL = f"""
    with genes as ({_GENES}),
    {_methods_cte("matched", "matched")},
    {_methods_cte("not_matched", "notmatched")},
    {_methods_cte("not_called", "notcalled")}
    select s.primaryexternalid, s.symbol, s.curie,
           b.primaryexternalid, b.symbol, b.curie,
           p.identity, p.similarity, p.length, p.rank,
           m.names, nm.names, nc.names
    from genetogeneparalogy p
    join genes s on s.id = p.subjectgene_id
    join genes b on b.id = p.objectgene_id
    left join matched m on m.id = p.id
    left join not_matched nm on nm.id = p.id
    left join not_called nc on nc.id = p.id
    where not coalesce(p.obsolete, false)
"""


def to_row(values: list[str]) -> dict:
    (subject, subject_symbol, subject_taxon, paralog, paralog_symbol, paralog_taxon,
     identity, similarity, length, rank, matched, not_matched, not_called) = values
    return {
        "subjectGeneId": subject,
        "subjectGeneSymbol": subject_symbol,
        "subjectGeneTaxon": subject_taxon,
        "paralogueGeneId": paralog,
        "paralogueGeneSymbol": paralog_symbol,
        "paralogueGeneTaxon": paralog_taxon,
        # The fetcher writes `value or ""`, so a 0 metric becomes empty.
        "identity": "" if identity == "0" else identity,
        "similarity": "" if similarity == "0" else similarity,
        "length": "" if length == "0" else length,
        "rank": "" if rank == "0" else rank,
        "predictionMethodsMatched": matched,
        "predictionMethodsNotMatched": not_matched,
        "predictionMethodsNotCalled": not_called,
    }


def seeded_rows(rows, seed_rank: dict[str, int]) -> list[dict]:
    """One row per unordered gene pair: the direction whose subject is the
    earliest seed gene, in seed order (the order fetch_paralogs.py writes)."""
    kept: dict[tuple[str, str], dict] = {}
    for row in rows:
        subject, paralog = row["subjectGeneId"], row["paralogueGeneId"]
        if not (subject and paralog) or subject not in seed_rank:
            continue
        pair = (min(subject, paralog), max(subject, paralog))
        held = kept.get(pair)
        if held is None or seed_rank[subject] < seed_rank[held["subjectGeneId"]]:
            kept[pair] = row
    return sorted(kept.values(), key=lambda r: (seed_rank[r["subjectGeneId"]], r["paralogueGeneId"]))


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance paralogs from the curation DB into TSV")
    add_args(parser)
    args = parser.parse_args()

    configure_logging(args.verbose)
    started = time.monotonic()

    release = get_current_release()
    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit)
    seed_rank = {gid: i for i, gid in enumerate(gene_ids)}
    rows = [to_row(values) for values in copy_rows(SQL)]
    log.info("paralogs: %d rows in the database", len(rows))

    covered = {r["subjectGeneTaxon"] for r in rows}
    gap = gap_seeds(gene_ids, mod_gene_taxa(selected_mods(None, args.mods)), covered)
    if gap and not args.no_api_fallback:
        log.info("paralogs: %d seeds of taxa with no database rows go to the API", len(gap))
        fetched = api_rows(gap, "paralogs", fetch_paralogs._row, args)
        rows.extend(fetched)
        log.info("paralogs: %d API rows", len(fetched))
    elif gap:
        log.warning("paralogs: %d seeds of taxa with no database rows skipped", len(gap))

    rows = seeded_rows(rows, seed_rank)
    with TsvWriter(
        Path(args.out_dir) / "paralogs.tsv",
        fetch_paralogs.COLUMNS,
        release=release,
        source="curation DB genetogeneparalogy",
    ) as writer:
        for row in rows:
            writer.write_row(**row)

    log.info("paralogs: rows=%d seeds=%d elapsed=%.0fs",
             len(rows), len(gene_ids), time.monotonic() - started)


if __name__ == "__main__":
    main()
