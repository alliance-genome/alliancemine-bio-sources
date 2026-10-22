#!/usr/bin/env python3
"""
Build disease-annotations-detail.tsv (fetch_disease_annotations.py's TSV)
from a restored curation database, plus the API for the orthology rows.

Direct rows (viaOrthologyOrder 0) come from SQL. One row per (gene, disease,
relation, qualifiers, with genes), as /disease/{id}/genes groups them: gene annotations
give their subject; AGM and allele annotations give their inferred and
asserted genes, with is_model_of (and its ameliorated/exacerbated forms)
reported as is_implicated_in and negated relations as is_not_*. ECO codes and
references are the union over the grouped annotations; parentSlimIds is the
disease and its DO ancestors; uniqueId is the lowest AGRKB annotation curie.

Rows inferred via orthology (viaOrthologyOrder 1, ~87% of the API's rows) are
computed by the Alliance indexer and are not in the curation database. With
--orthology api (default) they come from /disease/{id}/genes exactly as
fetch_disease_annotations.py fetches them, for every DOID the database
annotates; --orthology none writes the direct rows only.

pubmedIds lists the PubMed and MOD curies of the references
(reference_crossreference); referenceIds their AGRKB curies.

Output: data/disease-annotations-detail.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    MOD_TAXON,
    FetchStats,
    TsvWriter,
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    join_pipe,
    paginate,
    parse_mods,
    run_parallel,
)
from disease_db import RefMap, annotations, do_ancestors, gene_relation  # noqa: E402
from fetch_disease_annotations import COLUMNS, _row as api_row  # noqa: E402

log = logging.getLogger("alliance.extract.disease_annotations")


def group_key(gene_id: str, ann) -> tuple:
    """The API keeps annotations made with orthologous genes (SGD's
    ECO:0000250 with HGNC genes) apart from those without."""
    return (gene_id, ann.disease_id, gene_relation(ann.relation, ann.negated), "|".join(sorted(ann.qualifiers)),
            "|".join(sorted(ann.with_genes)))


def direct_rows(anns, ancestors: dict[str, list[str]], refs: RefMap,
                taxa: set[str] | None = None) -> list[dict]:
    """Gene-level rows of the curation annotations, grouped as the API groups them."""
    groups: dict[tuple, dict] = {}
    for ann in anns:
        for _link, gene in ann.genes():
            if taxa is not None and gene.taxon not in taxa:
                continue
            key = group_key(gene.id, ann)
            g = groups.get(key)
            if g is None:
                g = groups[key] = {"ann": ann, "gene": gene, "eco": set(), "refs": [], "curies": []}
            g["eco"].update(ann.evidence_codes)
            if ann.reference and ann.reference not in g["refs"]:
                g["refs"].append(ann.reference)
            g["curies"].append(ann.curie)
    rows = []
    for (gene_id, doid, relation, quals, _with), g in groups.items():
        ann, gene = g["ann"], g["gene"]
        pubs: list[str] = []
        for ref in g["refs"]:
            pubs.extend(c for c in refs.curies(ref) if c not in pubs)
        rows.append({
            "diseaseId": doid,
            "diseaseName": ann.disease_name,
            "subjectGeneId": gene_id,
            "subjectGeneTaxon": gene.taxon,
            "generatedRelationString": relation,
            "diseaseQualifiers": quals,
            "evidenceCodes": join_pipe(sorted(g["eco"])),
            "parentSlimIds": join_pipe([doid] + ancestors.get(doid, [])),
            "viaOrthologyOrder": "0",
            "pubmedIds": join_pipe(pubs),
            "referenceIds": join_pipe(g["refs"]),
            "uniqueId": min(g["curies"]),
        })
    rows.sort(key=lambda r: (r["diseaseId"], r["subjectGeneId"], r["generatedRelationString"],
                             r["diseaseQualifiers"]))
    return rows


def is_orthology_row(row: dict) -> bool:
    return row["viaOrthologyOrder"] not in ("", "0")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance disease annotation detail from the curation DB")
    parser.add_argument("--ids", help="Comma-separated DOID CURIEs")
    parser.add_argument("--limit", type=int, help="Process at most N diseases")
    parser.add_argument("--mods", help="Comma-separated MODs; keep genes of their taxa (default all)")
    parser.add_argument("--orthology", choices=["api", "none"], default="api",
                        help="api: add the via-orthology rows from /disease/{id}/genes (default); "
                             "none: direct rows from the database only")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    add_parallel_args(parser)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    release = get_current_release()
    mods = parse_mods(args.mods)
    taxa = {MOD_TAXON[m] for m in mods} if mods else None
    refs = RefMap.load()

    anns = annotations()
    disease_ids = sorted({a.disease_id for a in anns})
    if args.ids:
        disease_ids = [d.strip() for d in args.ids.split(",") if d.strip()]
    if args.limit:
        disease_ids = disease_ids[: args.limit]
    keep = set(disease_ids)
    rows = direct_rows((a for a in anns if a.disease_id in keep), do_ancestors(), refs, taxa)
    refs.log_misses("disease_annotations")
    log.info("disease_annotations: %d direct rows over %d diseases", len(rows), len(keep))

    out_dir = Path(args.out_dir)
    source = "curation DB" + (" + /disease/*/genes (via orthology)" if args.orthology == "api" else "")
    with TsvWriter(out_dir / "disease-annotations-detail.tsv", COLUMNS, release=release, source=source) as writer:
        for row in rows:
            writer.write_row(**row)
        if args.orthology == "none":
            return
        stats = FetchStats()
        seen_uids: set = set()

        def work(did, cache, stats):
            return [api_row(did, r) for r in paginate(f"/disease/{did}/genes", cache=cache, stats=stats)]

        def consume(_did, api_rows):
            for row in api_rows:
                if not is_orthology_row(row):
                    continue
                if taxa is not None and row["subjectGeneTaxon"] not in taxa:
                    continue
                uid = row["uniqueId"]
                if uid and uid in seen_uids:
                    continue
                if uid:
                    seen_uids.add(uid)
                writer.write_row(**row)

        report = run_parallel(disease_ids, work, consume, cache_name="disease_annotations", stats=stats,
                              workers=args.workers, label="diseases")
        finish_run(report, stats, out_dir, args.max_failure_rate, rows_written=writer.rows_written)
    log.info("disease_annotations: rows=%d (direct %d)", writer.rows_written, len(rows))


if __name__ == "__main__":
    main()
