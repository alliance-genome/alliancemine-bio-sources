#!/usr/bin/env python3
"""
Build experimental-disease.tsv (fetch_experimental_disease.py's TSV) from a
restored curation database.

The fetcher keeps the DAF gene rows that name primary genetic entities: gene
annotations a MOD backs with an allele or genotype. In the curation database
these are the inferred and asserted genes of AGM and allele annotations. Rows
group on (gene, disease, relation, ECO codes), the fetcher's key plus the
relation; references are the union over the group, evidenceCurie the first
annotation's AGRKB reference, evidenceAbbrs the ECO abbreviations (IMP, TAS).

evidencePmids lists the PubMed and MOD curies of the references
(reference_crossreference).

Output: data/experimental-disease.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    MOD_TAXON,
    TsvWriter,
    configure_logging,
    get_current_release,
    join_pipe,
    normalize_taxon,
    selected_mods,
)
from disease_db import SUBJECT, RefMap, annotations, eco_terms, gene_relation  # noqa: E402
from fetch_experimental_disease import COLUMNS, DEFAULT_MODS  # noqa: E402

log = logging.getLogger("alliance.extract.experimental_disease")


def experimental_rows(anns, refs: RefMap, eco: dict[str, tuple[str, str]],
                      taxa: set[str] | None = None, gene_ids: set[str] | None = None) -> list[dict]:
    groups: dict[tuple, dict] = {}
    for ann in anns:
        relation = gene_relation(ann.relation, ann.negated)
        codes = join_pipe(sorted(set(ann.evidence_codes)))
        for link, gene in ann.genes():
            if link == SUBJECT:
                continue
            if taxa is not None and gene.taxon not in taxa:
                continue
            if gene_ids is not None and gene.id not in gene_ids:
                continue
            key = (gene.id, ann.disease_id, relation, codes)
            g = groups.get(key)
            if g is None:
                g = groups[key] = {"ann": ann, "gene": gene, "refs": []}
            if ann.reference and ann.reference not in g["refs"]:
                g["refs"].append(ann.reference)
    rows = []
    for (gene_id, doid, relation, codes), g in groups.items():
        ann, gene = g["ann"], g["gene"]
        pubs: list[str] = []
        for ref in g["refs"]:
            pubs.extend(c for c in refs.curies(ref) if c not in pubs)
        rows.append({
            "geneId": gene_id,
            "geneTaxon": normalize_taxon(gene.taxon),
            "diseaseId": doid,
            "diseaseName": ann.disease_name,
            "relationName": relation,
            "evidenceCodes": codes,
            "evidenceAbbrs": join_pipe(eco.get(c, ("", ""))[1] for c in codes.split("|") if c),
            "evidencePmids": join_pipe(pubs),
            "evidenceCurie": g["refs"][0] if g["refs"] else "",
            "negated": "true" if ann.negated else "false",
            "dataProvider": ann.provider,
        })
    rows.sort(key=lambda r: (r["geneId"], r["diseaseId"], r["relationName"], r["evidenceCodes"]))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance experimental disease annotations from the curation DB")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--limit", type=int, help="Write at most N rows")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    release = get_current_release()
    mods = selected_mods(DEFAULT_MODS, args.mods)
    taxa = {MOD_TAXON[m] for m in mods}
    gene_ids = {g.strip() for g in args.ids.split(",") if g.strip()} if args.ids else None
    refs = RefMap.load()

    rows = experimental_rows(annotations(), refs, eco_terms(), taxa, gene_ids)
    refs.log_misses("experimental_disease")
    if args.limit:
        rows = rows[: args.limit]
    with TsvWriter(Path(args.out_dir) / "experimental-disease.tsv", COLUMNS, release=release,
                   source=f"curation DB (genes of allele/AGM annotations; {','.join(mods)})") as writer:
        for row in rows:
            writer.write_row(**row)
    log.info("experimental_disease: rows=%d", len(rows))


if __name__ == "__main__":
    main()
