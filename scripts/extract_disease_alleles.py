#!/usr/bin/env python3
"""
Build disease-alleles.tsv (fetch_disease_alleles.py's TSV) from a restored
curation database.

Rows are the allele-level disease annotations: an allele annotation's
subject, and the inferred and asserted alleles of an AGM annotation (MGI and
ZFIN infer one allele per genotype or fish; FB asserts the alleles of each
genotype). Rows group on (allele, disease, relation): ECO codes and
references are the union over the group, evidenceCurie is the first
annotation's AGRKB reference. relationName is is_implicated_in, or
is_not_implicated_in for a negated annotation.

The fetcher's FMS path reads the MODs' DAF submissions, which predate FB's
move to genotype (AGM) annotations; the database also has XB and the
negation flag the DAF rows drop.

evidencePmids lists the PubMed and MOD curies of the references
(reference_crossreference).

Output: data/disease-alleles.tsv
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
from disease_db import RefMap, annotations, gene_relation  # noqa: E402
from fetch_disease_alleles import COLUMNS, DEFAULT_MODS  # noqa: E402

log = logging.getLogger("alliance.extract.disease_alleles")


def allele_rows(anns, refs: RefMap, taxa: set[str] | None = None,
                disease_ids: set[str] | None = None) -> list[dict]:
    groups: dict[tuple, dict] = {}
    for ann in anns:
        if disease_ids is not None and ann.disease_id not in disease_ids:
            continue
        relation = gene_relation(ann.relation, ann.negated)
        for _link, allele in ann.alleles():
            if taxa is not None and allele.taxon not in taxa:
                continue
            key = (allele.id, ann.disease_id, relation)
            g = groups.get(key)
            if g is None:
                g = groups[key] = {"ann": ann, "allele": allele, "eco": set(), "refs": []}
            g["eco"].update(ann.evidence_codes)
            if ann.reference and ann.reference not in g["refs"]:
                g["refs"].append(ann.reference)
    rows = []
    for (allele_id, doid, relation), g in groups.items():
        ann, allele = g["ann"], g["allele"]
        pubs: list[str] = []
        for ref in g["refs"]:
            pubs.extend(c for c in refs.curies(ref) if c not in pubs)
        rows.append({
            "diseaseId": doid,
            "diseaseName": ann.disease_name,
            "alleleId": allele_id,
            "alleleSymbol": allele.symbol,
            "alleleTaxon": normalize_taxon(allele.taxon),
            "relationName": relation,
            "evidenceCodes": join_pipe(sorted(g["eco"])),
            "evidencePmids": join_pipe(pubs),
            "evidenceCurie": g["refs"][0] if g["refs"] else "",
            "negated": "true" if ann.negated else "false",
            "dataProvider": ann.provider,
        })
    rows.sort(key=lambda r: (r["alleleId"], r["diseaseId"], r["relationName"]))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance allele disease annotations from the curation DB")
    parser.add_argument(
        "--mods",
        help=f"Comma-separated MOD prefixes (default: {','.join(DEFAULT_MODS)}, env ALLIANCE_MODS)",
    )
    parser.add_argument("--ids", help="Comma-separated DOID CURIEs")
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
    disease_ids = {d.strip() for d in args.ids.split(",") if d.strip()} if args.ids else None
    refs = RefMap.load()

    rows = allele_rows(annotations(), refs, taxa, disease_ids)
    refs.log_misses("disease_alleles")
    if args.limit:
        rows = rows[: args.limit]
    with TsvWriter(Path(args.out_dir) / "disease-alleles.tsv", COLUMNS, release=release,
                   source=f"curation DB ({','.join(mods)})") as writer:
        for row in rows:
            writer.write_row(**row)
    log.info("disease_alleles: rows=%d", len(rows))


if __name__ == "__main__":
    main()
