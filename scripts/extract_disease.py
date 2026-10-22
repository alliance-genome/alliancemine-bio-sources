#!/usr/bin/env python3
"""
Build disease-alliance.tsv (fetch_disease.py's DISEASE-ALLIANCE COMBINED
TSV) from a restored curation database, plus FMS for the orthology rows.

MOD-curated rows come from SQL, one per (annotation, row subject, ECO code)
as the Alliance indexer writes them: the annotation's own subject (gene,
allele or affected_genomic_model), then a gene row for each inferred or
asserted gene and an allele row for each inferred or asserted allele of an
AGM or allele annotation, with InferredFromID/Symbol naming the annotated
subject. Reference is the PMID of the annotation's reference (else its MOD
curie, else its AGRKB curie). Source is the data provider,
"<provider> Via <secondary provider>" when there is one.

Rows inferred via orthology (Source "Alliance", ECO:0000501, 78% of the
file) are computed by the Alliance indexer and are not in the curation
database. With --orthology fms (default) they are copied from the release's
FMS DISEASE-ALLIANCE COMBINED file, which is what fetch_disease.py serves;
--orthology none writes the curated rows only.

Output: data/disease-alliance.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    MOD_TAXON,
    TsvWriter,
    configure_logging,
    fms_download_gz,
    fms_list_files,
    get_current_release,
    join_pipe,
    parse_mods,
)
from disease_db import (  # noqa: E402
    MODIFIER_CONDITIONS,
    SUBJECT,
    Annotation,
    RefMap,
    annotations,
    condition_label,
    eco_terms,
    gene_relation,
    negated_relation,
    species_names,
)

log = logging.getLogger("alliance.extract.disease")

COLUMNS = [
    "Taxon", "SpeciesName", "DBobjectType", "DBObjectID", "DBObjectSymbol", "AssociationType",
    "DOID", "DOtermName", "WithOrtholog", "InferredFromID", "InferredFromSymbol",
    "ExperimentalCondition", "Modifier", "EvidenceCode", "EvidenceCodeName", "Reference", "Date", "Source",
]
OBJECT_TYPE = {"gene": "gene", "allele": "allele", "agm": "affected_genomic_model"}
ORTHOLOGY_SOURCE = "Alliance"


def source_name(ann: Annotation) -> str:
    if ann.secondary_provider:
        return f"{ann.provider} Via {ann.secondary_provider}"
    return ann.provider


def condition_columns(conditions: list[tuple[str, str]]) -> tuple[str, str]:
    """(ExperimentalCondition, Modifier): 'Has Condition: <summary>' entries,
    with ameliorated_by/exacerbated_by conditions in Modifier."""
    exp, mod = [], []
    for rel_type, summary in conditions:
        text = f"{condition_label(rel_type)}: {summary}"
        (mod if rel_type in MODIFIER_CONDITIONS else exp).append(text)
    return join_pipe(exp), join_pipe(mod)


def subject_relation(ann: Annotation) -> str:
    return negated_relation(ann.relation) if ann.negated else ann.relation


def rows_for(ann: Annotation, refs: RefMap, eco: dict[str, tuple[str, str]],
             species: dict[str, str]) -> Iterator[dict]:
    exp, mod = condition_columns(ann.conditions)
    base = {
        "AssociationType": subject_relation(ann),
        "DOID": ann.disease_id,
        "DOtermName": ann.disease_name,
        "WithOrtholog": join_pipe(ann.with_genes),
        "ExperimentalCondition": exp,
        "Modifier": mod,
        "Reference": refs.preferred(ann.reference),
        "Date": ann.date,
        "Source": source_name(ann),
    }
    entities = [(SUBJECT, ann.subject)]
    if ann.kind != "gene":
        entities += ann.genes()
    if ann.kind == "agm":
        entities += ann.alleles()
    for link, ent in entities:
        row = dict(base, Taxon=ent.taxon, SpeciesName=species.get(ent.taxon, ""),
                   DBobjectType=OBJECT_TYPE[ent.kind], DBObjectID=ent.id, DBObjectSymbol=ent.symbol)
        if link != SUBJECT:
            row["InferredFromID"] = ann.subject.id
            row["InferredFromSymbol"] = ann.subject.symbol
            row["ExperimentalCondition"] = row["Modifier"] = ""
            row["AssociationType"] = gene_relation(ann.relation, ann.negated)
        for code in ann.evidence_codes or [""]:
            yield dict(row, EvidenceCode=code, EvidenceCodeName=eco.get(code, ("", ""))[0])


def fms_orthology_lines(release: str) -> Iterator[list[str]]:
    """Source=Alliance rows of the release's FMS DISEASE-ALLIANCE COMBINED file."""
    files = fms_list_files("DISEASE-ALLIANCE", release=release, sub_filter="COMBINED")
    if not files:
        raise SystemExit("No FMS DISEASE-ALLIANCE COMBINED file for the orthology rows")
    for fi in files:
        with gzip.open(fms_download_gz(fi["url"]), "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) == len(COLUMNS) and cols[-1] == ORTHOLOGY_SOURCE:
                    yield cols


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract the Alliance DISEASE-ALLIANCE TSV from the curation DB")
    parser.add_argument("--mods", help="Comma-separated MODs; keep rows of their taxa (default all)")
    parser.add_argument("--limit", type=int, help="Write at most N curated rows")
    parser.add_argument("--orthology", choices=["fms", "none"], default="fms",
                        help="fms: add the via-orthology rows of the FMS COMBINED file (default); "
                             "none: curated rows only")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    configure_logging(args.verbose)

    release = get_current_release()
    mods = parse_mods(args.mods)
    taxa = {MOD_TAXON[m] for m in mods} if mods else None
    refs = RefMap.load()
    eco = eco_terms()
    species = species_names()

    curated = 0
    orthology = 0
    source = "curation DB" + (" + FMS DISEASE-ALLIANCE COMBINED (Source=Alliance)" if args.orthology == "fms" else "")
    with TsvWriter(Path(args.out_dir) / "disease-alliance.tsv", COLUMNS, release=release, source=source) as writer:
        for ann in sorted(annotations(), key=lambda a: (a.subject.id, a.disease_id, a.curie)):
            for row in rows_for(ann, refs, eco, species):
                if taxa is not None and row["Taxon"] not in taxa:
                    continue
                writer.write_row(**row)
                curated += 1
                if args.limit and curated >= args.limit:
                    break
            if args.limit and curated >= args.limit:
                break
        refs.log_misses("disease")
        if args.orthology == "fms":
            for cols in fms_orthology_lines(release):
                if taxa is not None and cols[0] not in taxa:
                    continue
                writer.write_row(**dict(zip(COLUMNS, cols)))
                orthology += 1
    log.info("disease: rows=%d (curated %d, orthology %d)", curated + orthology, curated, orthology)


if __name__ == "__main__":
    main()
