#!/usr/bin/env python3
"""
Fetch molecular + genetic interactions from the Alliance API, emit two TSVs.

Each gene gets two API calls (paginated) — one for each endpoint. Rows are
written as interactions are discovered, deduplicated by the API's own
`uniqueId` field so the converter doesn't need to dedupe again.

Usage:
    python3 fetch_interactions.py
    python3 fetch_interactions.py --ids SGD:S000004103,SGD:S000002429
    python3 fetch_interactions.py --limit 50   # yeast genome head-50 smoke test

Output: data/molecular-interactions.tsv, data/genetic-interactions.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Make sibling common.py importable whether invoked as module or script
sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    join_pipe,
    open_cache,
    paginate,
)

log = logging.getLogger("alliance.fetch.interactions")


# Shared column schema for both interaction kinds. Any field missing from a
# particular endpoint is written as the empty string.
COMMON_COLUMNS = [
    "geneId",
    "geneSymbol",
    "geneTaxon",
    "partnerGeneId",
    "partnerSymbol",
    "partnerTaxon",
    "interactionType",
    "interactionTypeName",
    "interactorAType",
    "interactorATypeName",
    "interactorBType",
    "interactorBTypeName",
    "interactorARole",
    "interactorARoleName",
    "interactorBRole",
    "interactorBRoleName",
    "interactionSource",
    "interactionSourceName",
    "relation",
    "interactionId",
    "uniqueId",
    "pubmedId",
    "referenceId",
    "shortCitation",
    "crossRefs",
]
MOLECULAR_COLUMNS = COMMON_COLUMNS + [
    "detectionMethod",
    "detectionMethodName",
    "aggregationDatabase",
    "aggregationDatabaseName",
]
GENETIC_COLUMNS = COMMON_COLUMNS + [
    "phenotypesOrTraits",
]


def _extract_gene(node: dict, role: str = "subject") -> dict:
    """Pull partner-gene fields from one of the associationSubject/Object nodes."""
    return {
        f"{role}Id": get_nested(node, "primaryExternalId", default=""),
        f"{role}Symbol": get_nested(node, "geneSymbol", "displayText", default=""),
        f"{role}Taxon": get_nested(node, "taxon", "curie", default=""),
    }


def _extract_evidence(evidence_list: list) -> dict:
    """Take the first Reference in the evidence list (typical case)."""
    if not evidence_list:
        return {"pubmedId": "", "referenceId": "", "shortCitation": ""}
    first = evidence_list[0]
    return {
        "pubmedId": first.get("referenceID", ""),
        "referenceId": first.get("curie", ""),
        "shortCitation": first.get("shortCitation", ""),
    }


def _extract_term(node: dict | None) -> tuple[str, str]:
    """Return (curie, name) for an ontology-term-shaped API node."""
    if not node:
        return "", ""
    return node.get("curie", ""), node.get("name", "")


def _extract_crossrefs(node: dict) -> str:
    """Flatten the crossReferences list to a pipe-separated list of curies."""
    xrefs = node.get("crossReferences") or []
    return join_pipe(x.get("referencedCurie", "") for x in xrefs)


def _row_for_interaction(gene_id: str, kind: str, wrapper: dict) -> dict:
    """Convert one API result row into a flat dict keyed by our TSV columns.

    `kind` is 'molecular' or 'genetic'. The API wraps the payload in
    `geneMolecularInteraction` or `geneGeneticInteraction`.
    """
    inner_key = "geneMolecularInteraction" if kind == "molecular" else "geneGeneticInteraction"
    inner = wrapper.get(inner_key) or {}

    subject = inner.get("geneAssociationSubject") or {}
    object_ = inner.get("geneGeneAssociationObject") or {}

    # If the gene we queried is actually on the "object" side, swap so
    # geneId consistently refers to our subject gene.
    if subject.get("primaryExternalId") != gene_id and object_.get("primaryExternalId") == gene_id:
        subject, object_ = object_, subject

    subj_fields = _extract_gene(subject, "gene")
    obj_fields = _extract_gene(object_, "partnerGene")

    it_curie, it_name = _extract_term(inner.get("interactionType"))
    ia_t_curie, ia_t_name = _extract_term(inner.get("interactorAType"))
    ib_t_curie, ib_t_name = _extract_term(inner.get("interactorBType"))
    ia_r_curie, ia_r_name = _extract_term(inner.get("interactorARole"))
    ib_r_curie, ib_r_name = _extract_term(inner.get("interactorBRole"))
    is_curie, is_name = _extract_term(inner.get("interactionSource"))
    ev = _extract_evidence(inner.get("evidence") or [])

    row: dict = {
        "geneId": subj_fields["geneId"],
        "geneSymbol": subj_fields["geneSymbol"],
        "geneTaxon": subj_fields["geneTaxon"],
        "partnerGeneId": obj_fields["partnerGeneId"],
        "partnerSymbol": obj_fields["partnerGeneSymbol"],
        "partnerTaxon": obj_fields["partnerGeneTaxon"],
        "interactionType": it_curie,
        "interactionTypeName": it_name,
        "interactorAType": ia_t_curie,
        "interactorATypeName": ia_t_name,
        "interactorBType": ib_t_curie,
        "interactorBTypeName": ib_t_name,
        "interactorARole": ia_r_curie,
        "interactorARoleName": ia_r_name,
        "interactorBRole": ib_r_curie,
        "interactorBRoleName": ib_r_name,
        "interactionSource": is_curie,
        "interactionSourceName": is_name,
        "relation": get_nested(inner, "relation", "name", default=""),
        "interactionId": inner.get("interactionId", ""),
        "uniqueId": inner.get("uniqueId", ""),
        **ev,
        "crossRefs": _extract_crossrefs(inner),
    }

    if kind == "molecular":
        dm_curie, dm_name = _extract_term(inner.get("detectionMethod"))
        ag_curie, ag_name = _extract_term(inner.get("aggregationDatabase"))
        row["detectionMethod"] = dm_curie
        row["detectionMethodName"] = dm_name
        row["aggregationDatabase"] = ag_curie
        row["aggregationDatabaseName"] = ag_name
    else:
        phenotypes = inner.get("phenotypesOrTraits") or []
        row["phenotypesOrTraits"] = join_pipe(phenotypes)

    return row


def fetch_for_gene(
    gene_id: str,
    mol_writer: TsvWriter,
    gen_writer: TsvWriter,
    seen_mol: set,
    seen_gen: set,
    *,
    cache,
    stats: FetchStats,
) -> tuple[int, int]:
    """Fetch both endpoints for one gene; return (mol_count, gen_count)."""
    mol_count = 0
    gen_count = 0
    try:
        for result in paginate(f"/gene/{gene_id}/molecular-interactions", cache=cache, stats=stats):
            row = _row_for_interaction(gene_id, "molecular", result)
            uid = row.get("uniqueId")
            if uid and uid in seen_mol:
                continue
            if uid:
                seen_mol.add(uid)
            mol_writer.write_row(**row)
            mol_count += 1
    except Exception as e:
        log.warning("Molecular fetch failed for %s: %s", gene_id, e)

    try:
        for result in paginate(f"/gene/{gene_id}/genetic-interactions", cache=cache, stats=stats):
            row = _row_for_interaction(gene_id, "genetic", result)
            uid = row.get("uniqueId")
            if uid and uid in seen_gen:
                continue
            if uid:
                seen_gen.add(uid)
            gen_writer.write_row(**row)
            gen_count += 1
    except Exception as e:
        log.warning("Genetic fetch failed for %s: %s", gene_id, e)

    return mol_count, gen_count


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance interactions into TSVs")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--limit", type=int, help="Process at most N genes (for smoke tests)")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSVs (default: ../data)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    if args.ids:
        gene_ids = [g.strip() for g in args.ids.split(",") if g.strip()]
    else:
        gene_ids = enumerate_yeast_genes()
    if args.limit:
        gene_ids = gene_ids[: args.limit]

    log.info("Processing %d genes", len(gene_ids))
    release = get_current_release()
    out_dir = Path(args.out_dir)

    # Progress bar if available, otherwise fall through
    try:
        from tqdm import tqdm
        iterator = tqdm(gene_ids, desc="genes", unit="gene")
    except ImportError:
        iterator = gene_ids

    stats = FetchStats()
    seen_mol: set = set()
    seen_gen: set = set()

    with open_cache("interactions") as cache, \
         TsvWriter(
             out_dir / "molecular-interactions.tsv",
             MOLECULAR_COLUMNS,
             release=release,
             source=f"{getattr(cache, 'api_base', '')} /gene/*/molecular-interactions",
         ) as mol_writer, \
         TsvWriter(
             out_dir / "genetic-interactions.tsv",
             GENETIC_COLUMNS,
             release=release,
             source="/gene/*/genetic-interactions",
         ) as gen_writer:
        total_mol = total_gen = 0
        for gid in iterator:
            mol, gen = fetch_for_gene(gid, mol_writer, gen_writer, seen_mol, seen_gen, cache=cache, stats=stats)
            total_mol += mol
            total_gen += gen

    log.info("Done. Molecular rows: %d, genetic rows: %d", total_mol, total_gen)
    log.info("API stats: %s", stats)


if __name__ == "__main__":
    main()
