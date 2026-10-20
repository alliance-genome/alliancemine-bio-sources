#!/usr/bin/env python3
"""
Fetch molecular + genetic interactions from the Alliance API, emit two TSVs.

Each gene gets two API calls (paginated) — one for each endpoint. Rows are
deduplicated by the API's own `uniqueId` field so the converter doesn't need
to dedupe again. A gene whose either endpoint fails is recorded as failed and
contributes no rows; a rerun picks it up from the cache.

Usage:
    python3 fetch_interactions.py                  # every gene of the 9 MODs
    python3 fetch_interactions.py --mods SGD,HUMAN --workers 6
    python3 fetch_interactions.py --ids SGD:S000004103,SGD:S000002429
    python3 fetch_interactions.py --limit 50       # smoke test

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
    add_parallel_args,
    configure_logging,
    finish_run,
    get_current_release,
    get_nested,
    join_pipe,
    paginate,
    resolve_gene_ids,
    run_parallel,
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


def fetch_for_gene(gene_id: str, cache, stats: FetchStats) -> tuple[list[dict], list[dict]]:
    """Fetch both endpoints for one gene; return (molecular rows, genetic rows)."""
    mol = [
        _row_for_interaction(gene_id, "molecular", r)
        for r in paginate(f"/gene/{gene_id}/molecular-interactions", cache=cache, stats=stats)
    ]
    gen = [
        _row_for_interaction(gene_id, "genetic", r)
        for r in paginate(f"/gene/{gene_id}/genetic-interactions", cache=cache, stats=stats)
    ]
    return mol, gen


def _is_canonical(row: dict) -> bool:
    return row.get("geneId", "") <= row.get("partnerGeneId", "")


def _write_new(writer: TsvWriter, rows: list[dict], seen: set, pending: dict) -> int:
    """Write each uniqueId once, in the orientation geneId <= partnerGeneId.

    The API returns one interaction from both genes' endpoints, with the
    queried gene first. A non-canonical copy is held in `pending` until
    `_flush_pending`, so completion order never decides the kept direction.
    """
    written = 0
    for row in rows:
        uid = row.get("uniqueId")
        if not uid:
            writer.write_row(**row)
            written += 1
        elif uid in seen:
            continue
        elif _is_canonical(row):
            seen.add(uid)
            pending.pop(uid, None)
            writer.write_row(**row)
            written += 1
        else:
            pending.setdefault(uid, row)
    return written


def _flush_pending(writer: TsvWriter, pending: dict, seen: set) -> int:
    """Write the non-canonical copies whose canonical twin never arrived, sorted by uniqueId."""
    written = 0
    for uid in sorted(pending):
        if uid not in seen:
            seen.add(uid)
            writer.write_row(**pending[uid])
            written += 1
    pending.clear()
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance interactions into TSVs")
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (skips enumeration)")
    parser.add_argument("--mods", help="Comma-separated MODs to seed from (default: all 9, env ALLIANCE_MODS)")
    parser.add_argument("--limit", type=int, help="Process at most N genes (for smoke tests)")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSVs (default: ../data)",
    )
    add_parallel_args(parser)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    gene_ids = resolve_gene_ids(args.ids, args.mods, args.limit)
    out_dir = Path(args.out_dir)
    stats = FetchStats()
    seen_mol: set = set()
    seen_gen: set = set()
    pending_mol: dict = {}
    pending_gen: dict = {}
    totals = {"molecular": 0, "genetic": 0}

    with TsvWriter(
        out_dir / "molecular-interactions.tsv",
        MOLECULAR_COLUMNS,
        release=release,
        source="/gene/*/molecular-interactions",
    ) as mol_writer, TsvWriter(
        out_dir / "genetic-interactions.tsv",
        GENETIC_COLUMNS,
        release=release,
        source="/gene/*/genetic-interactions",
    ) as gen_writer:
        def consume(_gid, result):
            mol, gen = result
            totals["molecular"] += _write_new(mol_writer, mol, seen_mol, pending_mol)
            totals["genetic"] += _write_new(gen_writer, gen, seen_gen, pending_gen)

        report = run_parallel(gene_ids, fetch_for_gene, consume, cache_name="interactions",
                              stats=stats, workers=args.workers, label="genes")
        totals["molecular"] += _flush_pending(mol_writer, pending_mol, seen_mol)
        totals["genetic"] += _flush_pending(gen_writer, pending_gen, seen_gen)
        finish_run(report, stats, out_dir, args.max_failure_rate,
                   rows_written=mol_writer.rows_written + gen_writer.rows_written)

    log.info("interactions: molecular rows=%d  genetic rows=%d", totals["molecular"], totals["genetic"])


if __name__ == "__main__":
    main()
