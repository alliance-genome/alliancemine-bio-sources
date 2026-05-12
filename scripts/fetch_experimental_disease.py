#!/usr/bin/env python3
"""
Fetch per-gene experimental-evidence disease annotations from the Alliance API.

Endpoint: /gene/{id}/diseases-by-experiment

Disambiguates experimental disease annotations (gene-direct curation by a MOD)
from orthology-inferred disease annotations (already covered by
fetch_disease_annotations.py / alliance-disease-detail). Emits partial
DiseaseAnnotation items keyed on (experimentalGene, ontologyTerm, evidenceCode);
the integration engine merges these into rows already created by
alliance-disease-detail when keys collide, and otherwise inserts fresh rows
with evidenceType="experimental".

Usage:
    python3 fetch_experimental_disease.py
    python3 fetch_experimental_disease.py --mods MGI,ZFIN
    python3 fetch_experimental_disease.py --ids MGI:88276
    python3 fetch_experimental_disease.py --limit 100

Output: data/experimental-disease.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterator

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    enumerate_mod_genes,
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    iter_fms_datatype,
    join_pipe,
    normalize_taxon,
    open_cache,
    paginate,
)

DEFAULT_MODS = ["MGI", "ZFIN", "FB", "WB", "RGD", "SGD", "HUMAN", "XBXL", "XBXT"]


# MOD-prefix → NCBI taxon. Used because DAF rows do not carry taxon — must
# infer from data provider (MOD-curated DAF is always single-organism).
MOD_TAXON = {
    "MGI": "NCBITaxon:10090",
    "RGD": "NCBITaxon:10116",
    "ZFIN": "NCBITaxon:7955",
    "FB": "NCBITaxon:7227",
    "WB": "NCBITaxon:6239",
    "SGD": "NCBITaxon:559292",
    "HUMAN": "NCBITaxon:9606",
    "XBXL": "NCBITaxon:8355",
    "XBXT": "NCBITaxon:8364",
}


def _fms_rows(mods: list[str]) -> Iterator[dict]:
    """Yield TSV-row dicts for DAF records where objectType=='gene' AND
    primaryGeneticEntityIDs is non-empty — that combination marks
    experimental evidence (a gene-level disease annotation backed by an
    allele or AGM, vs orthology-inferred annotations which lack the
    genetic entity reference).
    """
    for mod in mods:
        log.info("FMS DAF extraction (experimental gene rows): %s", mod)
        taxon = MOD_TAXON.get(mod, "")
        for _sub, rec in iter_fms_datatype("DAF", sub_filter=mod):
            rel = rec.get("objectRelation") or {}
            if rel.get("objectType") != "gene":
                continue
            entities = rec.get("primaryGeneticEntityIDs") or []
            if not entities:
                continue
            gene_id = rec.get("objectId", "")
            disease_id = rec.get("DOid", "")
            if not gene_id or not disease_id:
                continue
            evidence = rec.get("evidence") or {}
            eco_codes = join_pipe(evidence.get("evidenceCodes") or [])
            pub = evidence.get("publication") or {}
            pmid = pub.get("publicationId", "")
            yield {
                "geneId": gene_id,
                "geneTaxon": taxon,
                "diseaseId": disease_id,
                "diseaseName": "",
                "relationName": rel.get("associationType", ""),
                "evidenceCodes": eco_codes,
                "evidenceAbbrs": "",          # not encoded in DAF; abbr is API-only
                "evidencePmids": pmid,
                "evidenceCurie": "",
                "negated": "",
                "dataProvider": mod,
            }

log = logging.getLogger("alliance.fetch.experimental_disease")


COLUMNS = [
    "geneId",
    "geneTaxon",
    "diseaseId",
    "diseaseName",
    "relationName",
    "evidenceCodes",
    "evidenceAbbrs",
    "evidencePmids",
    "evidenceCurie",
    "negated",
    "dataProvider",
]


def _pmids_from_xrefs(primary: dict) -> str:
    xrefs = get_nested(primary, "evidenceItem", "crossReferences", default=[]) or []
    return join_pipe(
        x.get("referencedCurie", "")
        for x in xrefs
        if isinstance(x, dict) and str(x.get("referencedCurie", "")).startswith("PMID:")
    )


def _row(seed_id: str, result: dict) -> dict | None:
    primary = (result.get("primaryAnnotations") or [{}])[0]
    obj = result.get("object") or {}
    disease_id = get_nested(obj, "curie", default="")
    if not disease_id:
        return None
    subj = primary.get("diseaseAnnotationSubject") or result.get("subject") or {}
    gene_taxon = normalize_taxon(get_nested(subj, "taxon", "curie", default=""))
    eco_codes = join_pipe(
        get_nested(c, "curie", default="")
        for c in (result.get("evidenceCodes") or [])
        if get_nested(c, "curie")
    )
    eco_abbrs = join_pipe(
        get_nested(c, "abbreviation", default="")
        for c in (result.get("evidenceCodes") or [])
        if get_nested(c, "abbreviation")
    )
    pmids = _pmids_from_xrefs(primary) or join_pipe(
        get_nested(r, "referenceID", default="")
        for r in (result.get("references") or [])
        if get_nested(r, "referenceID")
    )
    return {
        "geneId": seed_id,
        "geneTaxon": gene_taxon,
        "diseaseId": disease_id,
        "diseaseName": get_nested(obj, "name", default=""),
        "relationName": get_nested(result, "relation", "name", default=""),
        "evidenceCodes": eco_codes,
        "evidenceAbbrs": eco_abbrs,
        "evidencePmids": pmids,
        "evidenceCurie": get_nested(primary, "evidenceItem", "curie", default=""),
        "negated": str(primary.get("negated", "")).lower() if primary.get("negated") is not None else "",
        "dataProvider": get_nested(primary, "dataProvider", "abbreviation", default=""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Alliance experimental disease annotations into TSV")
    parser.add_argument(
        "--source",
        choices=["fms", "api"],
        default="fms",
        help="Backend: 'fms' (bulk DAF JSON.gz, filter to rows with primaryGeneticEntityIDs; "
             "default — the API endpoint /gene/{id}/diseases-by-experiment returns total:0 "
             "across all probed MODs as of 2026-05-11) or 'api' (per-gene fan-out).",
    )
    parser.add_argument("--ids", help="Comma-separated gene CURIEs (api path)")
    parser.add_argument(
        "--mods",
        default=",".join(DEFAULT_MODS),
        help="Comma-separated MOD prefixes",
    )
    parser.add_argument("--limit", type=int, help="Process at most N genes")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="ThreadPool size (default 4; keep low to avoid AGR WAF rate limit)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    out_dir = Path(args.out_dir)
    stats = FetchStats()
    seen_triples: set = set()

    if args.source == "fms":
        mods = [m.strip() for m in args.mods.split(",") if m.strip()]
        total = 0
        with TsvWriter(
            out_dir / "experimental-disease.tsv",
            COLUMNS,
            release=release,
            source=f"FMS DAF (gene rows w/ primaryGeneticEntityIDs; {','.join(mods)})",
        ) as writer:
            for row in _fms_rows(mods):
                triple = (row["geneId"], row["diseaseId"], row["evidenceCodes"])
                if triple in seen_triples:
                    continue
                seen_triples.add(triple)
                writer.write_row(**row)
                total += 1
                if args.limit is not None and total >= args.limit:
                    break
        log.info("FMS Done. rows=%d", total)
        return

    if args.ids:
        gene_ids = [g.strip() for g in args.ids.split(",") if g.strip()]
    else:
        mods = [m.strip() for m in args.mods.split(",") if m.strip()]
        if mods == ["SGD"]:
            gene_ids = enumerate_yeast_genes()
        else:
            gene_ids = enumerate_mod_genes(mods)
    if args.limit:
        gene_ids = gene_ids[: args.limit]

    log.info("Processing %d genes", len(gene_ids))

    def fetch_one(gid: str) -> tuple[str, list[dict], Exception | None]:
        rows: list[dict] = []
        try:
            with open_cache("experimental_disease") as cache:
                for result in paginate(
                    f"/gene/{gid}/diseases-by-experiment", cache=cache, stats=stats
                ):
                    row = _row(gid, result)
                    if row is not None:
                        rows.append(row)
            return gid, rows, None
        except Exception as e:
            return gid, [], e

    with TsvWriter(
        out_dir / "experimental-disease.tsv",
        COLUMNS,
        release=release,
        source="/gene/*/diseases-by-experiment",
    ) as writer:
        total_rows = 0
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(fetch_one, gid): gid for gid in gene_ids}
            try:
                from tqdm import tqdm
                pbar = tqdm(total=len(futures), desc="genes", unit="gene")
            except ImportError:
                pbar = None
            for fut in as_completed(futures):
                gid, rows, err = fut.result()
                if pbar:
                    pbar.update(1)
                if err is not None:
                    log.warning("Experimental-disease fetch failed for %s: %s", gid, err)
                    continue
                for row in rows:
                    triple = (row["geneId"], row["diseaseId"], row["evidenceCodes"])
                    if triple in seen_triples:
                        continue
                    seen_triples.add(triple)
                    writer.write_row(**row)
                    total_rows += 1
            if pbar:
                pbar.close()

    log.info("Done. rows=%d  API stats: %s", total_rows, stats)


if __name__ == "__main__":
    main()
