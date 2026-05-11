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

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    FetchStats,
    TsvWriter,
    configure_logging,
    enumerate_mod_genes,
    enumerate_yeast_genes,
    get_current_release,
    get_nested,
    join_pipe,
    normalize_taxon,
    open_cache,
    paginate,
)

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
    parser.add_argument("--ids", help="Comma-separated gene CURIEs")
    parser.add_argument(
        "--mods",
        default="MGI,ZFIN,FB,WB,RGD,SGD,XBXL,XBXT,HUMAN",
        help="Comma-separated MOD prefixes for seed enumeration",
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
    release = get_current_release()
    out_dir = Path(args.out_dir)

    stats = FetchStats()
    seen_triples: set = set()

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
