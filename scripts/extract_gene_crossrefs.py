#!/usr/bin/env python3
"""
Build gene-crossrefs.tsv from a restored Alliance curation database dump.

Same file, columns and value formats as fetch_gene_crossrefs.py, which reads
the FMS GENECROSSREFERENCEJSON file (last published for 9.0.0). Rows come
from each gene's crossReferences plus its gcrpCrossReference (the canonical
protein, the only UniProtKB id SGD genes carry; fetch_alphafold.py needs it).

    xrefId    the crossreference's referencedCurie
    xrefUrl   the page urlTemplate with [%s] set to the curie's local part
    xrefType  the page name in FMS spelling: default -> generic_cross_reference,
              expression_atlas -> gene/expressionAtlas, biogrid/orcs ->
              gene/biogrid_orcs, PANTHER's per-MOD tree pages -> gene/panther,
              and "gene" for the gene's own id (HGNC files it under default)

PANTHER families on the default page get the tree page named after the gene's
prefix (MGI, HGNC, ...) when one exists, the URL the FMS file carried.

One row per (gene, xrefId), as in the fetcher. When a curie has several pages
the row keeps "gene", else "default", else the first gene/* page by name.
Genes are those of the per-MOD gene lists (resolve_gene_ids), so the converter
never creates a Gene that alliance-genes.tsv lacks.

Usage:
    python3 extract_gene_crossrefs.py
    python3 extract_gene_crossrefs.py --mods SGD --limit 1000

Environment: CURATION_PSQL or CURATION_DB_URL (see curation_db.py).
Output: data/gene-crossrefs.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    MOD_TAXON,
    TsvWriter,
    configure_logging,
    get_current_release,
    normalize_taxon,
    resolve_gene_ids,
    selected_mods,
)
from curation_db import copy_rows  # noqa: E402
from extract_genes import gene_base_cte  # noqa: E402
from fetch_gene_crossrefs import COLUMNS  # noqa: E402

log = logging.getLogger("alliance.extract.gene_crossrefs")

# Page names whose FMS spelling differs; PANTHER pages named after a MOD are mapped in xref_type().
PAGE_TYPE = {
    "default": "generic_cross_reference",
    "expression_atlas": "gene/expressionAtlas",
    "biogrid/orcs": "gene/biogrid_orcs",
}

# One row per gene and curie, page chosen by rank then name (see module docstring).
XREF_SQL = """
WITH {gene_base},
x AS (
  SELECT gc.genomicentity_id AS gene_id, c.referencedcurie, c.resourcedescriptorpage_id
  FROM genomicentity_crossreference gc
  JOIN g ON g.id = gc.genomicentity_id
  JOIN crossreference c ON c.id = gc.crossreferences_id
  WHERE NOT coalesce(c.obsolete, false) AND NOT coalesce(c.internal, false)
  UNION ALL
  SELECT g.id, c.referencedcurie, c.resourcedescriptorpage_id
  FROM g JOIN crossreference c ON c.id = g.gcrpcrossreference_id
)
SELECT DISTINCT ON (g.primaryexternalid, x.referencedcurie)
       g.primaryexternalid, g.taxon, x.referencedcurie, rp.name, rp.urltemplate, rd.prefix
FROM x
JOIN g ON g.id = x.gene_id
LEFT JOIN resourcedescriptorpage rp ON rp.id = x.resourcedescriptorpage_id
LEFT JOIN resourcedescriptor rd ON rd.id = rp.resourcedescriptor_id
WHERE coalesce(x.referencedcurie, '') <> ''
ORDER BY g.primaryexternalid, x.referencedcurie,
         CASE WHEN rp.name = 'gene' THEN 0 WHEN rp.name = 'default' THEN 1
              WHEN rp.name LIKE 'gene/%' THEN 2 ELSE 3 END,
         rp.name
"""


PANTHER_PAGES_SQL = """
SELECT p.name, p.urltemplate
FROM resourcedescriptorpage p JOIN resourcedescriptor r ON r.id = p.resourcedescriptor_id
WHERE r.prefix = 'PANTHER' AND p.name <> 'default'
"""


def xref_sql(taxa: list[str], ids: list[str] | None = None) -> str:
    return XREF_SQL.format(gene_base=gene_base_cte(taxa, ids))


def xref_type(page: str, prefix: str) -> str:
    if prefix == "PANTHER" and page not in ("", "default"):
        return "gene/panther"
    return PAGE_TYPE.get(page, page)


def xref_url(curie: str, template: str, prefix: str, gene_id: str) -> str:
    """Fill [%s] with the curie's local part. PANTHER tree pages take the
    family as PAN_BOOK and the gene's local id as [%s], as in the FMS file."""
    if not template:
        return ""
    local = curie.split(":", 1)[-1]
    if prefix == "PANTHER" and "PAN_BOOK" in template:
        return (template.replace("PAN_BOOK", local.split(":", 1)[0])
                .replace("[%s]", gene_id.split(":", 1)[-1]))
    return template.replace("[%s]", local)


def row_from_db(gene_id: str, taxon: str, curie: str, page: str, template: str, prefix: str,
                panther_pages: dict[str, str] | None = None) -> dict:
    """panther_pages: PANTHER tree page name (a MOD prefix) -> urlTemplate."""
    gene_prefix = gene_id.split(":", 1)[0]
    if prefix == "PANTHER" and page == "default" and gene_prefix in (panther_pages or {}):
        page, template = gene_prefix, panther_pages[gene_prefix]
    if curie == gene_id and page == "default":
        page = "gene"
    return {
        "geneId": gene_id,
        "geneTaxon": normalize_taxon(taxon),
        "xrefId": curie,
        "xrefUrl": xref_url(curie, template, prefix, gene_id),
        "xrefType": xref_type(page, prefix),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Alliance gene cross-references from the curation DB")
    parser.add_argument("--mods", help="Comma-separated MODs (default: all 9, env ALLIANCE_MODS)")
    parser.add_argument("--limit", type=int, help="Cap rows emitted")
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).parent.parent / "data"),
        help="Directory for the emitted TSV",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    release = get_current_release()
    mods = selected_mods(None, args.mods)
    wanted = set(resolve_gene_ids(None, ",".join(mods), None))
    out_dir = Path(args.out_dir)
    start = time.time()
    genes: set[str] = set()
    unlisted = 0
    panther_pages = dict(copy_rows(PANTHER_PAGES_SQL))

    with TsvWriter(
        out_dir / "gene-crossrefs.tsv",
        COLUMNS,
        release=release,
        source="curation DB dump (crossreference, resourcedescriptorpage)",
    ) as writer:
        per_taxon: dict[str, int] = {}
        for gene_id, taxon, curie, page, template, prefix in copy_rows(
                xref_sql([MOD_TAXON[m] for m in mods])):
            if gene_id not in wanted:
                unlisted += 1
                continue
            genes.add(gene_id)
            writer.write_row(**row_from_db(gene_id, taxon, curie, page, template, prefix,
                                           panther_pages))
            per_taxon[taxon] = per_taxon.get(taxon, 0) + 1
            if args.limit is not None and writer.rows_written >= args.limit:
                log.info("Hit limit; stopping")
                break
        for mod in mods:
            log.info("%s: %d rows", mod, per_taxon.get(MOD_TAXON[mod], 0))
        if writer.rows_written == 0:
            log.error("no rows written")
            sys.exit(1)

    log.info("Done. rows=%d genes=%d skipped_unlisted=%d in %.0fs",
             writer.rows_written, len(genes), unlisted, time.time() - start)


if __name__ == "__main__":
    main()
