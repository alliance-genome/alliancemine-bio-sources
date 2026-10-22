"""
SQL and value mapping shared by the allele extractors (extract_allele_detail,
extract_alleles, extract_transgenic_alleles, extract_variants).

The allele universe mirrors the curation API's allele summary documents, the
source of both /api/allele/{id} and the FMS VARIANT-ALLELE files:
AlleleDAO.getAllAlleleSummaryIds (not obsolete, not internal, full name not
'wild type') joined as in findAllelesForSummaryByIds (symbol with formatText,
taxon, data provider and its cross reference). The fetchers seeded alleles
from the release's FMS VARIANT-ALLELE AlleleIds; seed_allele_ids() returns
that set so an extractor can keep the same scope, or None for every summary
allele of the selected MODs.
"""

from __future__ import annotations

import logging
from typing import Iterable

from common import MOD_TAXON, enumerate_allele_ids, normalize_taxon

log = logging.getLogger("alliance.extract.alleles")

# MODs the FMS VARIANT-ALLELE generator writes (FileGeneratorConfig.VariantsAlleles).
VARIANT_ALLELE_MODS = ["FB", "MGI", "RGD", "SGD", "WB", "ZFIN"]

# Organization abbreviation of the allele's data provider -> fetcher MOD name.
# XB covers both Xenopus species, so its MOD comes from the taxon.
_PROVIDER_MOD = {"FB": "FB", "MGI": "MGI", "RGD": "RGD", "SGD": "SGD", "WB": "WB", "ZFIN": "ZFIN"}
_TAXON_MOD = {taxon: mod for mod, taxon in MOD_TAXON.items()}


def mod_for(provider: str, taxon: str) -> str:
    """Fetcher MOD name of an entity from its data provider and taxon curie."""
    return _PROVIDER_MOD.get(provider) or _TAXON_MOD.get(normalize_taxon(taxon), "")


def provider_filter(mods: Iterable[str]) -> str:
    """SQL list of organization abbreviations covering the given MODs."""
    providers = sorted({"XB" if m in ("XBXL", "XBXT") else m for m in mods if m != "HUMAN"})
    return ", ".join(f"'{p}'" for p in providers) or "''"


def summary_alleles_sql(mods: Iterable[str]) -> str:
    """CTE body: one row per summary allele of the MODs' data providers
    (id, pid, provider, taxon, symbol, symbol_format, xref_display).
    xref_display is the summary document's crossReference displayName: the
    data provider cross reference, kept only when the provider's resource
    descriptor has the matching allele/ or transgene/references page."""
    return f"""
        SELECT a.id, b.primaryexternalid AS pid, org.abbreviation AS provider,
               ot.curie AS taxon, s.displaytext AS symbol, s.formattext AS symbol_format,
               CASE WHEN refpage.id IS NOT NULL THEN cr.displayname ELSE '' END AS xref_display
        FROM allele a
        JOIN biologicalentity b ON b.id = a.id AND b.obsolete = false AND b.internal = false
        JOIN slotannotation s ON s.singleallele_id = a.id
             AND s.slotannotationtype = 'AlleleSymbolSlotAnnotation' AND s.formattext IS NOT NULL
        JOIN ontologyterm ot ON ot.id = b.taxon_id
        JOIN crossreference cr ON cr.id = b.dataprovidercrossreference_id
        JOIN resourcedescriptorpage rd ON rd.id = cr.resourcedescriptorpage_id
        JOIN organization org ON org.id = b.dataprovider_id
        LEFT JOIN LATERAL (
            SELECT p.id FROM resourcedescriptorpage p
            JOIN resourcedescriptor rdesc ON rdesc.id = p.resourcedescriptor_id
            WHERE rdesc.prefix = org.abbreviation
              AND p.name = CASE WHEN rd.name = 'transgene' THEN 'transgene/references'
                                ELSE 'allele/references' END
            LIMIT 1) refpage ON true
        WHERE org.abbreviation IN ({provider_filter(mods)})
          AND a.id IN (
              SELECT a2.id FROM allele a2
              LEFT JOIN slotannotation fn ON fn.singleallele_id = a2.id
                   AND fn.slotannotationtype = 'AlleleFullNameSlotAnnotation'
              WHERE fn.displaytext IS NULL OR fn.displaytext <> 'wild type')
    """


FIRST_GENE_SQL = """
    SELECT DISTINCT ON (aga.alleleassociationsubject_id)
           aga.alleleassociationsubject_id, be.primaryexternalid, sa.displaytext
    FROM allelegeneassociation aga
    JOIN vocabularyterm rel ON rel.id = aga.relation_id AND rel.name = 'is_allele_of'
    JOIN gene g ON g.id = aga.allelegeneassociationobject_id
    JOIN biologicalentity be ON be.id = g.id
    JOIN slotannotation sa ON sa.singlegene_id = g.id
         AND sa.slotannotationtype = 'GeneSymbolSlotAnnotation'
    WHERE aga.internal = false AND aga.obsolete = false
    ORDER BY aga.alleleassociationsubject_id, aga.id
"""


def summary_alleles(copy, mods: Iterable[str]) -> list[list[str]]:
    """Rows of summary_alleles_sql ordered by DB id. Each extractor joins
    the other tables in Python: without planner statistics a single SQL join
    over the restored tables picks nested loops that run for hours."""
    return list(copy(f"SELECT * FROM ({summary_alleles_sql(mods)}) al ORDER BY id"))


def first_genes(copy) -> dict[str, tuple[str, str]]:
    """Allele DB id -> (gene id, gene symbol) of its first live is_allele_of
    association, by association id (the summary document's alleleOfGene)."""
    return {allele_id: (gene_id, symbol) for allele_id, gene_id, symbol in copy(FIRST_GENE_SQL)}


# Variants counted on an allele summary document (AlleleDAO variant query):
# live allele-variant association, variant type, a curated location on an
# assembly component.
LIVE_ALLELE_VARIANTS_SQL = """
    SELECT DISTINCT ava.alleleassociationsubject_id AS allele_id, v.id AS variant_id
    FROM allelevariantassociation ava
    JOIN variant v ON v.id = ava.allelevariantassociationobject_id
    JOIN ontologyterm vt ON vt.id = v.varianttype_id
    JOIN curatedvariantgenomiclocation cvg ON cvg.variantassociationsubject_id = v.id
    JOIN assemblycomponent ac ON ac.id = cvg.variantgenomiclocationassociationobject_id
    WHERE ava.obsolete = false AND ava.internal = false
"""


def seed_allele_ids(mods: list[str], all_alleles: bool) -> set[str] | None:
    """AlleleIds of the release's FMS VARIANT-ALLELE files (the fetchers'
    scope), or None when every summary allele is wanted."""
    if all_alleles:
        return None
    seed = set(enumerate_allele_ids(mods))
    log.info("Seed: %d FMS VARIANT-ALLELE allele ids for %s", len(seed), ",".join(mods))
    return seed


def alteration_type(variant_count: int) -> str:
    """AlleleSummaryDocumentBuilder.setAlterationType."""
    if variant_count <= 0:
        return "allele"
    if variant_count == 1:
        return "allele with one variant"
    return "allele with multiple variants"


def fms_category(variant_count: int) -> str:
    """VARIANT-ALLELE Category column."""
    if variant_count <= 0:
        return "allele"
    if variant_count == 1:
        return "allele with 1 known variant"
    return f"allele with {variant_count} known variants"


def yes_dash(flag: bool) -> str:
    return "yes" if flag else "-"


def variant_symbol(chromosome: str, hgvs: str) -> str:
    """'NT_033779.5:g.18844476_18844477ins' on 2L -> '2L:18844476_18844477ins'."""
    colon = hgvs.find(":")
    if not chromosome or colon <= 0:
        return hgvs
    rest = hgvs[colon + 1:]
    if rest.startswith("g."):
        rest = rest[2:]
    return f"{chromosome}:{rest}"


def padded_sequences(type_name: str, ref: str, var: str, padded: str) -> tuple[str, str]:
    """Reference / variant sequence as the VARIANT-ALLELE files print them:
    the padded base leads both, and an insertion of unknown sequence shows
    '.' after it."""
    if not padded:
        return ref, var
    if type_name == "insertion" and not var:
        var = "."
    return padded + ref, padded + var


_SEP = "\x1f"

VARIANTS_SQL = f"""
    WITH xr AS (
        SELECT gx.genomicentity_id AS id,
               string_agg(cr.referencedcurie, '{_SEP}' ORDER BY cr.id) AS curies
        FROM genomicentity_crossreference gx
        JOIN crossreference cr ON cr.id = gx.crossreferences_id
        WHERE gx.genomicentity_id IN (SELECT id FROM variant)
        GROUP BY gx.genomicentity_id)
    SELECT v.id, cvg.id, vt.curie, vt.name, cvg.hgvs, cvg.start, cvg."end", ac.name,
           ga_be.primaryexternalid, cvg.referencesequence, cvg.variantsequence,
           cvg.paddedbase, xr.curies
    FROM variant v
    JOIN ontologyterm vt ON vt.id = v.varianttype_id
    JOIN curatedvariantgenomiclocation cvg ON cvg.variantassociationsubject_id = v.id
    JOIN assemblycomponent ac ON ac.id = cvg.variantgenomiclocationassociationobject_id
    LEFT JOIN genomeassembly ga ON ga.id = ac.genomeassembly_id
    LEFT JOIN biologicalentity ga_be ON ga_be.id = ga.id
    LEFT JOIN xr ON xr.id = v.id
    ORDER BY v.id, cvg.id
"""

# Gene-level predicted consequences (one per affected gene) carry the
# VARIANT-ALLELE MostSevereConsequenceName / VariantAffectedGeneId values.
CONSEQUENCES_SQL = """
    SELECT pvc.variantgenomiclocation_id, pvc.id, otc.name, otc.severityorder,
           be_g.primaryexternalid, sa_g.displaytext
    FROM predictedvariantconsequence pvc
    LEFT JOIN predictedvariantconsequence_ontologyterm pvco ON pvco.predictedvariantconsequence_id = pvc.id
    LEFT JOIN ontologyterm otc ON otc.id = pvco.vepconsequences_id
    LEFT JOIN transcriptgeneassociation tga ON tga.transcriptassociationsubject_id = pvc.varianttranscript_id
    LEFT JOIN biologicalentity be_g ON be_g.id = tga.transcriptgeneassociationobject_id
    LEFT JOIN slotannotation sa_g ON sa_g.singlegene_id = be_g.id
         AND sa_g.slotannotationtype = 'GeneSymbolSlotAnnotation'
    WHERE pvc.genelevelconsequence = true
    ORDER BY pvc.variantgenomiclocation_id, pvc.id
"""


REFERENCES_SQL = """
    SELECT vr.variant_id, vr.references_id, cr.referencedcurie
    FROM variant_reference vr
    JOIN reference_crossreference rc ON rc.reference_id = vr.references_id
    JOIN crossreference cr ON cr.id = rc.crossreferences_id
    ORDER BY vr.variant_id, vr.references_id, cr.id
"""

# Prefixes that identify a paper (VariantAlleleFileGenerator.PAPER_PREFIXES).
PAPER_PREFIXES = {"PMID", "FB", "MGI", "RGD", "SGD", "WB", "XB", "ZFIN"}


def reference_curie(curies: Iterable[str], mod_prefix: str) -> str:
    """One curie per reference: PMID, else the file MOD's curie, else any
    paper curie (DOI, PMCID and others are skipped)."""
    pmid = mod_match = any_paper = ""
    for curie in curies:
        prefix = curie.split(":", 1)[0] if ":" in curie else ""
        if prefix == "PMID":
            pmid = pmid or curie
        elif prefix in PAPER_PREFIXES:
            if prefix == mod_prefix:
                mod_match = mod_match or curie
            any_paper = any_paper or curie
    return pmid or mod_match or any_paper


def information_references(per_reference: list[list[str]], mod_prefix: str) -> list[str]:
    """Distinct chosen curies of a variant's references, in reference order."""
    chosen = (reference_curie(curies, mod_prefix) for curies in per_reference)
    return list(dict.fromkeys(c for c in chosen if c))


def _references_by_variant(rows: Iterable[list[str]]) -> dict[str, list[list[str]]]:
    out: dict[str, dict[str, list[str]]] = {}
    for variant_id, reference_id, curie in rows:
        out.setdefault(variant_id, {}).setdefault(reference_id, []).append(curie)
    return {vid: list(refs.values()) for vid, refs in out.items()}


def severity_joined(consequences: Iterable[tuple[str, int | None]]) -> str:
    """Distinct consequence names, most severe (lowest order) first, comma-joined;
    terms without an order rank last."""
    ranked = sorted(((o if o is not None else 1 << 30), i, n)
                    for i, (n, o) in enumerate(consequences) if n)
    return ",".join(dict.fromkeys(n for _, _, n in ranked))


def _consequences_by_location(rows: Iterable[list[str]]) -> dict[str, dict]:
    """location id -> {consequences: one comma-joined entry per gene-level
    predicted consequence, genes: [(id, symbol)] in first-seen order}."""
    terms: dict[str, dict[str, list]] = {}
    genes: dict[str, dict[str, str]] = {}
    for loc_id, pvc_id, name, order, gene_id, gene_symbol in rows:
        terms.setdefault(loc_id, {}).setdefault(pvc_id, []).append((name, to_int(order)))
        if gene_id:
            genes.setdefault(loc_id, {}).setdefault(gene_id, gene_symbol)
    out = {}
    for loc_id, per_pvc in terms.items():
        entries = [severity_joined(cons) for cons in per_pvc.values()]
        out[loc_id] = {
            "consequences": [e for e in entries if e],
            "genes": list(genes.get(loc_id, {}).items()),
        }
    return out


def load_variants(copy) -> dict[str, dict]:
    """Variant DB id -> fields of its first curated genomic location, in the
    VARIANT-ALLELE value formats. copy is curation_db.copy_rows."""
    by_location = _consequences_by_location(copy(CONSEQUENCES_SQL))
    references = _references_by_variant(copy(REFERENCES_SQL))
    variants: dict[str, dict] = {}
    for (vid, loc_id, type_id, type_name, hgvs, start, end, chrom, assembly,
         ref, var, padded, xrefs) in copy(VARIANTS_SQL):
        if vid in variants:
            continue
        ref, var = padded_sequences(type_name, ref, var, padded)
        cons = by_location.get(loc_id, {"consequences": [], "genes": []})
        variants[vid] = {
            "hgvs": hgvs,
            "symbol": variant_symbol(chrom, hgvs),
            "type_id": type_id,
            "type_name": type_name,
            "assembly": assembly,
            "chromosome": chrom,
            "start": start,
            "end": end,
            "ref": ref,
            "var": var,
            "consequences": cons["consequences"],
            "gene_ids": [g for g, _ in cons["genes"]],
            "gene_symbols": [s for _, s in cons["genes"]],
            "xrefs": split_agg(xrefs),
            "references": references.get(vid, []),
        }
    log.info("Loaded %d variants", len(variants))
    return variants


def allele_variant_links(copy) -> dict[str, list[str]]:
    """Allele DB id -> variant DB ids of its live variant links, ascending."""
    links: dict[str, list[str]] = {}
    for allele_id, variant_id in copy(LIVE_ALLELE_VARIANTS_SQL + " ORDER BY 1, 2"):
        links.setdefault(allele_id, []).append(variant_id)
    return links


def split_agg(value: str, sep: str = _SEP) -> list[str]:
    """Split a string_agg column joined with sep; '' -> []."""
    return [v for v in value.split(sep) if v] if value else []


def to_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
