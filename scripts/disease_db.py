"""
Disease annotations read from a restored Alliance curation database dump,
shared by the extract_disease*.py scripts.

annotations() yields one DiseaseAnnotation per non-obsolete, non-internal
curation row (gene, allele and AGM annotations), with its subject, DO term,
relation, ECO codes, qualifiers, conditions and reference. Each carries the
genes and alleles the curation rows link to it: an AGM annotation's inferred
and asserted genes and alleles, an allele annotation's inferred and asserted
genes. The Alliance indexer builds its gene- and allele-level rows from these
links.

An annotation names its reference by AGRKB curie; RefMap.load() reads the
reference's PubMed and MOD curies from reference_crossreference and falls
back to the AGRKB curie for a reference without cross-references.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from curation_db import copy_rows

log = logging.getLogger("alliance.disease_db")

# Curation relation -> gene-level relation the API reports. Model relations of
# an AGM annotation become is_implicated_in on its genes.
_GENE_RELATION = {
    "is_model_of": "is_implicated_in",
    "is_ameliorated_model_of": "is_implicated_in",
    "is_exacerbated_model_of": "is_implicated_in",
}
_NEGATED = {
    "is_implicated_in": "is_not_implicated_in",
    "is_marker_for": "is_not_marker_for",
    "is_model_of": "does_not_model",
}
# Condition relation types the DISEASE-ALLIANCE TSV reports as Modifier; the
# others (has_condition, induced_by, ...) go to ExperimentalCondition.
MODIFIER_CONDITIONS = {"ameliorated_by", "exacerbated_by", "not_ameliorated_by", "not_exacerbated_by"}

SUBJECT = "subject"


@dataclass
class Entity:
    id: str
    taxon: str
    symbol: str
    kind: str  # gene | allele | agm


@dataclass
class Annotation:
    id: str
    kind: str  # gene | allele | agm
    subject: Entity
    disease_id: str
    disease_name: str
    relation: str
    negated: bool
    evidence_codes: list[str]
    reference: str
    provider: str
    secondary_provider: str
    date: str
    qualifiers: list[str]
    with_genes: list[str]
    conditions: list[tuple[str, str]]
    curie: str  # AGRKB annotation curie
    # (link kind, Entity): inferred_gene, asserted_gene, inferred_allele, asserted_allele
    links: list[tuple[str, Entity]] = field(default_factory=list)

    def genes(self) -> list[tuple[str, Entity]]:
        """Genes the annotation is about: its subject when it is a gene,
        else the inferred and asserted genes (each gene once)."""
        if self.kind == "gene":
            return [(SUBJECT, self.subject)]
        return _unique_links(self.links, ("inferred_gene", "asserted_gene"))

    def alleles(self) -> list[tuple[str, Entity]]:
        if self.kind == "allele":
            return [(SUBJECT, self.subject)]
        return _unique_links(self.links, ("inferred_allele", "asserted_allele"))


def _unique_links(links: list[tuple[str, Entity]], kinds: tuple[str, ...]) -> list[tuple[str, Entity]]:
    seen: set[str] = set()
    out = []
    for kind in kinds:
        for link_kind, ent in links:
            if link_kind == kind and ent.id not in seen:
                seen.add(ent.id)
                out.append((link_kind, ent))
    return out


def gene_relation(relation: str, negated: bool) -> str:
    """Relation string the API reports on the gene-level row of an annotation."""
    base = _GENE_RELATION.get(relation, relation)
    return negated_relation(base) if negated else base


def negated_relation(relation: str) -> str:
    return _NEGATED.get(relation, relation.replace("is_", "is_not_", 1) if relation.startswith("is_") else relation)


def yyyymmdd(timestamp: str) -> str:
    """'2012-01-06 18:16:47+00' -> '20120106'."""
    return timestamp[:10].replace("-", "") if timestamp else ""


def condition_label(relation_type: str) -> str:
    """'has_condition' -> 'Has Condition' as the DISEASE-ALLIANCE TSV prints it."""
    return " ".join(w.capitalize() for w in relation_type.split("_"))


def pmid_first(curies: list[str]) -> list[str]:
    """PubMed curies first, then MOD curies, each group sorted."""
    return sorted(curies, key=lambda c: (not c.startswith("PMID:"), c))


def split_pipe(value: str) -> list[str]:
    return [v for v in value.split("|") if v] if value else []


# ---------------------------------------------------------------------------
# References


class RefMap:
    """AGRKB reference curie -> PubMed / MOD publication curies, PMID first.
    Unmapped references fall back to the AGRKB curie."""

    def __init__(self, mapping: dict[str, list[str]] | None = None) -> None:
        self.mapping = mapping or {}
        self.misses: set[str] = set()

    @classmethod
    def load(cls) -> "RefMap":
        mapping = {ref: pmid_first(split_pipe(curies)) for ref, curies in copy_rows(REFERENCE_CURIES_SQL)}
        log.info("references with publication curies: %d", len(mapping))
        return cls(mapping)

    def curies(self, agrkb: str) -> list[str]:
        if not agrkb:
            return []
        hit = self.mapping.get(agrkb)
        if hit:
            return hit
        self.misses.add(agrkb)
        return [agrkb]

    def preferred(self, agrkb: str) -> str:
        """One curie for a reference: its PMID, else its first MOD curie."""
        curies = self.curies(agrkb)
        return curies[0] if curies else ""

    def log_misses(self, label: str) -> None:
        if self.misses:
            log.warning("%s: %d references have no publication cross-reference; written as AGRKB curies",
                        label, len(self.misses))


# ---------------------------------------------------------------------------
# SQL

# Entities (genes, alleles, AGMs) used by the annotations, with taxon and the
# symbol (displayText) each TSV prints: gene symbol, allele symbol, AGM full name.
_ENTITY_SQL = """
    SELECT b.id, b.primaryexternalid, t.curie AS taxon,
           CASE WHEN g.id IS NOT NULL THEN 'gene' WHEN a.id IS NOT NULL THEN 'allele'
                WHEN m.id IS NOT NULL THEN 'agm' ELSE '' END AS kind
    FROM biologicalentity b
    LEFT JOIN ontologyterm t ON t.id = b.taxon_id
    LEFT JOIN gene g ON g.id = b.id
    LEFT JOIN allele a ON a.id = b.id
    LEFT JOIN affectedgenomicmodel m ON m.id = b.id
    WHERE b.id IN (SELECT id FROM used) AND NOT b.obsolete AND NOT b.internal
"""

_USED_SQL = """
    used AS (
        SELECT diseaseannotationsubject_id AS id FROM genediseaseannotation
        UNION SELECT diseaseannotationsubject_id FROM allelediseaseannotation
        UNION SELECT diseaseannotationsubject_id FROM agmdiseaseannotation
        UNION SELECT inferredgene_id FROM agmdiseaseannotation
        UNION SELECT inferredallele_id FROM agmdiseaseannotation
        UNION SELECT inferredgene_id FROM allelediseaseannotation
        UNION SELECT assertedgenes_id FROM agmdiseaseannotation_gene
        UNION SELECT assertedalleles_id FROM agmdiseaseannotation_allele
        UNION SELECT assertedgenes_id FROM allelediseaseannotation_gene
        UNION SELECT with_id FROM diseaseannotation_gene
    )
"""

ENTITIES_SQL = f"""
    WITH {_USED_SQL},
    ent AS ({_ENTITY_SQL}),
    sym AS (
        SELECT s.singlegene_id AS id, s.displaytext AS symbol FROM slotannotation s
        WHERE s.slotannotationtype = 'GeneSymbolSlotAnnotation' AND NOT s.obsolete
          AND s.singlegene_id IN (SELECT id FROM used)
        UNION ALL
        SELECT s.singleallele_id, s.displaytext FROM slotannotation s
        WHERE s.slotannotationtype = 'AlleleSymbolSlotAnnotation' AND NOT s.obsolete
          AND s.singleallele_id IN (SELECT id FROM used)
        UNION ALL
        SELECT s.singleagm_id, s.displaytext FROM slotannotation s
        WHERE s.slotannotationtype = 'AgmFullNameSlotAnnotation' AND NOT s.obsolete
          AND s.singleagm_id IN (SELECT id FROM used)
    )
    SELECT ent.id, ent.primaryexternalid, ent.taxon, ent.kind, min(sym.symbol)
    FROM ent LEFT JOIN sym ON sym.id = ent.id
    GROUP BY 1, 2, 3, 4
"""

ANNOTATIONS_SQL = """
    WITH subj AS (
        SELECT id, diseaseannotationsubject_id AS subject, 'gene' AS kind FROM genediseaseannotation
        UNION ALL SELECT id, diseaseannotationsubject_id, 'allele' FROM allelediseaseannotation
        UNION ALL SELECT id, diseaseannotationsubject_id, 'agm' FROM agmdiseaseannotation
    ),
    eco AS (
        SELECT x.diseaseannotation_id AS id, string_agg(o.curie, '|' ORDER BY o.curie) AS codes
        FROM diseaseannotation_ontologyterm x JOIN ontologyterm o ON o.id = x.evidencecodes_id
        GROUP BY 1
    ),
    qual AS (
        SELECT x.diseaseannotation_id AS id, string_agg(v.name, '|' ORDER BY v.name) AS names
        FROM diseaseannotation_vocabularyterm x JOIN vocabularyterm v ON v.id = x.diseasequalifiers_id
        GROUP BY 1
    ),
    withg AS (
        SELECT x.diseaseannotation_id AS id, string_agg(x.with_id::text, '|' ORDER BY x.with_id) AS ids
        FROM diseaseannotation_gene x GROUP BY 1
    ),
    cond AS (
        SELECT x.annotation_id AS id,
               string_agg(coalesce(v.name, '') || '=' || coalesce(ec.conditionsummary, ''), '|'
                          ORDER BY v.name, ec.conditionsummary) AS conds
        FROM diseaseannotation_conditionrelation x
        JOIN conditionrelation cr ON cr.id = x.conditionrelations_id AND NOT cr.obsolete
        LEFT JOIN vocabularyterm v ON v.id = cr.conditionrelationtype_id
        JOIN conditionrelation_experimentalcondition y ON y.conditionrelation_id = cr.id
        JOIN experimentalcondition ec ON ec.id = y.conditions_id
        GROUP BY 1
    )
    SELECT d.id, s.kind, s.subject, obj.curie, obj.name, rel.name, d.negated::text,
           eco.codes, ref.curie, org.abbreviation, sorg.abbreviation,
           coalesce(d.datecreated, d.dateupdated)::text, qual.names, withg.ids, cond.conds, d.curie
    FROM diseaseannotation d
    JOIN subj s ON s.id = d.id
    JOIN ontologyterm obj ON obj.id = d.diseaseannotationobject_id
    JOIN vocabularyterm rel ON rel.id = d.relation_id
    LEFT JOIN informationcontententity ref ON ref.id = d.evidenceitem_id
    LEFT JOIN organization org ON org.id = d.dataprovider_id
    LEFT JOIN organization sorg ON sorg.id = d.secondarydataprovider_id
    LEFT JOIN eco ON eco.id = d.id
    LEFT JOIN qual ON qual.id = d.id
    LEFT JOIN withg ON withg.id = d.id
    LEFT JOIN cond ON cond.id = d.id
    WHERE NOT d.obsolete AND NOT d.internal
"""

LINKS_SQL = """
    SELECT id, 'inferred_gene', inferredgene_id FROM agmdiseaseannotation WHERE inferredgene_id IS NOT NULL
    UNION ALL SELECT id, 'inferred_allele', inferredallele_id FROM agmdiseaseannotation
        WHERE inferredallele_id IS NOT NULL
    UNION ALL SELECT id, 'inferred_gene', inferredgene_id FROM allelediseaseannotation
        WHERE inferredgene_id IS NOT NULL
    UNION ALL SELECT agmdiseaseannotation_id, 'asserted_gene', assertedgenes_id FROM agmdiseaseannotation_gene
        WHERE assertedgenes_id IS NOT NULL
    UNION ALL SELECT agmdiseaseannotation_id, 'asserted_allele', assertedalleles_id FROM agmdiseaseannotation_allele
        WHERE assertedalleles_id IS NOT NULL
    UNION ALL SELECT allelediseaseannotation_id, 'asserted_gene', assertedgenes_id FROM allelediseaseannotation_gene
        WHERE assertedgenes_id IS NOT NULL
"""

# PubMed and MOD publication prefixes, the curies the API lists in pubmedIds;
# DOI, PMCID and other cross-references are left out.
PUBLICATION_PREFIXES = ", ".join(f"'{p}'" for p in ("PMID", "MGI", "RGD", "FB", "SGD", "ZFIN", "WB", "Xenbase"))
REFERENCE_CURIES_SQL = f"""
    SELECT ref.curie,
           string_agg(DISTINCT x.referencedcurie, '|')
    FROM reference_crossreference rx
    JOIN informationcontententity ref ON ref.id = rx.reference_id
    JOIN crossreference x ON x.id = rx.crossreferences_id AND NOT x.obsolete
    WHERE rx.reference_id IN (SELECT evidenceitem_id FROM diseaseannotation)
      AND split_part(x.referencedcurie, ':', 1) IN ({PUBLICATION_PREFIXES})
    GROUP BY 1
"""

ECO_SQL = "SELECT curie, name, abbreviation FROM ontologyterm WHERE ontologytermtype = 'ECOTerm'"
SPECIES_SQL = "SELECT t.curie, s.fullname FROM species s JOIN ontologyterm t ON t.id = s.taxon_id"
DO_ANCESTORS_SQL = """
    SELECT s.curie, string_agg(DISTINCT o.curie, '|')
    FROM ontologytermclosure c
    JOIN ontologyterm s ON s.id = c.closuresubject_id AND s.ontologytermtype = 'DOTerm'
    JOIN ontologyterm o ON o.id = c.closureobject_id
    WHERE c.distance > 0
    GROUP BY 1
"""


def load_entities() -> dict[str, Entity]:
    ents = {}
    for db_id, pid, taxon, kind, symbol in copy_rows(ENTITIES_SQL):
        if pid:
            ents[db_id] = Entity(pid, taxon, symbol, kind)
    log.info("entities: %d", len(ents))
    return ents


def annotations(entities: dict[str, Entity] | None = None) -> list[Annotation]:
    """Every live disease annotation whose subject is a live entity."""
    ents = entities if entities is not None else load_entities()
    links: dict[str, list[tuple[str, Entity]]] = {}
    for ann_id, kind, ent_id in copy_rows(LINKS_SQL):
        ent = ents.get(ent_id)
        if ent:
            links.setdefault(ann_id, []).append((kind, ent))
    out = []
    dropped = 0
    for row in copy_rows(ANNOTATIONS_SQL):
        (ann_id, kind, subj_id, doid, doname, relation, negated, eco, ref, provider, sprovider,
         date, quals, with_ids, conds, curie) = row
        subject = ents.get(subj_id)
        if subject is None:
            dropped += 1
            continue
        out.append(Annotation(
            id=ann_id, kind=kind, subject=subject, disease_id=doid, disease_name=doname,
            relation=relation, negated=negated == "true", evidence_codes=split_pipe(eco),
            reference=ref, provider=provider, secondary_provider=sprovider, date=yyyymmdd(date),
            qualifiers=split_pipe(quals),
            with_genes=[ents[i].id for i in split_pipe(with_ids) if i in ents],
            conditions=[tuple(c.split("=", 1)) for c in split_pipe(conds)],
            curie=curie, links=links.get(ann_id, []),
        ))
    log.info("disease annotations: %d (%d with an obsolete or internal subject dropped)", len(out), dropped)
    return out


def eco_terms() -> dict[str, tuple[str, str]]:
    """ECO curie -> (name, abbreviation)."""
    return {curie: (name, abbr) for curie, name, abbr in copy_rows(ECO_SQL)}


def species_names() -> dict[str, str]:
    return dict(copy_rows(SPECIES_SQL))


def do_ancestors() -> dict[str, list[str]]:
    """DOID -> its ancestors over is_a and part_of."""
    return {doid: split_pipe(anc) for doid, anc in copy_rows(DO_ANCESTORS_SQL)}
