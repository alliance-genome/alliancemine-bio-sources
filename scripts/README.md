# Alliance data fetchers

Python scripts that pull data from upstream sources — the Alliance public API
(`https://www.alliancegenome.org/api`), the FMS bulk dumps at
`download.alliancegenome.org`, external services (STRING, ChEMBL, AlphaFold,
NCBI), and sibling InterMine instances via PathQuery — and emit TSV files in
the column shape the matching Java converter under `../alliance-*/src/main/
java/...` already knows how to read.

The Java layer stays thin: each `fetch_*.py` handles pagination, retry, rate
limiting, caching, and source-format quirks; the Java converter just parses a
stable TSV.

## Running

```
python3 scripts/fetch_all.py                              # every fetcher in order
python3 scripts/fetch_all.py --only genes,phenotypes      # subset
python3 scripts/fetch_all.py --skip crispr_screens,string # exclude slow ones
python3 scripts/fetch_all.py --limit 100                  # smoke test

# direct invocation
python3 scripts/fetch_interactions.py
python3 scripts/fetch_interactions.py --ids 'SGD:S000004103,SGD:S000002429'
python3 scripts/fetch_disease_alleles.py --source api    # force API fallback
python3 scripts/fetch_string.py --taxa 9606 --min-score 700
python3 scripts/fetch_crispr_screens.py --mods MGI,FB    # skip 753MB HUMAN
```

Each fetcher writes TSVs under `data/` relative to the repo root (override
with `--out-dir`). HTTP responses cache in `scripts/.cache/<fetcher>.sqlite`;
FMS gz payloads cache in `scripts/.fms_cache/`; STRING downloads cache in
`scripts/.string_cache/`.

`fetch_all.py` is the Docker pipeline's pre-build step — it runs every fetcher
in dependency order and exits non-zero if any fails.

## Environment

| Variable | Meaning | Default |
|---|---|---|
| `ALLIANCE_API_BASE` | API root | `https://www.alliancegenome.org/api` |
| `ALLIANCE_FMS_BASE` | FMS root | `https://fms.alliancegenome.org/api` |
| `ALLIANCE_FETCH_CACHE` | sqlite HTTP cache | `scripts/.cache/` |
| `ALLIANCE_FMS_CACHE` | FMS gz cache | `scripts/.fms_cache/` |

## Fetchers

### Primary

| Script | Source | Output | Notes |
|---|---|---|---|
| `fetch_genes.py` | API `/gene/{id}` | `data/alliance-genes.tsv` | Primary identifier source; everything downstream joins here |

### API-driven per-gene enrichment

| Script | Source | Output | Notes |
|---|---|---|---|
| `fetch_interactions.py` | API `/gene/{id}/molecular-interactions` + `/genetic-interactions` | `data/molecular-interactions.tsv` + `data/genetic-interactions.tsv` | API per-entity; WAF-rate-limited at ≥20 workers |
| `fetch_orthologs.py` | API `/gene/{id}/orthologs` | `data/orthologs.tsv` | |
| `fetch_paralogs.py` | API `/gene/{id}/paralogs` | `data/paralogs.tsv` | |
| `fetch_allele_detail.py` | API `/allele/{id}` | `data/allele-detail.tsv` | Enrichment-merge into alliance-alleles |
| `fetch_disease_annotations.py` | API `/disease/{id}/genes` | `data/disease-annotations-detail.tsv` | |
| `fetch_disease_models.py` | API `/gene/{id}/models` | `data/disease-models.tsv` | |
| `fetch_phenotypes.py` | API `/gene/{id}/phenotypes` | `data/phenotypes.tsv` | |

### FMS bulk (default `--source fms`)

These all accept `--source {fms,api}`. FMS path reads bulk JSON.gz dumps;
API path is the per-entity fallback. FMS is faster + immune to the WAF.

| Script | FMS datatype | Output | Notes |
|---|---|---|---|
| `fetch_transgenic_alleles.py` | `CONSTRUCT` + `ALLELE` (`contains` relation) | `data/transgenic-alleles.tsv` | Filter alleles with construct ref |
| `fetch_disease_alleles.py` | `DAF` (`objectType=allele`) | `data/disease-alleles.tsv` | |
| `fetch_experimental_disease.py` | `DAF` (`objectType=gene` + non-empty `primaryGeneticEntityIDs`) | `data/experimental-disease.tsv` | API endpoint perpetually `total:0`; FMS rescues |
| `fetch_variants.py` | `VARIANT-ALLELE-JSON` per taxon | `data/variants.tsv` | Filter `VariantId` non-null |
| `fetch_gene_descriptions.py` | `GENE-DESCRIPTION-JSON` | `data/gene-descriptions.tsv` | AGR auto-curated description text per gene |
| `fetch_agms.py` | `AGM` | `data/agms.tsv` | DiseaseModel + Strain (when `subtype=strain`) |
| `fetch_allele_phenotypes.py` | `PHENOTYPE` (filter `primaryGeneticEntityIDs` non-empty) | `data/allele-phenotypes.tsv` | |
| `fetch_htp.py` | `HTPDATASET` + `HTPDATASAMPLE` | `data/htp-datasets.tsv` + `data/htp-samples.tsv` | Populates Sample + SampleCharacteristic |
| `fetch_molecules.py` | `MOLECULE` | `data/molecules.tsv` | Stale at FMS 8.x; schema-ready for fresh upstream |

### Cross-references + external (dependency-ordered)

| Script | Source | Output | Notes |
|---|---|---|---|
| `fetch_gene_crossrefs.py` | FMS `GENECROSSREFERENCEJSON` | `data/gene-crossrefs.tsv` | Filters to MOD-subject rows; feeds alphafold |
| `fetch_alphafold.py` | URL template against `alphafold.ebi.ac.uk/entry/{acc}` | `data/alphafold.tsv` | **Depends on** `gene-crossrefs.tsv` (reads UniProtKB rows) |
| `fetch_string.py` | STRING DB v12 (`stringdb-downloads.org`) | `data/string-interactions.tsv` | Per-taxon downloads; default `--min-score 700` |
| `fetch_crispr_screens.py` | FMS `BIOGRID-ORCS` (streaming tar-extract) | `data/crispr-screens.tsv` | 753MB HUMAN tarball; HIT=YES filter |
| `fetch_chembl.py` | EBI FTP `chembl_uniprot_mapping.txt` | `data/chembl-targets.tsv` | 17k Protein cross-refs to ChEMBL target IDs |

### Cross-mine PathQuery REST

| Script | Source | Output |
|---|---|---|
| `fetch_mousemine_strains.py` | MouseMine PathQuery | `data/mouse-strains.tsv` |
| `fetch_wormmine_rnai.py` | WormMine PathQuery | `data/worm-rnai.tsv` |

## Enrichment-merge pattern

Many bio-source modules (the `-detail` suffix is a tell, but not all
enrichment modules end in `-detail`) emit **partial** InterMine items that
the integration engine merges into items from the primary source via shared
integration keys. The fetcher emits the same primary identifier so the merge
works without conflict.

Examples:
- `alliance-allele-detail` enriches `alliance-alleles` via `Allele.key_alleleid`
- `alliance-disease-detail` enriches `alliance-disease` via `DiseaseAnnotation.key_subject_term`
- `alliance-gene-descriptions` enriches existing Gene via `Gene.key_primaryidentifier`
- `alliance-allele-phenotypes` enriches existing PhenotypeAnnotation via `alleleSubject + ontologyTerm`
- `alliance-alphafold` enriches existing Protein via `primaryAccession`
- `alliance-chembl` emits CrossReference on Protein keyed by `subject + identifier`

## Column-schema coupling

Each fetcher's `COLUMNS` list and each Java converter's `COL_*` constants
must stay in lockstep. When adding a column, update both. The Java converter
silently skips rows with too few columns, so a mismatch is easy to miss —
verify with a smoke `--limit 5` run after any change.

## Shared helpers (`common.py`)

- `http_get_json(path, *, cache, stats)` — GET with retry + sqlite cache
- `paginate(path, *, limit, cache, stats)` — Alliance API `?limit=N&page=P`
- `iter_fms_datatype(datatype, *, sub_filter)` — stream FMS bulk shards
- `fms_list_files(datatype, *, sub_filter)` — manifest lookup
- `fms_download_gz(url)` — file-cache gz payloads
- `fms_iter_json_records(path, *, key)` — stream records from cached gz
- `intermine_paginate(base_url, query_xml)` — PathQuery REST against
  sibling InterMine instances (MouseMine / WormMine)
- `TsvWriter(path, columns, *, release, source)` — atomic header-commented
  TSV writer
- `open_cache(name)` — sqlite cache opener (WAL + busy_timeout for
  ThreadPool safety)
- `FetchStats`, `configure_logging`, `get_nested`, `join_pipe`,
  `normalize_taxon`, `get_current_release`,
  `enumerate_yeast_genes`, `enumerate_mod_genes`

## Outputs cumulative inventory (full integrate, 2026-05-13)

| TSV | Rows | Class populated |
|---|---|---|
| `alliance-genes.tsv` | ~360k | Gene (primary) |
| `molecular-interactions.tsv` | varies | Interaction |
| `genetic-interactions.tsv` | varies | Interaction |
| `orthologs.tsv` | varies | Homologue |
| `paralogs.tsv` | varies | Paralogue |
| `phenotypes.tsv` | varies | PhenotypeAnnotation (gene-keyed) |
| `disease-annotations-detail.tsv` | varies | DiseaseAnnotation (gene-keyed) |
| `disease-models.tsv` | varies | DiseaseModel + AlleleInteraction |
| `transgenic-alleles.tsv` | 13k MGI | TransgenicAllele |
| `disease-alleles.tsv` | 6k MGI | DiseaseAnnotation.alleleSubject |
| `experimental-disease.tsv` | 4.5k MGI | DiseaseAnnotation.experimentalGene |
| `variants.tsv` | 97k 9-taxa | Variant + VariantConsequence |
| `gene-descriptions.tsv` | 915k all MODs | Gene.autoDescription + 7 more attrs |
| `agms.tsv` | 159k all MODs | DiseaseModel + Strain (subtype=strain) |
| `allele-phenotypes.tsv` | 355k MGI | PhenotypeAnnotation (allele-keyed) |
| `htp-*.tsv` | 21k + 204k | DataSet + Sample + SampleCharacteristic |
| `gene-crossrefs.tsv` | 1.58M | CrossReference (on Gene) |
| `alphafold.tsv` | 262k | Protein.alphaFoldUrl |
| `string-interactions.tsv` | 208k yeast | STRINGInteraction |
| `crispr-screens.tsv` | 1.9M HUMAN+FB+MGI | CRISPRScreenResult |
| `molecules.tsv` | 6k WB | Molecule |
| `chembl-targets.tsv` | 17k | CrossReference (ChEMBL on Protein) |
