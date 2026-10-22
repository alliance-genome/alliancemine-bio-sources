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
ALLIANCE_RELEASE=9.2.0 python3 scripts/fetch_all.py --out-dir /root/data/api   # build run
python3 scripts/fetch_all.py --only genes,phenotypes      # subset
python3 scripts/fetch_all.py --skip crispr_screens,string # exclude slow ones
python3 scripts/fetch_all.py --mods SGD,MGI --limit 100 --allow-release-mismatch  # smoke test

# direct invocation
python3 scripts/fetch_interactions.py                     # all 9 MODs, 6 workers
python3 scripts/fetch_interactions.py --mods FB,WB --workers 4
python3 scripts/fetch_interactions.py --ids 'SGD:S000004103,HGNC:1100'
python3 scripts/fetch_disease_alleles.py --source api    # force API fallback
python3 scripts/fetch_string.py --taxa 9606 --min-score 700
python3 scripts/fetch_crispr_screens.py --mods MGI,FB    # skip 753MB HUMAN
python3 -m unittest discover -s scripts/tests             # unit tests (stdlib only)
```

Each fetcher writes TSVs under `data/` relative to the repo root (override
with `--out-dir`). HTTP responses cache in `scripts/.cache/<fetcher>.sqlite`;
FMS gz payloads cache in `scripts/.fms_cache/`; STRING downloads cache in
`scripts/.string_cache/`.

Cache keys include the Alliance release and API base, so a new release never
reads old responses. A run without `ALLIANCE_RELEASE` makes one extra
`/releaseInfo` call to learn the release (cached per process).

`fetch_all.py` is the Docker pipeline's pre-build step — it runs every fetcher
in dependency order and exits non-zero if any fails.

### Release pinning

Every fetcher works against one release: `ALLIANCE_RELEASE` (or
`fetch_all.py --release`). It selects the FMS files (most recent upload per
dataSubType in that release's manifest), the BGI files the gene seeds come
from, and the `# Alliance release:` line of each TSV header. Unpinned, a
fetcher follows `/api/releaseInfo`, then FMS `releaseversion/current` — never
the newest release by date, which FMS lists weeks before it is live.
`fetch_all.py` refuses to start (exit 2) when the pin differs from
`/api/releaseInfo`; `--allow-release-mismatch` overrides for smoke tests.

### Seeding and throughput

The per-entity API fetchers (`genes`, `interactions`, `orthologs`,
`paralogs`, `phenotypes`, `disease_models`, `allele_detail`,
`disease_annotations`, and the `--source api` paths) seed from the pinned
release: genes from the 9 MODs' BGI files (369K genes in 9.1.0; per-MOD lists
cached in `$ALLIANCE_FETCH_CACHE/gene-ids/<release>/`), alleles from the
per-taxon `VARIANT-ALLELE` TSVs (538K), DOIDs from `DISEASE-ALLIANCE`
COMBINED (6.7K). `--mods` / `ALLIANCE_MODS` narrow the MODs; `--ids` and
`--limit` still work.

They run `--workers` threads (default 6, `ALLIANCE_FETCH_WORKERS`); the AGR
WAF answers HTTP 405 under heavier load and those are retried with backoff.
Progress (done/total, entities/s, req/s, failures, throttled, ETA) is logged
every 1000 entities or 60 s. A failed entity does not stop the run: failed
ids (including HTTP 400/404) are written to
`<out-dir>/<fetcher>.failures.tsv`, and the fetcher exits non-zero when
failures exceed `--max-failure-rate` (default 0.5%) or no row was written;
in both cases the TSV is not written. Successful responses are
committed to the cache as they arrive, so a rerun resumes from it.

## Environment

| Variable | Meaning | Default |
|---|---|---|
| `ALLIANCE_API_BASE` | API root | `https://www.alliancegenome.org/api` |
| `ALLIANCE_FMS_BASE` | FMS root | `https://fms.alliancegenome.org/api` |
| `ALLIANCE_RELEASE` | pinned Alliance release | `/api/releaseInfo` |
| `ALLIANCE_MODS` | MOD subset (narrows each fetcher's default) | all 9 |
| `ALLIANCE_FETCH_WORKERS` | default `--workers` | `6` |
| `ALLIANCE_MAX_FAILURE_RATE` | default `--max-failure-rate` | `0.005` |
| `ALLIANCE_FETCH_CACHE` | sqlite HTTP cache + gene-id lists | `scripts/.cache/` |
| `ALLIANCE_FMS_CACHE` | FMS gz cache | `scripts/.fms_cache/` |
| `ALLIANCE_PATHQUERY_MAX_AGE` | seconds a cached PathQuery page stays valid | `604800` |

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

### Curation DB extractors (`fetch_all.py --source db`)

`extract_orthologs.py` and `extract_paralogs.py` write the same
`orthologs.tsv` and `paralogs.tsv` from SQL against a restored curation dump
(`CURATION_PSQL`, see `curation_db.py`). Each runs one query over the whole
table and keeps the rows whose subject is a seed gene.

- Orthologs: `genetogeneorthologygenerated.strictfilter` rows (the API labels
  them `stringent`); obsolete genes dropped.
- Paralogs: one row per unordered pair, the direction whose subject is the
  earliest seed (identity, similarity, length and rank are per subject).
- Method lists keep the join-table row order, the order the API returns.
- The 9.2.0 dump has no SGD orthology or paralogy rows. Seeds of a species
  with no database rows go to the API (fetcher row mapping and cache); SGD
  ortholog rows are mirrored for the other species. `--no-api-fallback`
  skips this.

### FMS bulk (default `--source fms`)

These all accept `--source {fms,api}`. FMS path reads bulk JSON.gz dumps;
API path is the per-entity fallback. FMS is faster + immune to the WAF.
A MOD with no file in the release manifest fails the fetcher; narrow
`--mods` or `ALLIANCE_MODS` to skip it.

| Script | FMS datatype | Output | Notes |
|---|---|---|---|
| `fetch_alleles.py` | `VARIANT-ALLELE-JSON` per taxon → legacy 27-col TSV | `data/alliance-alleles.tsv` | Replaces frozen 4.0.0 `VARIANT-ALLELE_COMBINED.tsv`; converter patched to skip `#` comment headers instead of fixed 15 lines |
| `fetch_disease.py` | `DISEASE-ALLIANCE/COMBINED` (14 shards) | `data/disease-alliance.tsv` | Plain TSV concatenation; drops per-shard headers |
| `fetch_expression.py` | `EXPRESSION-ALLIANCE/COMBINED` (13 shards) | `data/expression-alliance.tsv` | COMBINED shards populated (≥1.77M rows/shard verified 2026-05-13) |
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

### Curation DB extractors (`fetch_all.py --source db`)

`extract_<name>.py` writes the same TSV as `fetch_<name>.py` from SQL against
a restored curation dump (`CURATION_PSQL` or `CURATION_DB_URL`, see
`curation_db.py`). The restore must include `resourcedescriptor`, `resourcedescriptorpage` (URL
templates) and `reference_crossreference` (publication ids).

| Script | Replaces | Notes |
|---|---|---|
| `extract_genes.py` | `fetch_genes.py` | One query for all MODs; gene set = per-MOD gene lists; listed ids the DB lacks go through the API. `automatedDescription` comes from FMS GENE-DESCRIPTION-JSON (the dump has no automated_gene_description notes). `dateProduced` is filled from dateUpdated (the API leaves it empty). |
| `extract_gene_crossrefs.py` | `fetch_gene_crossrefs.py` | crossReferences + gcrpCrossReference per listed gene; page names mapped to the FMS `xrefType` spelling. Own id on `default` becomes `gene`; PANTHER families get the MOD tree page. |
| `extract_expression.py` | `fetch_expression.py` | `geneexpressionannotation`, `expressionpattern`, `anatomicalsite` and their crossreference links. SourceURL lists every MGI assay and WB expression pattern of a row (FMS keeps one). FB rows carry FBcv anatomy qualifiers. |
| `extract_htp.py` | `fetch_htp.py` | `htpexpressiondatasetannotation*`, `htpexpressiondatasetsampleannotation*`. anatomyStatement is the anatomy term name (the DB keeps no statement). dateAssigned is blank. MGI samples have no sample id and are skipped. |
| `extract_phenotypes.py` | gene/allele/AGM phenotype annotations grouped per (gene, statement); allele and AGM rows join on `inferredgene_id`; `reference_crossreference` for publications | `uniqueId` is `<gene>\|<statement>` (the API value is a random index id). 3,161 sampled rows (MGI, ZFIN, FB, WB, SGD, RGD, Xenbase, HGNC genes) match on every other column. |
| `extract_allele_phenotypes.py` | allele annotations (gene = inferred gene, else the allele's only gene) and AGM annotations (inferred allele, no gene) | FB, WB, ZFIN annotations have no phenotype terms in the curation DB, so those MODs still read FMS through `fetch_allele_phenotypes._iter_mod`. AGM asserted-allele links are empty (244 MGI rows lost). |
| `extract_agms.py` | `affectedgenomicmodel`, AGM full-name slot annotation (display text), subtype term | ZFIN subtype is `fish` (FMS: empty). The DB holds 566K AGMs vs 159K in FMS, including 49.6K more strains (WB 47.6K, FB 1.8K, RGD 184). |
| `extract_disease_annotations.py` | direct rows (viaOrthologyOrder 0) | via-orthology rows from `/disease/{id}/genes` (`--orthology api`, default; `none` skips them) |
| `extract_disease.py` | MOD-curated rows | Source=Alliance orthology rows from FMS DISEASE-ALLIANCE COMBINED (`--orthology fms`, default) |
| `extract_disease_models.py` | all rows | none |
| `extract_disease_alleles.py` | all rows | none |
| `extract_experimental_disease.py` | all rows | none |

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
