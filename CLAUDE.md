# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

AllianceMine bio-sources: a collection of **InterMine** data loaders that integrate biological data (genes, alleles, orthologs, ontologies, expression, interactions, etc.) from the Alliance of Genome Resources (AGR), SGD, and homology databases into a unified genomic data warehouse.

Each top-level directory (e.g. `alliance-genes`, `sgd`, `mpo`, `hpo`, `sgd-gff`) is an independent Gradle subproject — a **"bio-source"** that parses one dataset and produces InterMine `Item` XML consumed by a downstream mine build.

## Build / Test

Version pins live in `gradle.properties` (`imVersion`, `bioVersion` — currently `5.1.0`). These are read via `System.getProperty` in the root `build.gradle`, so Gradle invocations already pick them up from `gradle.properties`.

- **Build everything**: `./gradlew build`
- **Build one module**: `./gradlew :bio-source-alliance-genes:build` (note the `bio-source-` prefix — see `settings.gradle`)
- **Run all tests**: `./gradlew test` (runs with `ignoreFailures = true` and `forkEvery = 1`)
- **Run one module's tests**: `./gradlew :bio-source-sgd:test`
- **Aggregate report**: `./gradlew testReport` → `build/reports/allTests`
- **Clean build**: `./gradlew clean build`

Tests use `scanForTestClasses = false` and only include `**/*Test.class`, so test class filenames *must* end in `Test`.

## Architecture

### How a bio-source is wired up

Every module follows the same InterMine convention — understanding it means you can navigate any module quickly:

1. **`{module}.properties`** (at module root) declares the source type to the InterMine integration engine:
   - `have.file.custom.tgt = true` → custom file parser (most Alliance sources, e.g. `alliance-genes`)
   - `have.db.tgt = true` → reads from a SQL database (e.g. `sgd`, which has both `SgdConverter` and `SgdProcessor`)
   - `have.file.obo = true` → OBO ontology loader with metadata like `obo.ontology.name` (e.g. `mpo`, `hpo`, `pato`)
   - `converter.class` → fully qualified Java class under `org.intermine.bio.dataconversion`
2. **`src/main/java/org/intermine/bio/dataconversion/{Name}Converter.java`** implements the parser. Custom file converters extend `BioFileConverter`; DB converters follow a `Converter` + `Processor` split (see `sgd/`). Parsers read input with `FormattedTextParser`, emit `Item`s via the inherited `ItemWriter`, and typically cache chromosomes/genes/synonyms in local `Map`s to de-duplicate refs.
3. **`src/main/resources/{module}_additions.xml`** extends the genomic data model with classes/attributes/references specific to this source. This file is merged with the root `alliancemine-global_additions.xml` by the `biosource-dbmodel` Gradle plugin (see `bioSourceDBModelConfig.globalAdditionsFile` in root `build.gradle`).
4. **`src/main/resources/{module}_keys.properties`** defines the **primary keys** the InterMine integration engine uses to merge records across sources (e.g. `Gene.key_primaryidentifier = primaryIdentifier, organism`). Getting these wrong silently produces duplicate objects.
5. **`src/test/java/.../{Name}ConverterTest.java`** extends `ItemsTestCase` and compares produced items against a fixture XML in `src/test/resources`.

### Root-level wiring

- **`settings.gradle`** — registers every module with a `bio-source-` prefix (e.g. directory `alliance-genes/` → project `:bio-source-alliance-genes`). **Directories on disk without a matching entry here (e.g. `diopt-orthologs/`, `sgd-complementation/`, `psi-complexes/`, `psi-mi-ontology/`) are NOT in the build** — treat them as archived/WIP until added to `settings.gradle`.
- **`build.gradle`** — `subprojects {}` block applies the `dbmodel` + `biosource-dbmodel` InterMine Gradle plugins, pins Java 1.8, declares shared dependencies (`bio-core`, `intermine-integrate`, `commons-collections`, `junit 4.8.2`, `xmlunit`), and configures the global additions merge. Each module's own `build.gradle` is usually minimal — just a `sourceSets {}` override (including `build/gen` for generated code) and `processResources { from('.') { include ("*.properties")} }` to pick up `{module}.properties`.
- **`alliancemine-global_additions.xml`** — mine-wide model extensions (custom classes like `ARS`, `Allele`, `AlleleInteraction`, plus attribute additions to standard InterMine types). Edit with care: changes affect every module's generated model.

### Adding a new bio-source

1. Create the directory and mirror the layout of an existing similar module (file-based → copy `alliance-genes`, DB-based → copy `sgd`, ontology → copy `mpo`).
2. Add the module to `settings.gradle` using the `bio-source-{name}` naming convention.
3. Write the converter, `_additions.xml`, and `_keys.properties`. Declare the target type in `{module}.properties`.

### API-backed sources (`scripts/` + fetcher pattern)

The Alliance's FMS file-distribution service is being wound down. New data feeds — and existing feeds that need richer shapes than FMS ever surfaced — come via the **live API at `https://www.alliancegenome.org/api`**. The API is per-entity JSON (no bulk endpoint), so a direct Java-side HTTP fan-out would be the wrong place for retry / pagination / caching logic.

Pattern: **Python fetcher → intermediate TSV → existing Java `BioFileConverter`**. All the messy HTTP plumbing lives under `scripts/`; the Java converters keep the same TSV-reader shape they always had.

```
[Alliance API] ---> [scripts/fetch_*.py]  ---> [data/*.tsv]  ---> [Java BioFileConverter]
    per-entity        pagination, retry,         stable column     reads TSV, same pattern
    JSON              sqlite cache, TSV writer   schema per feed   as FMS-era converters
```

Shared helpers in `scripts/common.py`: `http_get_json` (retry + sqlite cache), `paginate`, `intermine_paginate` (PathQuery REST against sister InterMine instances), `TsvWriter` (atomic, header-commented), `open_cache`, `FetchStats`, `configure_logging`, `get_nested`, `join_pipe`, `get_current_release`, `enumerate_yeast_genes`, `enumerate_mod_genes`. Each fetcher is a thin layer on top — see `scripts/README.md` for the current catalogue and how to run them.

`fetch_all.py` is the orchestrator (`python3 scripts/fetch_all.py`, with `--only <names>` and `--limit N` for smoke testing) that the Docker pipeline should invoke as a single pre-build step. It runs the API-backed fetchers in order and exits non-zero on any failure.

**Not every fetcher hits the Alliance API.** Some MOD-curated entities (mouse strains, worm RNAi screens) live only in their MOD's own InterMine instance, so the corresponding fetchers run PathQuery REST against MouseMine / WormMine via `intermine_paginate`. They share the same TSV-output contract as the Alliance-API fetchers.

Today's fetcher-backed sources (each driven by a paired `fetch_*.py` script):

| Source module | Fetcher script | API endpoint(s) | Emitted TSV |
|---|---|---|---|
| `alliance-genes` (enriched) | `fetch_genes.py` | `/gene/{id}` | `alliance-genes.tsv` |
| `alliance-genetic-interactions` | `fetch_interactions.py` | `/gene/{id}/genetic-interactions` | `genetic-interactions.tsv` |
| `alliance-molecular-interactions` | `fetch_interactions.py` | `/gene/{id}/molecular-interactions` | `molecular-interactions.tsv` |
| `alliance-paralogs` | `fetch_paralogs.py` | `/gene/{id}/paralogs` | `paralogs.tsv` |
| `alliance-phenotypes` | `fetch_phenotypes.py` | `/gene/{id}/phenotypes` | `phenotypes.tsv` |
| `alliance-disease-models` | `fetch_disease_models.py` | `/gene/{id}/models` | `disease-models.tsv` |
| `alliance-allele-detail` | `fetch_allele_detail.py` | `/allele/{id}` | `allele-detail.tsv` |
| `alliance-ortholog-detail` | `fetch_orthologs.py` | `/gene/{id}/orthologs` | `orthologs.tsv` |
| `alliance-disease-detail` | `fetch_disease_annotations.py` | `/disease/{id}/genes` | `disease-annotations-detail.tsv` |
| `alliance-mouse-strains` | `fetch_mousemine_strains.py` | MouseMine PathQuery (`/service/query/results`) | `mouse-strains.tsv` |
| `alliance-worm-rnai` | `fetch_wormmine_rnai.py` | WormMine PathQuery (`/service/query/results`) | `worm-rnai.tsv` |

**Column-schema coupling**: each Python fetcher's `COLUMNS` list and each Java converter's `COL_*` constants must stay in lockstep. When adding a new column, update both.

**Enrichment-merge pattern**: the `*-detail` bio-source modules don't own their target class. They emit *partial* items keyed on a shared integration key (e.g. `Allele.key_alleleid`, `Homologue.key_pair`, `DiseaseAnnotation.key_subject_term`) so InterMine's integration engine merges them into the items produced by the primary source. This keeps each fetcher + converter focused on one API endpoint without forcing the primary converter to know about all enrichment columns.

**Seed gene-ID lists**: `scripts/common.enumerate_yeast_genes()` currently pulls from the FMS BGI SGD export (still functional). `enumerate_mod_genes(mods)` is the cross-MOD generalisation — same FMS-BGI shape, takes a list of MOD prefixes. Once FMS is fully deprecated both need a replacement source — a plausible future home is the API itself if a paginated `/gene` list lands, or MOD-specific exports (SGD publishes a `chromosomal_feature.tab`).

## Code Conventions

- **Java 8** (`sourceCompatibility = 1.8`), 4-space indentation, braces on same line.
- License header: `Copyright (C) 2002-YYYY AllianceMine` LGPL block at the top of every converter (copy from a neighboring file).
- Logging: Log4j (`intermine-resources` provides it), not `System.out`. The sanitization pass migrated every converter; new code should follow.
- Null/empty checks: use `StringUtils` from `commons-lang` (already imported in most converters).
- Identifier naming in model: `primaryIdentifier` + `secondaryIdentifier` (case-sensitive, matters for key resolution).

## Repository Hygiene

- `bin/`, `build/`, `.gradle/`, `scripts/.cache/` are generated — root `.gitignore` covers them; do not commit.
- `intermine.log` and `query.log` are runtime artifacts from local mine builds.
- The global git rule (from `~/.claude/CLAUDE.md`): **never add "Claude" references or co-authorship to commit messages**.

## Known cross-repo TODOs

These belong in sibling repos, not here, but future sessions working on AllianceMine-wide issues should be aware:

- **`alliancemine/project.xml` references `EXPRESSION-ALLIANCE_COMBINED.tsv`** which is header-only in FMS 9.0.0 (the populated data lives in the per-MOD files at `EXPRESSION-ALLIANCE_{MOD}_*.tsv`). Until that project.xml is updated, AllianceMine ingests zero expression data. Fix: point at the per-MOD glob, or (post-FMS) switch to a Python expression fetcher that emits the same schema.
- **`agr_intermine_builder/docker/alliancemine/`** Dockerfile currently drives the FMS pull via `legacy/alliancemine-unified/download_data.py`. Add a step to run `python3 scripts/fetch_all.py` before the gradle build so the new TSVs under `data/` exist by the time the converters run.
- **AllianceGenesConverter's input path**: `project.xml` currently points at `/root/data/genes/` (a locally-built TSV). The new `scripts/fetch_genes.py` emits `data/alliance-genes.tsv` with the same 14-column schema + 3 enrichment columns. `project.xml` needs updating when the migration lands.
