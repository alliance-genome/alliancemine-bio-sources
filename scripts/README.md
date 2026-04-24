# Alliance API fetchers

Python scripts that pull per-entity JSON from the Alliance API
(`https://www.alliancegenome.org/api`) and emit TSV files in the shape the
InterMine Java converters under `../alliance-*/src/main/java/...` already
know how to read.

The goal is to keep the Java layer thin — each `fetch_*.py` handles the
pagination, retry, rate limiting, and caching; the matching Java converter
just parses a stable TSV column schema. This is the replacement strategy as
FMS (the old file-distribution service) is retired.

## Running

```
python3 scripts/fetch_interactions.py          # yeast genome by default
python3 scripts/fetch_interactions.py --ids 'SGD:S000004103,SGD:S000002429'
python3 scripts/fetch_interactions.py --limit 100     # smoke test first 100 genes
```

Each fetcher writes its TSVs under `data/` relative to the repo root (override
with `--out-dir`). A SQLite cache lives at `scripts/.cache/<fetcher>.sqlite`
so reruns skip API calls that already succeeded.

## Environment

| Variable | Meaning | Default |
|---|---|---|
| `ALLIANCE_API_BASE` | API root | `https://www.alliancegenome.org/api` |
| `ALLIANCE_FETCH_CACHE` | where sqlite cache goes | `scripts/.cache/` |

## Fetchers

| Script | Endpoint(s) | Output |
|---|---|---|
| `fetch_interactions.py` | `/gene/{id}/molecular-interactions` + `/gene/{id}/genetic-interactions` | `data/molecular-interactions.tsv` + `data/genetic-interactions.tsv` |

More fetchers will land in subsequent phases (paralogs, phenotypes,
gene-metadata enrichment).
