# Drug Target Explorer

A data engineering solution for exploring OpenTargets drug-target associations.
Given a gene target of interest, returns clinical and marketed compounds acting on that target,
ranked by specificity and adverse-effect burden.

## Architecture

```
mlcs-coding-interview-main/
├── data/
│   ├── molecules.tsv          ← drug/molecule data with linked targets
│   └── adverseEffects.tsv     ← adverse event LLR data per drug
├── src/
│   ├── etl.py                 ← ETL pipeline: parse TSVs → SQLite
│   ├── api.py                 ← FastAPI REST API
│   └── query_client.py        ← CLI user script with table + chart output
├── results/                   ← output charts saved here
├── drug_targets.db            ← generated SQLite database
├── requirements.txt
└── README.md
```

## Quickstart

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Run ETL

Parses both TSV files, resolves Ensembl IDs to HGNC gene symbols via
biothings_client, and writes a normalised SQLite database.

```bash
python src/etl.py
```

Output — `drug_targets.db` with three tables:

- `drugs` — molecule metadata (name, approval status, clinical phase, etc.)
- `drug_targets` — M:N drug × target links with resolved gene symbols
- `adverse_effects` — one row per drug × adverse event (with LLR score)

### 3. Start the API

```bash
uvicorn src.api:app --reload --port 8000
```

Interactive docs available at: **http://127.0.0.1:8000/docs**

### 4. Query a target

In a second terminal (keep the API terminal open):

```bash
python src/query_client.py AR
python src/query_client.py DRD2 --top 10 --out results/
python src/query_client.py ENSG00000169083 --top 20 --out results/
```

Prints a ranked table in the terminal and saves a bubble chart to `results/`.

## API Reference

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/health` | Health check |
| GET | `/targets?limit=N` | All targets with drug counts |
| GET | `/drugs/{target}?top_n=N` | Ranked drugs (Ensembl ID or gene symbol) |
| GET | `/drugs/{target}/top` | Single best drug for a target |

Both Ensembl IDs (e.g. `ENSG00000169083`) and HGNC gene symbols (e.g. `AR`) are accepted interchangeably.

## Scoring Method

Each drug is ranked by a **composite score** (lower = better):

1. **Specificity** (`n_targets`) — number of distinct targets the drug acts on.
   Drugs hitting fewer targets are more specific to the queried gene.

2. **Adverse-effect burden** (`median_llr`) — median log-likelihood ratio (LLR)
   across all reported adverse events for that drug. Median is chosen as a robust
   aggregator that resists inflation from outlier high-LLR events.

3. **Clinical adjustments** — small bonuses for approved drugs and higher clinical
   phase; penalties for black box warnings and withdrawn compounds.

```
composite_score = rank(n_targets) + rank(median_llr)
                + approval_bonus  + phase_bonus
                + blackbox_penalty + withdrawn_penalty
```

Lower composite score = better drug recommendation.

## Visualization

`query_client.py` produces a bubble chart for each queried target saved to `results/`:

- **X-axis** — number of targets hit (specificity; fewer = more specific)
- **Y-axis** — median LLR (safety; lower = safer)
- **Bubble size** — inversely proportional to composite score (larger = better)
- **Color** — composite rank (green = best, red = worst)
- **Dashed lines** — 25th percentile thresholds marking the ideal zone (bottom-left)

## Design Decisions

- **SQLite** was chosen for simplicity and zero-dependency deployment. The schema
  is normalised into three tables with indices on all join keys for fast lookup.
- **FastAPI** provides automatic OpenAPI docs, Pydantic request validation,
  and clean separation between query logic and HTTP routing.
- **biothings_client** handles Ensembl to HGNC mapping in batched API calls,
  respecting rate limits automatically.
- **Median LLR** is preferred over mean or max because it is robust to the long
  tail of rare adverse events that have disproportionately high LLR values.
- **Column auto-detection** in the ETL means the pipeline adapts gracefully to
  minor schema changes in the source TSV files without requiring code edits.
