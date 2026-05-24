# ClinVar pipeline (MAGI)

ClinVar / NT-v3 analysis pipeline for cohort preparation, impact scoring, LLM evaluation, and publication figures. Intended as a subfolder of the [magi](https://github.com/ddofer/magi) repository.

## Quick start

```bash
cd clinvar_pipeline
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Populate data/ (see below), then run stages:
python3 scripts/01_prepare_variants.py
python3 scripts/02_extract_signals.py --snp --indel
python3 scripts/06_compute_impact_scores.py
python3 tools/gather_canonical_llm.py
python3 scripts/run_figures.py --fig all
```

See [ARTIFACTS.md](ARTIFACTS.md) for the full stage list and directory layout.

## What is in git

This folder is **source code + small config only** (~few MB). The following are **not** committed (see `.gitignore`):

| Excluded | Typical size | How to obtain |
|----------|--------------|---------------|
| `data/parquet/*.parquet` | GB-scale | NT delta predictions + `tools/sync_inputs.py` |
| `data/clinvar_tables/*.gz` | hundreds of MB | ClinVar release files |
| `data/metadata/*.gff.gz`, processed CSVs | tens of MB | MANE GFF + `scripts/00_build_mane_metadata.py` |
| `data/llm/*.parquet` | varies | `tools/gather_canonical_llm.py` |
| `output/` | multi-GB when built | Run pipeline scripts |
| `.venv/` | ~650 MB | Local `pip install` |

Threshold CSVs (`data/thresholds_*.csv`) are kept in git when present.

## Data setup

Place inputs under `data/` following [data/README.md](data/README.md).

If you have a local ClinVar / NT data tree (e.g. from a private analysis checkout), sync into `data/`:

```bash
python3 tools/sync_inputs.py --source /path/to/clinvar/data/root --force
```

Without `--source`, `sync_inputs.py` looks for files in the **parent of `clinvar_pipeline/`** (the magi repo root or a sibling ClinVar checkout).

Required large files (minimum to run impact + figures):

- `data/parquet/clinvar_new_deltas.parquet`
- `data/parquet/clinvar_indel_deltas.parquet`
- `data/clinvar_tables/variant_summary.txt.gz`
- `data/clinvar_tables/submission_summary.txt.gz`

## Before pushing to GitHub

From the **magi repo root** (or inside `clinvar_pipeline/`):

```bash
bash clinvar_pipeline/tools/check_large_files.sh
git add clinvar_pipeline
git status
```

The check script fails if any tracked file exceeds 50 MB. Your local `data/` and `output/` should stay untracked thanks to `.gitignore`.

## Uploading to `ddofer/magi`

```bash
git clone https://github.com/ddofer/magi.git
cd magi
# copy or rsync your clinvar_pipeline/ folder here (exclude .venv and output/)
rsync -a --exclude '.venv' --exclude 'output' --exclude '__pycache__' \
  /path/to/clinvar_pipeline/ ./clinvar_pipeline/

bash clinvar_pipeline/tools/check_large_files.sh
git add clinvar_pipeline
git commit -m "Add clinvar_pipeline analysis submodule"
git push
```

## Layout

```
clinvar_pipeline/
├── clinvar/          # Core library
├── scripts/          # Pipeline stages 00–07 + run_figures.py
├── figures/          # Publication figure modules
├── tools/            # Sync, validation, LLM gather
├── vendor/           # Bundled nt_mechanism_assignment.py
├── data/             # Inputs (not in git; see data/README.md)
├── output/           # Generated artifacts (not in git)
├── config.py
├── requirements.txt
└── ARTIFACTS.md
```

## Docs

- [ARTIFACTS.md](ARTIFACTS.md) — pipeline overview
- [clinvar/README.md](clinvar/README.md) — library modules
- [figures/README.md](figures/README.md) — figure scripts
- [data/README.md](data/README.md) — input layout
- [output/README.md](output/README.md) — generated outputs
