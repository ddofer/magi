# ClinVar pipeline — layout overview

Refactored ClinVar / NT-v3 analysis pipeline. Lives as `clinvar_pipeline/` inside the [magi](https://github.com/ddofer/magi) repository (or standalone). Inputs go in `data/`; generated files go to `output/`.

## Repository tree

```
clinvar_pipeline/
├── config.py           # PATHS, impact columns, thresholds
├── data/               # Inputs + canonical LLM parquets (see data/README.md)
├── output/             # Generated parquets, figures, validation (see output/README.md)
├── clinvar/            # Python library (see clinvar/README.md)
├── scripts/            # Numbered pipeline stages (00–07) + run_figures.py
├── figures/            # Figure modules (see figures/README.md)
├── tools/              # sync_inputs, gather LLM, validation, compare_outputs
├── requirements.txt
└── ARTIFACTS.md        # This file
```

The parent ClinVar analysis checkout (or `--source` path for `tools/sync_inputs.py`) supplies large parquets and ClinVar tables. Mechanism assignment is bundled under `vendor/`.

## `data/` vs `output/`

| Directory | Contents | Typical refresh |
|-----------|----------|-----------------|
| `data/` | Symlinks/copies of deltas, ClinVar tables, thresholds, MANE metadata, `data/llm/*.parquet` | `tools/sync_inputs.py`, `tools/gather_canonical_llm.py` |
| `output/` | Strict/annotated/signaled parquets, impact scores, mechanisms, figures, validation JSON | Pipeline `scripts/` and `tools/run_validation.py` |

Keeping them separate allows refreshing inputs without deleting regenerated results.

## Pipeline stages (scripts)

| Step | Script | Main outputs |
|------|--------|----------------|
| 0 | `00_build_mane_metadata.py` | `data/metadata/MANE_processed.csv`, `Promoter_processed.csv` |
| 1 | `01_prepare_variants.py` | `output/parquet/*_strict.parquet`, `*_annotated.parquet`, `validation/cohort_funnel.json` |
| 2 | `02_extract_signals.py` | `output/parquet/*_signaled.parquet` |
| 3 | `03_build_prompts.py` | `output/parquet/*_with_prompts.parquet` |
| 4–5 | `04_run_llm_parallel.py`, `05_run_llm_batch.py` | Optional local LLM re-runs |
| 6 | `06_compute_impact_scores.py` | `output/impact/*_deltas_impact.parquet` |
| 7 | `07_assign_mechanisms.py` | `output/mechanisms/*_assigned_mechanisms.csv` |
| — | `run_figures.py` | `output/figures/fig1c.png`, `fig2a/b.png`, `fig3a`–`fig3e*.png` |

Stage 1 is memory-intensive (deferred BED/BW loading). Stages 6–7 and figures assume prior cohort and LLM artifacts exist.

## Cohort model (parquet overlap)

The rationaled ClinVar cohort flows through nested parquets with **the same variant rows** after the rationale filter:

1. **Quality** — gold-star, deduped variants with merged rationales (count only; `cohort_funnel.json`)
2. **Strict** — quality minus placeholder rationales; full track columns attached
3. **Annotated** — strict + MANE regions
4. **Signaled** — annotated + top-k signal summaries for LLM / mechanisms

Separate files exist so expensive steps can be re-run independently (e.g. re-extract signals without reloading deltas).

Impact parquets (`output/impact/`) cover the **full delta benchmark** for density plots; figure code merges scores onto the signaled cohort when needed.

## Publication figures

Flat outputs under `output/figures/`:

| File | Description |
|------|-------------|
| `fig1c.png` | Cohort label + rationale coverage (2×2 panel) |
| `fig2a.png` | MAGI impact score vs gnomAD AF (Global_z_sum_log scatter) |
| `fig2b.png` | SNP / indel / combined Global_z density (KDE) |
| `fig3a.png` | Concordance by variant type |
| `fig3b.png` | SNP concordance by Global_z quartile |
| `fig3c.png` | Indel concordance by Global_z quartile |
| `fig3d.png` | Indel concordance by \|size\| mod 3 |
| `fig3e.png` | SNP concordance by NT-assigned mechanism |
| `fig3e_indel.png` | Indel mechanism concordance |

All publication PNGs live directly under `output/figures/` (no subfolders).

## Tools

| Tool | Purpose |
|------|---------|
| `sync_inputs.py` | Copy/symlink parent-repo inputs into `data/` |
| `gather_canonical_llm.py` | Build `data/llm/snp_llm_results.parquet`, `indel_llm_results.parquet` |
| `compute_cohort_funnel.py` | Regenerate `cohort_funnel.json` without full track attach |
| `validate_mane.py` | MANE/Promoter validation vs parent metadata |
| `validate_figures.py` | Check figure inputs/outputs |
| `run_validation.py` | Chained sync → LLM gather → MANE → impact → mechanisms → figures |
| `compare_outputs.py` | Diff pipeline parquets vs legacy parent outputs |

## Environment

```bash
conda activate vep_env   # or equivalent env with requirements.txt
cd clinvar_pipeline
pip install -r requirements.txt
```

## Further reading

- [clinvar/README.md](clinvar/README.md) — library modules and data flow  
- [figures/README.md](figures/README.md) — figure scripts and inputs  
- [output/README.md](output/README.md) — generated paths and regeneration commands  
- [data/README.md](data/README.md) — input layout and sync
