# ClinVar pipeline

End-to-end analysis pipeline for ClinVar variants scored with **NT-v3**: cohort preparation, signal extraction, impact scoring, LLM-based concordance evaluation, and publication figures. This module is part of the [MAGI](https://github.com/ddofer/magi) project.

## Overview

The pipeline processes labelled SNP and indel variants through a staged parquet workflow:

1. **Prepare** — ClinVar mapping, review-status filtering, MANE annotation  
2. **Signal** — threshold-gated top-k genomic track summaries  
3. **Evaluate** — optional LLM rationale judging (parallel or batch)  
4. **Score** — benign-referenced global impact metrics  
5. **Figure** — reproducible matplotlib outputs (cohort pies, impact densities, concordance panels)

Core logic lives in the `clinvar/` Python package; orchestration scripts are under `scripts/`.

## Requirements

- Python 3.10+
- Dependencies in [`requirements.txt`](requirements.txt) (NumPy, Pandas, PyArrow, Matplotlib, scikit-learn, SciPy, Seaborn; optional Google GenAI for LLM stages)

## Installation

```bash
cd clinvar_pipeline
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Data

Large inputs are **not** stored in this repository (see [`.gitignore`](../.gitignore)). Populate `data/` before running the pipeline:

| Input | Description |
|-------|-------------|
| `data/parquet/clinvar_new_deltas.parquet` | NT-v3 SNP delta predictions |
| `data/parquet/clinvar_indel_deltas.parquet` | NT-v3 indel delta predictions |
| `data/clinvar_tables/variant_summary.txt.gz` | ClinVar variant summary (GRCh38) |
| `data/clinvar_tables/submission_summary.txt.gz` | ClinVar submission summaries |
| `data/metadata/` | MANE GFF, processed region tables, track metadata |
| `data/AF/` | gnomAD allele-frequency subsets |
| `data/thresholds_*.csv` | Signal threshold tables (small; may ship with the repo) |
| `data/llm/` | Canonical LLM judge results (built by `tools/gather_canonical_llm.py`) |

Sync from an existing data tree:

```bash
python3 tools/sync_inputs.py --source /path/to/data/root --force
```

If `--source` is omitted, the tool looks for files in the parent directory of `clinvar_pipeline/`. See [data/README.md](data/README.md) for the full layout.

Generated artifacts are written to `output/` (also gitignored).

## Usage

```bash
# Metadata (once)
python3 scripts/00_build_mane_metadata.py

# Cohort + annotation
python3 scripts/01_prepare_variants.py
python3 scripts/02_extract_signals.py

# Impact scores + figures
python3 scripts/06_compute_impact_scores.py
python3 tools/gather_canonical_llm.py
python3 scripts/07_assign_mechanisms.py
python3 scripts/run_figures.py --fig all
```

Optional LLM re-runs: `scripts/04_run_llm_parallel.py`, `scripts/05_run_llm_batch.py`.

**Rule table (supplement).** Run `python3 scripts/export_rule_table.py [--out PATH]`; `scripts/07_assign_mechanisms.py` also writes it. Both files are generated from the rule definitions in `vendor/nt_mechanism_assignment.py`, so a rule change only needs a rerun.
- `output/mechanisms/rule_table_paper.csv`: the supplement table. One row per rule, in the order the rules are checked, with conditions written out ("Δ splice donor ≤ −0.2"). Rules that are the same for SNVs and indels appear once. The last row is the label for variants where no rule fires.
- `output/mechanisms/rule_table.csv`: one row per condition, for scripts.
- A copy of both files is tracked in `supplement/` (`--out supplement/rule_table.csv`), so coauthors can open them on GitHub. Rerun the export after any rule change.
- `can_be_assigned` is False for a rule that a rule checked earlier always pre-empts. The check compares conditions on the same feature only, so it can miss blocking through other features.

Draft caption for the paper table:

> **Rule layer.** Each rule tests Δ = Alt − Ref in NTv3's predicted score for an annotation element at the variant (a 0–1 sigmoid score, not a calibrated probability; for indels, the mean over the letters of the longer allele). Tier 5 rules use sequence-model metrics instead. A rule fires when at least one "Any of" condition and every "All of" condition hold; limits are inclusive, and a missing value counts as 0, so it never meets a condition. Rules are checked from the top and the first one that fires is assigned; all rules that fire are kept in NT_All_Mechanisms. A variant where no rule fires is labelled UNCERTAIN_SIGNIFICANCE, which is not a ClinVar VUS. LLR and the reference-letter probability come from the model's prediction at the variant letter in the reference sequence, without masking. The indel log-probability change is a z-score over all scored indels, before the ClinVar quality filter. UTR element names (5UTR+, 5UTR−, 3UTR+, 3UTR−) are NTv3's own. Mechanism names describe patterns in the model's predicted annotations, not observed molecular events.

**Rule provenance.** Mechanism assignment adds three columns next to `NT_Mechanism`:
- `NT_Tier`: tier of the assigned rule, or −1 when no rule fired.
- `NT_Rule`: the conditions that fired, with the observed values. A missing input prints as 0, as the rule sees it.
- `NT_Scored`: False when all of a variant's rule inputs are missing.

`analyses/nt_mechanism_assignment.py` is a symlink to `vendor/` (rebuttal scripts import it from there).

**Frozen benign background.** `compute_impact_frame(df, background="gradio_app/data/magi_baseline_stats.csv")` z-scores against the stored per-type benign mean and SD. Use it instead of the cohort's own benign rows for sets without labels, such as gnomAD or saturation mutagenesis.
- `tools/build_background.py SNV.parquet INDEL.parquet -o out.csv` rebuilds that file (population SD, NaN ignored). Run it through srun with about 4 CPUs and 48 GB.
- `scripts/06_compute_impact_scores.py --background PATH` uses it and writes `*_deltas_impact_frozen_bg.parquet`, so the paper figures keep reading the default files. Without the flag, script 06 uses the cohort's own benign rows, as in the paper.
- Rebuilt from the paper deltas, it matches the tracked file byte for byte.
- Background keys are lower-case ClinVar types: snp, deletion, insertion and indel.

Regenerate impact scores, mechanisms, and figures from existing cohorts, then validate:

```bash
python3 tools/run_validation.py
```

This command refreshes inputs and rewrites derived outputs. Run individual numbered stages to prepare cohorts or extract signals.

## Outputs

Publication figures are written to `output/figures/`:

| Figure | Description |
|--------|-------------|
| `fig1c.png` | Cohort label and rationale coverage |
| `fig2a.png` | MAGI impact score vs gnomAD AF |
| `fig2b.png` | Global impact density (SNP / indel / combined) |
| `fig3a`–`fig3e*.png` | LLM concordance and mechanism panels |

See [figures/README.md](figures/README.md) and [ARTIFACTS.md](ARTIFACTS.md).

## Repository layout

```
clinvar_pipeline/
├── clinvar/          # Core library
├── scripts/          # Pipeline stages 00–10, run_figures.py
├── figures/          # Publication figure modules
├── tools/            # Input sync, validation, LLM gather
├── vendor/           # Bundled mechanism-assignment logic
├── data/             # Inputs (local; not in git)
├── output/           # Generated artifacts (local; not in git)
├── config.py         # Paths and settings
└── ARTIFACTS.md      # Detailed stage reference
```

## Documentation

- [ARTIFACTS.md](ARTIFACTS.md) — pipeline stages and artifact model  
- [clinvar/README.md](clinvar/README.md) — library modules  
- [figures/README.md](figures/README.md) — figure scripts and inputs  
- [data/README.md](data/README.md) — input directory layout  

## License

See the license file in the [MAGI repository root](https://github.com/ddofer/magi).
