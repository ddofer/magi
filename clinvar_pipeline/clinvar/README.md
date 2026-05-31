# `clinvar` — pipeline library

Python package for ClinVar variant preparation, NT delta integration, signal extraction, impact scoring, and LLM evaluation helpers. CLI entry points live in `../scripts/`; path configuration is in `../config.py`.

## Module map

| Module | Role |
|--------|------|
| `pipeline.py` | Orchestrates stage 1 (prep) and stage 2 (signal extraction) |
| `variant_prep.py` | Load delta parquets, ClinVar ID mapping, rationale merge, strict filter, deferred BED/BW attach |
| `clinvar_enrichment.py` | Submission-summary rationale aggregation and quality filtering |
| `rationale_lookup.py` | Ad-hoc ClinVar rationale lookup by locus (no cohort filter) |
| `region_annotation.py` | MANE-based genomic region labels on strict cohorts |
| `mane_build.py` | Build `MANE_processed.csv` / `Promoter_processed.csv` from MANE GFF |
| `signal_extraction_threshold.py` | Threshold-gated top-k BED/BW/MLM signals (default stage 2 path) |
| `signal_extraction.py` | Alternate signal extraction utilities |
| `impact_scoring.py` | Benign-referenced z-scores and global impact columns |
| `prompts.py` | Render per-variant LLM prompts from signaled rows |
| `llm_results.py` | Canonical column lists for gathered LLM parquets |
| `constants.py` | Shared column lists, review-status maps, feature priority |
| `llm/` | Parallel and batch LLM evaluation, response parsing |
| `prompts/` | System prompt templates (`system_prompt_snp.txt`, `system_prompt_indel.txt`, `system_prompt_animals_*.txt`) |
| `animals_pipeline.py` | OMIA animal deltas → signaled parquets (stage 8) |
| `animals_prompts.py` | OMIA animal per-variant prompts + system prompt loader |
| `animals_data.py` | Load merged OMIA animal LLM eval tables for Fig S1 |

## Data flow (stages 1–2)

```
NT delta parquet
    → map ClinVar IDs + gold stars
    → merge submission rationales
    → dedupe + drop placeholder rationales  →  strict parquet
    → MANE region annotation                →  annotated parquet
    → top-k signal extraction             →  signaled parquet
    → optional prompt build                 →  with_prompts parquet
```

Strict filtering removes variants whose combined rationale is a placeholder (`No rationale provided.`, etc.). Row counts typically drop from the **quality** cohort to the **rationaled** cohort; funnel stats are written to `../output/validation/cohort_funnel.json`.

## Memory model

Large delta parquets use **deferred track loading**: slim rows are filtered first, then BED/BW columns are attached in batches (`variant_prep.attach_track_columns`). This avoids loading thousands of track columns for variants that fail early filters.

## Impact scoring

`impact_scoring.py` implements the logic from `impact_score.ipynb`:

- Z-score BED/BW deltas against **Benign** rows, separately per `variant_type`
- Aggregate to `Global_z_sum_log`, `Composite_*`, etc.
- `attach_impact_to_signaled()` joins score columns onto signaled parquets at figure time (no duplicate merged parquet)

## External dependency

Mechanism assignment (`../scripts/07_assign_mechanisms.py`) imports `nt_mechanism_assignment.py` from the parent ClinVar repository root (hypothesis rules over BED/BW/MLM features).

## OMIA animals (stages 8–10)

```
data/parquet/deltas_animals_{snp,indel}.parquet
    → vectorized signal extraction (`animals_pipeline.extract_animals_signals`)
    → output/parquet/animals_*_signaled.parquet
    → prompt build (`animals_prompts.build_prompts_table`)
    → output/parquet/animals_*_with_prompts.parquet
    → optional Gemini eval (`scripts/10_animals_run_llm_parallel.py`)
    → data/llm/animals_*_evaluation_results*.parquet
    → Fig S1 per-species SNP vs indel accuracy (`figures/figs1_animals_concordance.py`)
```

System prompts: `prompts/system_prompt_animals_snp.txt`, `system_prompt_animals_indel.txt`.

## Related docs

- [../ARTIFACTS.md](../ARTIFACTS.md) — repository layout overview
- [../output/README.md](../output/README.md) — generated parquet and figure paths
- [../figures/README.md](../figures/README.md) — publication figure scripts
