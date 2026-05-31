# `data` — inputs and canonical LLM results

Local copy of ClinVar / NT analysis inputs. **Not committed to git** except `README.md`, `.gitkeep` placeholders, and small threshold CSVs (see `../.gitignore`).

Populated by `../tools/sync_inputs.py` or manual copy.

## Layout

```
data/
├── metadata/              MANE GFF, MANE_processed, Promoter_processed, tracks metadata
├── clinvar_tables/        variant_summary, submission_summary (gz)
├── parquet/               clinvar_*_deltas.parquet + deltas_animals_{snp,indel}.parquet
├── AF/                    gnomAD allele-frequency subsets
├── thresholds_snps.csv
├── thresholds_indels.csv
├── uncertain_ids.csv      # VUS VariationIDs for fig1
└── llm/
    ├── snp_llm_results.parquet
    ├── indel_llm_results.parquet
    ├── animals_snp_evaluation_results.parquet
    └── animals_indel_evaluation_results_2.parquet
```

## Sync

```bash
python3 tools/sync_inputs.py --source /path/to/data/root --force

# Default: parent directory of clinvar_pipeline/ (symlinks large parquets ≥200 MB)
python3 tools/sync_inputs.py --force
```

## Canonical LLM parquets

Built by `tools/gather_canonical_llm.py` from parent exports (`parquet/snp_evaluation_results_v2.csv`, `llm_results/indel_results.parquet`, etc.). These are trimmed to columns used by figure scripts (`clinvar/llm_results.py`), not raw API responses.

**OMIA animal LLM results** (`animals_*_evaluation_results*.parquet`) are produced by [`notebooks/analysis_v2 animals.ipynb`](../notebooks/analysis_v2%20animals.ipynb) and consumed by supplementary Fig S1 (`figures/figs1_animals_concordance.py`).

```bash
python3 tools/gather_canonical_llm.py --force
```

## Related docs

- [../ARTIFACTS.md](../ARTIFACTS.md)
- [../output/README.md](../output/README.md)
