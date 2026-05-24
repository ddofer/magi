# `data` — inputs and canonical LLM results

Local copy of ClinVar / NT analysis inputs. **Not committed to git** except `README.md`, `.gitkeep` placeholders, and small threshold CSVs (see `../.gitignore`).

Populated by `../tools/sync_inputs.py` or manual copy.

## Layout

```
data/
├── metadata/              MANE GFF, MANE_processed, Promoter_processed, tracks metadata
├── clinvar_tables/        variant_summary, submission_summary (gz)
├── parquet/               clinvar_*_deltas.parquet (often symlinked)
├── AF/                    gnomAD allele-frequency subsets
├── thresholds_snps.csv
├── thresholds_indels.csv
├── uncertain_ids.csv      # VUS VariationIDs for fig1
└── llm/
    ├── snp_llm_results.parquet
    └── indel_llm_results.parquet
```

## Sync

```bash
python3 tools/sync_inputs.py --source /path/to/data/root --force

# Default: parent directory of clinvar_pipeline/ (symlinks large parquets ≥200 MB)
python3 tools/sync_inputs.py --force
```

## Canonical LLM parquets

Built by `tools/gather_canonical_llm.py` from parent exports (`parquet/snp_evaluation_results_v2.csv`, `llm_results/indel_results.parquet`, etc.). These are trimmed to columns used by figure scripts (`clinvar/llm_results.py`), not raw API responses.

```bash
python3 tools/gather_canonical_llm.py --force
```

## Related docs

- [../ARTIFACTS.md](../ARTIFACTS.md)
- [../output/README.md](../output/README.md)
