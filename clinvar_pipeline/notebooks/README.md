# Notebooks

| Notebook | Description |
|----------|-------------|
| [`variant_rationale_lookup.ipynb`](variant_rationale_lookup.ipynb) | Look up ClinVar submission rationales for any variant by locus |
| [`analysis_v2 animals.ipynb`](analysis_v2%20animals.ipynb) | Legacy reference for OMIA animal analysis (superseded by scripts 08–10) |
| [`animals_result_analysis.ipynb`](animals_result_analysis.ipynb) | OMIA LLM concordance plots (source for supplementary Fig S1) |

Run from this directory or from `clinvar_pipeline/` root.

## OMIA animals (use scripts, not the notebook)

The reproducible workflow lives in `../scripts/`:

```bash
python3 scripts/08_animals_extract_signals.py --variant-type snp
python3 scripts/08_animals_extract_signals.py --variant-type indel
python3 scripts/09_animals_build_prompts.py --variant-type snp
python3 scripts/09_animals_build_prompts.py --variant-type indel
# optional re-run:
# GEMINI_API_KEY=... python3 scripts/10_animals_run_llm_parallel.py --variant-type snp
python3 scripts/run_figures.py --fig s1
```

See [../ARTIFACTS.md](../ARTIFACTS.md) and [../clinvar/README.md](../clinvar/README.md) for inputs/outputs.
