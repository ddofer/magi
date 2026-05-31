# `figures` — publication plots

Matplotlib modules that reproduce analysis figures from notebooks and `../scripts/run_figures.py`. Figure **scripts** live here; data loaders for OMIA animals live in `../clinvar/animals_data.py`.

## Layout

| File | Role |
|------|------|
| `common.py` | Shared paths, label normalization, annotated parquet loaders |
| `plotting.py` | Core plot functions (KDE, concordance bars, mechanisms) |
| `merge_eval.py` | Merge canonical LLM results with signaled cohort + impact scores |
| `mechanism_canonicalize.py` | Regex buckets for free-text mechanism labels |
| `fig1_pies.py` | Label + rationale panel → `fig1c.png` |
| `benchmark_data.py` | Load combined SNP+indel benchmark deltas / impact frame |
| `fig2a_impact_vs_af.py` | Impact score vs gnomAD AF scatter → `fig2a.png` |
| `fig2_densities.py` | Global impact density panel → `fig2b.png` |
| `fig3_abc.py` | Variant-type concordance + Global_z quartiles → `fig3a`–`fig3c` |
| `fig3d_indel_frame.py` | Indel \|size\| mod 3 concordance → `fig3d.png` |
| `fig3e_mechanism_concordance.py` | NT mechanism vs LLM concordance → `fig3e.png`, `fig3e_indel.png` |
| `figs1_animals_concordance.py` | Supplementary Fig S1 — per-species SNP vs indel accuracy → `figs1.png` |

## Outputs

| Output | Source module |
|--------|----------------|
| `fig1c.png` | `fig1_pies.py` |
| `fig2a.png` | `fig2a_impact_vs_af.py` |
| `fig2b.png` | `fig2_densities.py` |
| `fig3a.png`–`fig3e*.png` | fig3 modules |
| `figs1.png` | `figs1_animals_concordance.py` (supplementary S1) |

Optional CSV sidecars are saved beside the matching PNG basename (e.g. `figs1.csv`).

## Inputs (via `config.PATHS`)

- **ClinVar cohort / LLM:** see [../clinvar/README.md](../clinvar/README.md)
- **Fig S1:** `data/llm/animals_*_evaluation_results*.parquet` + species labels from `omia_all_species_inferred_full.parquet` (notebook) or row-aligned `deltas_animals_{snp,indel}.parquet` fallback

## Run

```bash
cd clinvar_pipeline
python3 scripts/run_figures.py --fig s1
python3 scripts/run_figures.py --fig all
```

## Related docs

- [../ARTIFACTS.md](../ARTIFACTS.md)
- [../clinvar/README.md](../clinvar/README.md)
- [../output/README.md](../output/README.md)
