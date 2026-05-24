# `figures` — publication plots

Matplotlib modules that reproduce analysis figures from the parent notebooks (`fig1.ipynb`, `fig2_density_label.ipynb`, `fig3_a_b_c.ipynb`, `fig3d_indel_frame.ipynb`, `signal_mechanism_concordance_fig.ipynb`). Invoked by `../scripts/run_figures.py`.

## Layout

| File | Role |
|------|------|
| `common.py` | Shared paths, label normalization, annotated parquet loaders |
| `plotting.py` | Core plot functions ported from notebooks (KDE, concordance bars, mechanisms) |
| `merge_eval.py` | Merge canonical LLM results with signaled cohort + impact scores |
| `mechanism_canonicalize.py` | Regex buckets for free-text mechanism labels |
| `fig1_pies.py` | Label + rationale panel → `fig1c.png` |
| `benchmark_data.py` | Load combined SNP+indel benchmark deltas / impact frame |
| `fig2a_impact_vs_af.py` | Impact score vs gnomAD AF scatter → `fig2a.png` |

Reference notebook `fig2a_impact_vs_af.ipynb` (not in git; may contain large embedded outputs) lives locally beside the Python module.
| `fig2_densities.py` | Global impact density panel → `fig2b.png` |
| `fig3_abc.py` | Variant-type concordance + Global_z quartiles → `fig3a`–`fig3c` |
| `fig3d_indel_frame.py` | Indel \|size\| mod 3 concordance → `fig3d.png` |
| `fig3e_mechanism_concordance.py` | NT mechanism vs LLM concordance → `fig3e.png`, `fig3e_indel.png` |

## Outputs

Primary PNGs are written flat under `../output/figures/`:

| Output | Source module |
|--------|----------------|
| `fig1c.png` | `fig1_pies.py` |
| `fig2a.png` | `fig2a_impact_vs_af.py` |
| `fig2b.png` | `fig2_densities.py` |
| `fig3a.png` | `fig3_abc.py` |
| `fig3b.png` | `fig3_abc.py` (SNP Global_z quartiles) |
| `fig3c.png` | `fig3_abc.py` (indel Global_z quartiles) |
| `fig3d.png` | `fig3d_indel_frame.py` |
| `fig3e.png` | `fig3e_mechanism_concordance.py` (SNP) |
| `fig3e_indel.png` | `fig3e_mechanism_concordance.py` (indel) |

Optional CSV sidecars (proportions / counts) are saved beside the matching PNG basename (e.g. `fig3b.csv`, `fig3e_proportions.csv`).

## Inputs (via `config.PATHS`)

- **Raw benchmark deltas:** `data/parquet/clinvar_*_deltas.parquet` (fig2a impact scoring input)
- **Impact:** `output/impact/*_deltas_impact.parquet` (fig2a, fig2b)
- **Cohort:** `output/parquet/*_annotated.parquet`, `*_signaled.parquet` (fig1, mechanisms)
- **LLM:** `data/llm/snp_llm_results.parquet`, `indel_llm_results.parquet` (fig3)
- **Mechanisms:** `output/mechanisms/*_assigned_mechanisms.csv` (fig3e)
- **Funnel:** `output/validation/cohort_funnel.json` (fig1 rationale pies; auto-computed if missing)

`merge_eval.load_merged_eval_frame()` builds the combined SNP+indel frame used by fig3a–d: LLM judgments left-joined to impact-enriched signaled rows on `#VariationID` / merge keys.

## Run

```bash
cd clinvar_pipeline
python3 scripts/run_figures.py --fig 2a    # impact vs AF scatter
python3 scripts/run_figures.py --fig all
python3 scripts/run_figures.py --fig 3e --assign-mechanisms   # refresh mechanism CSVs first
```

Paths are defined in `../config.py` (`fig1c`, `fig2a`, `fig2b`, `fig3a`, …).

## Related docs

- [../ARTIFACTS.md](../ARTIFACTS.md)
- [../output/README.md](../output/README.md)
