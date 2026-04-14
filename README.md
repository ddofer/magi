# MAGI Genome Analysis

This repository contains the research code used to run MAGI analyses built on the NTv3 foundation model for genomic variant interpretation. It is organized for paper support and collaborator use. The repository can remain private while still being cleanly structured for reproducible review.

## Repository Layout

```text
inference.py                    # Shared human ClinVar inference pipeline

inference_indel.py              # INDEL wrapper around inference.py
inference_omia.py               # OMIA multi-species inference pipeline
analysis.py                     # Downstream statistical analysis and plotting
smoke_test.py                   # Reduced-size validation run for inference.py
genomic_featImp_experiment_v3.ipynb
eda_consequence_score_distributions.ipynb
fig2b_impact_vs_af.ipynb
data/
  genomes/                      # Downloaded non-human genomes for OMIA runs
  results/                      # OMIA outputs
  stav_data/                    # Stav analysis inputs and derived files
ntv3_gradio_app/                # Gradio demo application and its bundled assets
```

## Main Entry Points

### 1. Human SNP inference

Run the main ClinVar SNP pipeline:

```bash
python inference.py
```

This pipeline:

- loads the NTv3 model and tokenizer
- extracts reference and alternate genomic sequence windows
- computes BED and BigWig delta features
- computes sequence-model features such as LLR, KL, and optional embeddings
- writes a parquet output table with metadata and model-derived features

Default input and output:

- input: `clinvar_input.parquet`
- output: `clinvar_new_deltas.parquet`

### 2. Human INDEL inference

Run the INDEL wrapper:

```bash
python inference_indel.py
```

This reuses the shared pipeline in `inference.py` but switches the input/output files and enables embeddings by default.

Default input and output:

- input: `clinvar_indel.parquet`
- output: `clinvar_indel_deltas.parquet`

### 3. OMIA multi-species inference

Run the animal-variant pipeline:

```bash
python inference_omia.py --validate-only
python inference_omia.py --species dog cat --sample-size 50
python inference_omia.py --full
```

This pipeline:

- fetches and caches OMIA single-locus variants
- normalizes them into the same ref/alt format used by the human pipelines
- downloads and indexes species genomes under `data/genomes/`
- runs NTv3 inference for supported species
- writes outputs under `data/results/`

Currently supported NTv3 animal species in the script include dog, cat, chicken, zebrafish, rat, and mouse.

### 4. Downstream analysis

Run the analysis helpers on an inference parquet:

```bash
python analysis.py --input clinvar_new_deltas.parquet
```

The analysis module provides:

- normalized impact scoring
- pathogenic vs benign comparisons
- per-feature significance testing
- fingerprint plots for single variants
- helper functions for case-study inspection

### 5. Gradio demo app

The web app lives in `ntv3_gradio_app/` and is intentionally kept separate from the top-level research pipelines.

Typical local startup:

```bash
cd ntv3_gradio_app
python app.py
```

See `ntv3_gradio_app/README.md` and `ntv3_gradio_app/SETUP_GUIDE.md` for app-specific details.

### 6. Manuscript case studies

To reproduce the paper-ready figures and captions for the manuscript case studies:

```bash
python analyses/generate_case_study_assets.py --case-id mat1a_gly336arg --case-id alox15b_rs9895916 --case-id col4a2_chr13_110492070_ga
python analyses/assemble_case_study_manuscript_package.py
```

This pipeline reads the targets from `analyses/case_study_manifest.csv`, runs MAGI inference, and exports high-resolution (600 DPI) fingerprint/region-track panels alongside JSON metadata, CSV tables, and review artifacts into `analyses/paper_case_studies/`. The assembly script then composites these into the final manuscript panels. Use `analyses/paper_case_studies/manuscript_figure_package.md` for caption templates and asset mappings.

## Data Layout Notes

- Human reference assets are kept at the repository root and inside `ntv3_gradio_app/data/` for the app-specific workflow.
- OMIA artifacts live under `data/`.
- Stav-related files now live under `data/stav_data/`.
- Some large data files are intentionally ignored in git through `.gitignore`.

## Expected Inputs

The human inference scripts expect parquet inputs with at least these columns:

- `chrom`
- `pos`
- `ref`
- `alt`
- `label`

Additional columns are preserved into the output when possible.

## Lightweight Validation

For a quick structural check after changes:

```bash
python smoke_test.py --mode snp --model-size 100M
python smoke_test.py --mode indel --model-size 100M
```

These are still model runs, but they are smaller than the default pipeline settings. For an even lighter check during refactoring, use module import and syntax validation without launching inference.

## Dependencies

Core Python dependencies used across the pipelines include:

- `torch`
- `transformers`
- `pandas`
- `numpy`
- `pyfaidx`
- `tqdm`
- `requests`

Some workflows also require:

- `samtools` for reference validation and FASTA indexing in OMIA workflows
- GPU access for practical full NTv3 inference speed

## Notes on NTv3 Outputs

- BED outputs are treated as multilabel probabilities after sigmoid conversion.
- BigWig outputs are also converted with sigmoid and are interpreted as independent tracks.
- For non-human OMIA runs, BigWig outputs are intentionally disabled because the NTv3 architecture does not provide meaningful non-human BigWig deltas.

## References

- NTv3 model: [InstaDeepAI on Hugging Face](https://huggingface.co/InstaDeepAI)
- ClinVar dataset source used in this repo: [songlab/clinvar](https://huggingface.co/datasets/songlab/clinvar)
- OMIA: [omia.org](https://www.omia.org/)

## Citation

If you use this repository for manuscript review or follow-up analysis, cite the associated MAGI manuscript and the NTv3 model paper.
