# MAGI: Mechanistic Consequences of Genetic Variants via Genomic Foundation Models

Code for **MAGI**, which interprets genetic variants by comparing reference and alternate predictions from the NTv3 genomic foundation model (annotation elements, epigenomic/transcriptional tracks, sequence-model metrics) and mapping the signals to named molecular mechanisms.

Paper: [bioRxiv 2026.05.31.729117](https://www.biorxiv.org/content/10.64898/2026.05.31.729117v1.full)

## Repository layout

```text
inference.py            # Human ClinVar SNV inference (NTv3)
inference_indel.py      # Indel wrapper around inference.py
inference_omia.py       # OMIA multi-species (animal) inference
analysis.py             # Impact scoring, statistics, fingerprint plots
smoke_test.py           # Reduced-size validation run
analyses/               # Notebooks, figure scripts, generated tables, case studies
clinvar_pipeline/       # ClinVar cohort -> signals -> LLM concordance -> figures (see its README)
gradio_app/             # MAGI web app (Hugging Face Space; see its README)
paper_figures/          # Print-quality regeneration of manuscript panels (see RUNBOOK.md)
data/                   # OMIA outputs and small analysis inputs
LLM pipeline inputs/    # Inputs for the LLM rationale pipeline
```

Large files (reference genomes, delta parquets, model checkpoints) are not tracked; see `.gitignore`.

## Usage

### Human SNV inference

```bash
python inference.py          # clinvar_input.parquet -> clinvar_new_deltas.parquet
```

Loads NTv3, extracts reference/alternate windows, computes annotation-element and track deltas plus sequence-model metrics (LLR, KL, optional embeddings), and writes a parquet with metadata and features.

Input parquet needs at least `chrom`, `pos`, `ref`, `alt`, `label`; other columns are carried through.

### Indel inference

```bash
python inference_indel.py    # clinvar_indel.parquet -> clinvar_indel_deltas.parquet
```

### OMIA multi-species inference

```bash
python inference_omia.py --validate-only
python inference_omia.py --species dog cat --sample-size 50
python inference_omia.py --full
```

Supported species: dog, cat, chicken, zebrafish, rat, mouse. Genomes are downloaded to `data/genomes/`, outputs written to `data/results/`. Track (BigWig) outputs are disabled for non-human species.

### Downstream analysis

```bash
python analysis.py --input clinvar_new_deltas.parquet
```

### Web app

```bash
cd gradio_app
pip install -r requirements.txt
python app.py
```

### Case-study figures

```bash
python analyses/generate_case_study_assets.py --case-id mat1a_gly336arg --case-id alox15b_rs9895916 --case-id col4a2_chr13_110492070_ga
python analyses/assemble_case_study_manuscript_package.py
```

Targets come from `analyses/case_study_manifest.csv`; outputs go to `analyses/paper_case_studies/`.

### Quick check

```bash
python smoke_test.py --mode snp --model-size 100M
python smoke_test.py --mode indel --model-size 100M
```

## Dependencies

`torch`, `transformers`, `pandas`, `numpy`, `pyfaidx`, `tqdm`, `requests`. OMIA workflows also need `samtools`. A GPU is recommended for full NTv3 inference. LLM stages read `GEMINI_API_KEY` from the environment.

## Citation

If you use MAGI, please cite:

> Dan Ofer, Stav Zok, Michal Linial. MAGI: Mechanistic Consequences of Genetic Variants via Genomic Foundation Models. *bioRxiv* 2026.05.31.729117; doi: [10.64898/2026.05.31.729117](https://doi.org/10.64898/2026.05.31.729117)

```bibtex
@article{ofer2026magi,
  title   = {MAGI: Mechanistic Consequences of Genetic Variants via Genomic Foundation Models},
  author  = {Ofer, Dan and Zok, Stav and Linial, Michal},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.05.31.729117},
  url     = {https://www.biorxiv.org/content/10.64898/2026.05.31.729117v1}
}
```

Please also cite the NTv3 model ([InstaDeepAI on Hugging Face](https://huggingface.co/InstaDeepAI)).

## Data sources

- ClinVar benchmark: [songlab/clinvar](https://huggingface.co/datasets/songlab/clinvar)
- OMIA: [omia.org](https://www.omia.org/)
