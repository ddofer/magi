# MAGI: Mechanistic Consequences of Genetic Variants via Genomic Foundation Models

Code for **MAGI**, which interprets genetic variants by comparing reference and alternate predictions from the NTv3 genomic foundation model (annotation elements, epigenomic/transcriptional tracks, sequence-model metrics) and mapping the signals to named molecular mechanisms.

Paper: [bioRxiv 2026.05.31.729117](https://www.biorxiv.org/content/10.64898/2026.05.31.729117v1.full)

## Repository layout

```text
inference.py            # Human ClinVar SNV inference (NTv3)
inference_indel.py      # Indel wrapper around inference.py
inference_omia.py       # OMIA multi-species (animal) inference
analysis.py             # Impact scoring, statistics, fingerprint plots
smoke_test.py           # Two-variant model check (writes to scratch/)
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

Input parquet needs at least `chrom`, `pos`, `ref`, `alt`; `label` and other columns are carried through.

Other variant sets use the same script with flags (`python inference.py --help`):

```bash
python inference.py --input my_variants.parquet --output my_deltas.parquet --compile --share-windows
```

- `--compile`: `torch.compile` the model. About 1.5x faster after a one-minute compile. Values differ from eager runs by about 0.3% on the largest deltas.
- `--share-windows`: nearby variants share one reference pass. Only variants at the same position within a 128-letter model token can share, because a 1-letter window shift changes NTv3's track deltas by about 37%. The saving is about 25% of passes for saturation mutagenesis, 14% for ClinVar indels and 2% for ClinVar SNVs. Shared variants sit off the window centre by whole tokens. That changed SNV features by about 1-4%, and moved 2 test indels by about 10%. Do not mix the two modes in one cohort.
- Results are written in chunks to `<output>.parts/`. A run that is stopped resumes from there; the folder is removed when the output is written.
- `inference.sbatch` runs this on a SLURM GPU, for example `sbatch inference.sbatch --input my_variants.parquet --output my_deltas.parquet`.

Two column groups come after the paper's features. They use the model's standard scales, not MAGI's sigmoid, and their prefixes keep them out of the impact scores:

- `LFC_*`: track change as a log2 fold change of the raw predicted signal, log2(alt+1) − log2(ref+1), over the selected tracks.
  - `LFC_max_abs` is the largest change, keeping its sign; `LFC_top_track` is its track.
  - `LFC_mean_abs` is the mean absolute change.
  - `LFC_max_abs_<group>` is the largest change per assay group (each histone mark, DNase, CAGE, RNA).
  - For an existing SNV output, `python lfc_from_deltas.py in.parquet out.parquet` adds the same columns without a GPU. Track values stored as 1.0 are too saturated to recover and are skipped.
- `WIN_<element>_gain` / `_loss`: the largest rise and drop in splice donor, splice acceptor, start codon and stop codon probability within 50 letters of the variant. `_gain_pos` / `_loss_pos` give the offset in reference letters (NaN when there is no rise or drop). This catches, for example, a splice site created a few letters away. For indels, letters after the variant are compared shifted by the length difference.

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

### Checks

Run the CPU regression tests without loading models or calling external services:

```bash
python -m unittest discover -s tests -v
```

Run two variants through a model (requires inference inputs and model access):

```bash
python smoke_test.py --mode snp --model-size 100M
python smoke_test.py --mode indel --model-size 100M
```

Smoke outputs stay under `scratch/smoke_snp/` or `scratch/smoke_indel/`. These checks use an 8K window; the indel run also exercises embeddings. They check output structure and probability bounds; they do not validate biological accuracy.

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
