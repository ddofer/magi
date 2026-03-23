# Setup and Usage Guide

## Quick Start

MAGI is a variant interpretation app that scores genomic variants using the NTv3 foundation model.

Current runtime configuration:

- Model: `InstaDeepAI/NTv3_650M_post`
- Sequence window: `32768` bp (32 kb)
- Batch limit: `10` variants
- Reference genome: local `data/hg38.fa` when available for human, otherwise UCSC / Ensembl fallback

## Hugging Face Spaces

1. Create a new Space at [Hugging Face Spaces](https://huggingface.co/spaces).
2. Select the `Gradio` SDK and `zero-a10` hardware.
3. Upload the contents of this directory.
   - **Important:** `data/hg38.fa` is ~915 MB. Do **not** commit it to the Space repo without Git LFS — instead rely on the UCSC API fallback (the app does this automatically when the file is absent).
   - `data/MANE_processed.csv` (~235 MB) and `data/MANE_processed.parquet` (~24 MB) should be committed; the parquet is auto-generated from the CSV on first run.
4. Set the `HF_TOKEN` Space secret so the app can access the gated `NTv3_650M_post` model weights.
5. Start the Space.

Notes:

- The app downloads model weights from Hugging Face on first launch and caches them.
- If `data/hg38.fa` is not present, sequence retrieval falls back to the UCSC API automatically — no extra configuration is needed.

## Local Installation

### Prerequisites

- Python 3.8+
- Enough disk space for the model cache and optional local genome
- Internet access for model download and, if needed, UCSC fallback sequence retrieval
- A Hugging Face account with access to `InstaDeepAI/NTv3_650M_post`

### Install dependencies

```bash
cd ntv3_gradio_app
pip install -r requirements.txt
```

### Optional: download the local genome

Automatic:

```bash
python download_hg38.py
```

Manual:

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mv hg38.fa data/hg38.fa
```

If the local genome is absent, the app still runs by querying UCSC for sequence windows.

### Configure Hugging Face access

The active NTv3 model is gated. If model loading fails with an authentication error:

1. Accept the model terms at [InstaDeepAI/NTv3_650M_post](https://huggingface.co/InstaDeepAI/NTv3_650M_post)
2. Set `HF_TOKEN` in your environment

Example:

```bash
export HF_TOKEN=your_token_here
```

### Validate installation

```bash
python test_installation.py
```

This checks:

- Python dependencies
- Required annotation and metadata files
- Sequence access through the local genome or UCSC fallback
- Model loading

### Run the app

```bash
python app.py
```

Or:

```bash
gradio app.py
```

## Input Format

### Single variant input

- Human chromosome: `chr1`-`chr22`, `chrX`, `chrY`, `chrM`
- Non-human chromosome: bare names such as `1`, `X`, `MT` or `chr`-prefixed names
- Position: 1-based coordinate for the selected species assembly
- `ref`: reference allele
- `alt`: alternate allele

Example:

- Chromosome: `chr17`
- Position: `7675088`
- Ref: `C`
- Alt: `T`

### Batch CSV input

Required columns:

```csv
chrom,pos,ref,alt
chr17,7675088,C,T
chr7,117559593,ATCT,A
chr13,32340300,G,A
```

Batch runs are limited to 10 variants.

## Output Overview

The single-variant view includes:

- Variant summary with gene, region class, and core metrics
- Ranked BED and BigWig signals
- Region track plot
- Full BED table
- Sequence-model metrics
- Rule-based interpretation panel

Important notes:

- `Global_z_sum_log` is a MAGI ranking score, not a pathogenicity probability.
- The simple high/moderate/low tier shown in the summary card is currently based on `Impact_Score_BED`.
- The interpretation panel is heuristic and deterministic. It summarizes existing outputs; it does not add a new predictive model.

## Troubleshooting

### Model download or authentication failed

- Confirm that you accepted the model terms for `InstaDeepAI/NTv3_650M_post`
- Confirm that `HF_TOKEN` is set
- Retry after verifying network access to Hugging Face

### No local genome found

- The app can run without `data/hg38.fa`
- In that case, sequence windows are requested from the UCSC API
- For faster local runs, download `hg38.fa` into `data/`

### Slow inference

- First launch is slower because weights may need to be downloaded
- CPU inference is much slower than GPU inference
- The 650M model requires more memory than smaller NTv3 variants

### Batch upload rejected

- Ensure the CSV contains `chrom`, `pos`, `ref`, and `alt`
- Ensure the file has no more than 10 rows

### Coordinates look wrong

- Human predictions expect GRCh38/hg38 coordinates
- Non-human predictions expect coordinates from the selected Ensembl assembly
- Non-human chromosome names can be bare or `chr`-prefixed

## Development Notes

Source-of-truth files for runtime behavior:

- `app.py`: UI, summaries, and batch handling
- `inference.py`: model loading, context length, and sequence retrieval
- `annotation.py`: region classes and annotation flags
- `analysis.py`: ranking and MAGI score computation
- `interpretation.py`: rule-based summary text

## References

- NTv3 model page: [InstaDeepAI/NTv3_650M_post](https://huggingface.co/InstaDeepAI/NTv3_650M_post)
- NTv3 paper: [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.01.11.523679v2)
- Gradio docs: [gradio.app/docs](https://www.gradio.app/docs)

**Last updated:** March 2026
