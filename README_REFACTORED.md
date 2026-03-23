# NTv3 Genomic Variant Analysis

This repository contains a refactored implementation of the MAGI,  NTv3 (Nucleotide Transformer v3) based genomic variant analysis pipeline, split into two modular components for better organization and reusability.

## 📁 File Structure

```
nt3_inference.py      # Deep learning inference pipeline
nt3_analysis.py       # Downstream statistical analysis and visualization
nt3_genomic_featImp_experiment_v4.py  # Original combined script (kept for reference)
```

## 🔧 Components

### 1. `nt3_inference.py` - Inference Pipeline

Handles all deep learning model operations:

- Loading the NTv3 model and tokenizer
- Extracting genomic sequences (from local genome file or UCSC API)
- Running model inference on variants
- Computing prediction deltas (Alt - Ref)
- Calculating MLM (Masked Language Model) features:
  - **LLR** (Log-Likelihood Ratio): log(P(alt) / P(ref))
  - **MLM_Prior**: P(ref) from model
  - **MLM_Delta**: P(alt) - P(ref)
- Filtering BigWig tracks to meaningful subsets
- Saving results to parquet format

**Key Features:**

- ✅ BigWig track filtering (reduces ~7000 to ~3358 meaningful tracks)
- ✅ BED element predictions (21 genomic features)
- ✅ MLM features for variant effect prediction
- ✅ Fallback to UCSC API if genome file unavailable

**Usage:**

```bash
# Single GPU
python nt3_inference.py

# Multi-GPU with Accelerate
CUDA_VISIBLE_DEVICES=0,1 accelerate launch --num_processes 2 --mixed_precision bf16 nt3_inference.py

# Configuration via environment or editing script:
# DEBUG_MODE = True  (process first 70 variants)
# INPUT_DATA_FILE = "clinvar_input.parquet"
# OUTPUT_RESULTS_FILE = "clinvar_new_deltas.parquet"
```

**Output:**

- Parquet file with columns:
  - `chrom`, `pos`, `ref`, `alt`, `label`
  - `REF_BED_*` and `D_BED_*` (21 elements, probabilities)
  - `REF_BW_*` and `D_BW_*` (3358 tracks, logits)
  - `LLR`, `MLM_Prior`, `MLM_Delta`

### 2. `nt3_analysis.py` - Analysis Pipeline

Performs downstream analysis on inference results:

- Computing normalized impact scores
- Statistical comparisons (pathogenic vs benign)
- Per-element significance testing
- Visualization functions
- Case study analysis

**Key Features:**

- ✅ Z-score normalization for proper BED/BW scale handling
- ✅ Cohen's d effect sizes
- ✅ Bonferroni-corrected significance testing
- ✅ Variant fingerprint visualization
- ✅ Combined BED + BigWig plots
- ✅ Safe case study functions (handles missing variants)

**Usage:**

```bash
# Run analysis on inference results
python nt3_analysis.py --input clinvar_new_deltas.parquet

# Debug mode (first 70 variants)
python nt3_analysis.py --input clinvar_new_deltas.parquet --debug
```

**Output:**

- Statistical summaries printed to console
- Normalized impact scores computed
- Per-element analysis showing which genomic features differ between pathogenic/benign
- Example variant demonstrations

## 🚀 Complete Workflow

### Step 1: Prepare Data

```python
# Create clinvar_input.parquet with columns: chrom, pos, ref, alt, label
# label: True = Pathogenic, False = Benign
```

### Step 2: Run Inference

```bash
# Single GPU (faster for small datasets)
python nt3_inference.py

# Multi-GPU (recommended for full ClinVar ~40K variants)
CUDA_VISIBLE_DEVICES=0,1 accelerate launch --num_processes 2 --mixed_precision bf16 nt3_inference.py
```

### Step 3: Analyze Results

```bash
python nt3_analysis.py --input clinvar_new_deltas.parquet
```

## 📊 Key Outputs

### MLM Features

- **LLR (Log-Likelihood Ratio)**:
  - Negative = alt less likely than ref (potential loss-of-function)
  - Positive = alt more likely than ref (potential gain-of-function)
- **MLM_Prior**: Model's confidence in reference allele
- **MLM_Delta**: Direct change in probability

### BigWig Track Filtering

Reduces ~7000 tracks to ~3358 by keeping only:

- **Catlas** (kai*): 222 tracks
- **CAGE** (CNhs*): 1276 tracks  
- **GTEx**: 84 tissue-specific tracks
- **ENCODE** (ENCSR*): 1774 experiments

Filtered tracks focus on:

- Histone marks (H3K4me3, H3K27ac, H3K36me3, etc.)
- DNase/ATAC-seq (chromatin accessibility)
- Transcription factor binding
- Tissue-specific expression

### BED Elements (21 genomic features)

Probabilities (0-1 scale) for:

- protein_coding_gene, lncRNA
- exon, intron
- splice_donor, splice_acceptor
- promoter (tissue-specific/invariant)
- enhancer (tissue-specific/invariant)
- 5'UTR, 3'UTR
- CTCF-bound, polyA_signal
- start_codon, stop_codon, skipped_exon

## 🔬 Example Analysis Workflow

```python
import pandas as pd
from nt3_analysis import (
    compute_normalized_impact_scores,
    compare_pathogenic_benign,
    plot_variant_fingerprint,
    plot_combined_fingerprint
)

# Load results
df = pd.read_parquet("clinvar_new_deltas.parquet")

# Compute impact scores
df = compute_normalized_impact_scores(df, normalization='z_score')

# Statistical comparison
stats = compare_pathogenic_benign(df, score_type='Impact_Score_BW')
print(f"Cohen's d: {stats['cohens_d']:.3f}")
print(f"P-value: {stats['p_value']:.2e}")

# Visualize variant
plot_variant_fingerprint(df, variant_index=0, delta_prefix="D_BED_", top_k=15)
plot_combined_fingerprint(df, variant_index=0, top_k_each=8)
```

## ⚙️ Configuration

### Inference Configuration (`nt3_inference.py`)

```python
DEBUG_MODE = False  # Process first 70 variants only
CONTEXT_LEN = 64*1024  # Sequence context window
BATCH_SIZE = 4  # Inference batch size
USE_BED = True  # Include BED element predictions
USE_BIGWIGS = True  # Include BigWig track predictions

# File paths
INPUT_DATA_FILE = "clinvar_input.parquet"
OUTPUT_RESULTS_FILE = "clinvar_new_deltas.parquet"
GENOME_FILE = "hg38.2bit"  # Local genome (optional)
METADATA_FILE = "functional_tracks_metadata_human.csv"
```

### Analysis Configuration (`nt3_analysis.py`)

```python
# Normalization methods:
# - 'z_score': Z-score normalize BW tracks using benign distribution
# - 'mad': Use MAD (Median Absolute Deviation) for robust normalization
# - 'separate': Keep BED and BW scores separate

df = compute_normalized_impact_scores(df, normalization='z_score')
```

## 📈 Performance

- **ClinVar (SongLab dataset ~40K variants)**: ~8-10 hours on dual GPUs, 64K context window

## 🐛 Troubleshooting

### "species_ids must be >= 6" Error

The model requires species ID >= 6. The script automatically handles this using `model.encode_species(["human"])`.

### Missing BED Tracks

Some NTv3 model versions may not include BED track predictions. The script safely handles this with `getattr(out, "bed_tracks_logits", None)`.

### Out of Memory

Reduce `BATCH_SIZE` in `nt3_inference.py` or use gradient checkpointing.

### UCSC API Slow

Download the hg38.2bit genome file for much faster sequence extraction:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
```

## 📚 References

- **NTv3 Model**: [InstaDeepAI/NTv3_100M_post](https://huggingface.co/InstaDeepAI/NTv3_100M_post)
- **ClinVar Dataset**: [songlab/clinvar](https://huggingface.co/datasets/songlab/clinvar)
- **ENCODE Project**: [encode-project.org](https://www.encodeproject.org/)
- **GTEx Portal**: [gtexportal.org](https://gtexportal.org/)

**New way:**

```bash
# Step 1: Inference
python nt3_inference.py

# Step 2: Analysis  
python nt3_analysis.py
```
