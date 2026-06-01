---
title: MAGI Variant Interpreter
emoji: 🧬
colorFrom: purple
colorTo: blue
sdk: gradio
sdk_version: 6.13.0
app_file: app.py
pinned: false
license: mit
suggested_hardware: zero-a10g
---

## 🧬 MAGI Variant Interpreter

### Variant impact scoring using MAGI and the NTv3 foundation model

MAGI is a lightweight demo for reviewing variant-associated signals in a local genomic window. It uses the external NTv3 model for sequence predictions and adds MAGI scoring, annotation, ranking, and short rule-based summaries.

MAGI was developed by Dan Ofer, Stav Zok, and Michal Linial. The NTv3 model was developed by InstaDeepAI and collaborators.

## Features

- **Single Variant Analysis**: Manual input of SNPs and indels with a compact summary
- **Batch Processing**: Upload CSV files with up to 10 variants
- **Multi-species support**: Human plus supported animals and plants via Ensembl sequence retrieval
- **Gene Annotation**: Automatic annotation using MANE Select RefSeq transcripts
- **Ranked Signal Review**:
  - BED outputs from the NTv3 configuration
  - Filtered BigWig context tracks
  - Sequence-model metrics such as LLR, KL divergence, and embedding distances
- **Region Track View**: Zoomable probability tracks for the top disrupted BED and BigWig outputs
- **Rule-Based Signal Interpretation**: Short deterministic summary of the strongest ranked signals
- **Impact Scoring**: Quantitative metrics for variant prioritization, including MAGI `Global_z_sum_log`

## Usage

### Single Variant

1. Select species and chromosome.
2. Enter a 1-based genomic position.
  Human uses GRCh38/hg38 coordinates. Non-human species use the selected species' current Ensembl assembly, and chromosome names can be bare (`1`, `X`, `MT`) or `chr`-prefixed.
3. Enter reference and alternate alleles (for example `C` → `T` for a SNP, `ATCT` → `A` for a deletion).
4. Click **Predict Impact**.

**Example:**

- Chromosome: `chr17`  
- Position: `7675088`  
- Ref: `C`  
- Alt: `T`  
- (TP53 pathogenic missense variant)

### Batch Upload

Upload a CSV file with these columns:

```csv
chrom,pos,ref,alt
chr17,7675088,C,T
chr7,117559593,ATCT,A
chr13,32332771,AGAGA,AGA
```

**Limit:** 10 variants per batch

## Output Interpretation

### Impact Scores

- **MAGI `Global_z_sum_log`**: Burden score computed as `Σ log(1 + |z_j|)` across z-scored BED and BigWig delta tracks using bundled baseline statistics

  - Higher values indicate broader or stronger deviation from the baseline set
  - This is a ranking score, not a calibrated pathogenicity probability

- **Impact_Score_BED**: Mean of top-3 largest absolute BED element deltas

  - This score is also used for the simple high/moderate/low summary tier shown in the single-variant card
  - Tier thresholds: `HIGH > 0.10`, `MODERATE > 0.05`, otherwise `LOW`

- **LLR (Log-Likelihood Ratio)**: For SNPs only, log(P(ALT)/P(REF))

  - Positive: Alternate allele more likely
  - Negative: Reference allele more likely

- **KL Divergence**: Measures how much the predicted token distribution changes around the variant

  - Higher values indicate a larger local sequence-model shift

### Rule-Based Signal Interpretation

The single-variant view includes a short interpretation block that:

- summarizes the most disrupted BED and BigWig signals in plain language
- highlights whether the strongest evidence is coding-related, splice-related, promoter-related, or context-dependent
- adds sequence-model context from LLR and KL divergence

This panel is deterministic and uses the same ranked signals already displayed in the table and plots. It is a heuristic summary, not a calibrated pathogenicity assessment.

The summary card also surfaces the top raw BED, BigWig, and MLM signals first, using a minimum absolute magnitude threshold of 0.03.

**MAGI baseline note:** the app bundles baseline statistics under `data/magi_baseline_stats.csv`, so `Global_z_sum_log` is computed locally.

### Functional Annotations

**Region classes:**

- `CODING`
- `CODING, SPLICE`
- `SPLICE`
- `UTR_5` / `UTR_3`
- `PROMOTER`
- `INTRONIC`
- `GENIC_OTHER`
- `OTHER`

**BigWig tracks:**

- Histone modifications: H3K4me3, H3K27ac, H3K36me3, H3K27me3, etc.
- Chromatin accessibility: ATAC-seq, DNase-seq
- Gene expression or transcription-linked assays such as CAGE

**Direction:**

- **Gain of Function (GOF)**: Δ > 0 → predicted signal increased
- **Loss of Function (LOF)**: Δ < 0 → predicted signal decreased

### Region Classification

Variants are automatically classified into:

- `CODING`: Overlaps coding sequence
- `CODING, SPLICE`: Coding region near splice junction
- `SPLICE`: Intronic splice site (±2 bp from exon boundary)
- `UTR_5` / `UTR_3`: 5' or 3' untranslated region
- `PROMOTER`: Within 2 kb upstream of TSS
- `INTRONIC`: Intronic (not splice site)
- `GENIC_OTHER`: Within gene boundaries (other)
- `OTHER`: Intergenic

## Data Sources

- **Model:** InstaDeepAI/NTv3_650M_post (HuggingFace)
- **Annotation:** MANE Select v1.3 (RefSeq + Ensembl)
- **BigWig Metadata:** ENCODE v3, GTEx, FANTOM5, GEO, CATLAS
- **Reference Genome:** GRCh38/hg38 (local `hg38.fa` preferred; UCSC API used as fallback when local sequence is unavailable)

## Species Support

- **Human:** full app support, including BigWig context tracks and MANE transcript annotation
- **Non-human animals and plants:** BED outputs and sequence-model metrics via Ensembl / Ensembl Plants sequence retrieval
- **Important:** MAGI baseline z-scores are derived from human variants, so `Global_z_sum_log` is less directly comparable for non-human species

## Runtime Behavior

- A Hugging Face token may be required to access the gated NTv3 model weights.
- No external LLM calls are made.
- No secret-management flow is used by the Gradio app.
- Ordinary external network access can still occur when downloading model assets from Hugging Face or when falling back to the UCSC sequence API.

## Limitations

1. Predictions are computational and require experimental validation
2. The app uses a 32 kb local sequence window and does not model long-range chromatin effects
3. No phasing information is used
4. Human uses GRCh38/hg38; non-human coordinates must match the selected species assembly available through Ensembl
5. Batch processing limited to 10 variants
6. BigWig context tracks and MANE transcript annotation are currently human-only in this app
7. Region Track View zoom is limited by the available NTv3 track-profile span for the current prediction

## Citation

If you use MAGI in your research, cite the MAGI manuscript. If you rely on the underlying foundation model, also cite NTv3.
## Links

- 📄 **Paper:** *MAGI: Mechanistic Consequences of Genetic Variants via Genomic Foundation Models* (Ofer, Zok & Linial — preprint forthcoming)
- 💻 **GitHub:** https://github.com/ddofer/magi<!-- 
## License

MIT License - see LICENSE file for details -->

## Acknowledgments

- InstaDeepAI and collaborators for developing the NTv3 model
- ENCODE, GTEx, FANTOM5 consortia for functional genomics data
- NCBI RefSeq and Ensembl for transcript annotations

---

*Developed for research. Feedback and contributions welcome!*
