# Figure 2d / 2e — regeneration validation

Generated 2026-07-25 · `regen_fig2_de_panels.py` · statistics verbatim from
`MAGI/analyses/genomic_featImp_experiment_v3.ipynb` **cell 24**.

## TL;DR — negative result on "find the right subset"

**There is no subset.** The published panels were computed from a results file
that no longer exists on this machine: **`MAGI/clinvar_full_deltas.parquet`**,
produced by **`InstaDeepAI/NTv3_100M_post` at `CONTEXT_LEN = 2048`**. The file
present today, `clinvar_new_deltas.parquet`, holds the **same 40,976 input
variants** re-scored by **`InstaDeepAI/NTv3_650M_post` at
`CONTEXT_LEN = 65536`**. The two runs differ in the model, not in the rows, so
no row filter over `clinvar_new_deltas.parquet` can reproduce the published bar
heights.

The regenerated panels are therefore **structurally identical** (same code, same
population, same ordering rule, same colours, all 21 bars green) but the bar
**heights are ~1.4x larger** on average because they come from the larger model.

## Which subset was identified, and the evidence

### 1. The published run used *all* 40,976 labelled variants — no filtering
The notebook's own stored output for cell 24 reads:

```
   N Pathogenic: 22252
   N Benign: 18724
   Bonferroni-corrected α: 0.0024
   Significant elements: 20/21
```

22,252 + 18,724 = 40,976 = every labelled row. Cell 8's stored output confirms
the load: `Loaded 40976 variants / Columns: 14771 / Labels: {'Pathogenic':
22252, 'Benign': 18724}`. Between cell 8 (`final_df = pd.read_parquet(
OUTPUT_RESULTS_FILE)`) and cell 24 the only assignment to `df` is cell 11's
`df = final_df.copy()`, which merely appends an `Impact_Score` column. So the
manuscript's Fig 2b filtering cascade (27,853 → 11,240 ≥2-gold-star → 9,786 with
rationale) is **not** applied to panels d/e. Hypothesis 1 is refuted.

### 2. The published run used a different results file
Cell 5 of the notebook:

```python
MODEL_NAME          = "InstaDeepAI/NTv3_100M_post"
CONTEXT_LEN         = 2048
OUTPUT_RESULTS_FILE = "clinvar_full_deltas.parquet"
```

`clinvar_full_deltas.parquet` is absent from the whole of `/mnt/d/Research`
(searched by name and by size; the only ≥1 GB delta table is
`clinvar_new_deltas.parquet`). Meanwhile `MAGI/inference.py`, whose documented
default output *is* `clinvar_new_deltas.parquet`, sets

```python
CONTEXT_LEN = 64 * 1024
MODEL_NAME  = "InstaDeepAI/NTv3_650M_post"  # InstaDeepAI/NTv3_100M_post
```

Four independent fingerprints separate the two files:

| property | published `clinvar_full_deltas.parquet` (from notebook stdout) | on-disk `clinvar_new_deltas.parquet` (measured) |
|---|---|---|
| rows | 40,976 | 40,976 |
| label split | 22,252 P / 18,724 B | 22,252 P / 18,724 B |
| total columns | 14,771 | 7,270 |
| `D_BW_*` tracks | 7,362 (222 `kai*`, 5,609 `ENCSR*`, 1,276 `CNhs*`, 84 `GTEX*`) | 3,602 (0 `kai*`, 2,240 `ENCSR*`, 1,276 `CNhs*`, 84 `GTEX*`) |
| `D_BED_*` min / max | −0.961818 / +0.946396 | −0.774623 / +0.714373 |
| `D_BED_*` mean / std | −0.000968 / 0.032843 | −0.000535 / 0.029819 |
| `D_BED_*` mean\|Δ\| / median\|Δ\| | 0.010791 / 0.002336 | 0.013053 / 0.006385 |
| all-NaN `D_BED_*` rows | 0 (all Welch p-values finite) | 4 |

The decisive one is the range: **a subset can only shrink a range, never widen
it.** The published data reach ±0.96 while every row on disk lies inside ±0.78,
so the published rows are not a subset of `clinvar_new_deltas.parquet`.

### 3. Other candidate sources were checked and excluded
* `indel_size` in `clinvar_new_deltas.parquet` is **0 for all 40,976 rows** — the
  file is already SNP-only, so "SNPs vs SNPs+indels" (hypothesis 2) cannot be the
  difference.
* `clinvar_input.parquet` (40,976 × 13) has no delta columns and no review
  status.
* `data/stav_data/annotated_snps_signaled.parquet` (9,786 rows, with
  `gold_stars`) and `snps_annotated_test.parquet`,
  `annotated_snps_signaled_impacted.parquet` all carry `D_BED_*` values that are
  **bit-identical** to `clinvar_new_deltas.parquet` on the shared
  chrom/pos/ref/alt key (max abs difference = 0.0). They are downstream of the
  650M run, so they contain no archived copy of the published numbers.
* `clinvar_indel_deltas.parquet` is the indel run (13,064 rows) — wrong
  population.
* `analyses/LLM Judgement results/snp_evaluation_results_v2.csv` holds LLM-judge
  output keyed by `#VariationID`; no delta columns.
* No `*deltas*` file anywhere under `/mnt/d/Research/NT_genomics/` or
  `/mnt/d/Research/OpenTargetsTransfer/` matches the published fingerprint.

### 4. What *did* get corrected
The 4 rows whose 21 `D_BED_*` values are all NaN (chr9:14889 A>G, chr9:15126
G>C, chr9:15978 C>T, chr9:16020 A>G — all Benign; a failed inference batch) are
now dropped, since `scipy.stats.ttest_ind` propagates NaN and was returning NaN
p-values for every element (0/21 "significant"). With them removed the analysis
population is 40,972 (22,252 Pathogenic / 18,720 Benign) and **21/21 elements
are Bonferroni-significant**, matching the all-green published panel e. The
published run had no such rows (its N Benign was the full 18,724).

## Published vs regenerated per-element means (all 21 elements)

Published values: 10 are the **exact** figures printed in the notebook's stored
cell-24 stdout (`results_df.head(10)`); the remaining 11 were recovered by
pixel-digitising the notebook's own stored panel-d PNG (identical to
`analyses/analyses/figures/fig-per_element_disruption_by_class.png`). The
digitiser was calibrated on the 7 y-gridlines and validated against the 10 exact
values: **max residual 4·10⁻⁵**, RMS 2·10⁻⁵. (The `intron` Pathogenic bar top is
hidden behind the legend box and was recovered from the legend's alpha-blended
pixels.)

Regenerated values are from `clinvar_new_deltas.parquet`, 40,972 rows.

| # | element | pub Path | regen Path | ΔPath | pub Benign | regen Benign | ΔBenign | pub diff (B−P) | regen diff | source of pub value |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `polyA_signal` | 0.00637 | 0.01968 | +0.01331 | 0.00859 | 0.02291 | +0.01432 | +0.00222 | +0.00323 | notebook stdout (exact) |
| 2 | `protein_coding_gene` | 0.00300 | 0.00799 | +0.00499 | 0.00408 | 0.00742 | +0.00334 | +0.00108 | -0.00057 | notebook stdout (exact) |
| 3 | `CTCF-bound` | 0.00313 | 0.00562 | +0.00249 | 0.00323 | 0.00498 | +0.00175 | +0.00009 | -0.00064 | notebook stdout (exact) |
| 4 | `3UTR+` | 0.00405 | 0.00699 | +0.00294 | 0.00362 | 0.00587 | +0.00225 | -0.00043 | -0.00112 | notebook stdout (exact) |
| 5 | `3UTR-` | 0.00435 | 0.00676 | +0.00241 | 0.00387 | 0.00483 | +0.00096 | -0.00048 | -0.00193 | notebook stdout (exact) |
| 6 | `lncRNA` | 0.00639 | 0.00816 | +0.00178 | 0.00570 | 0.00701 | +0.00131 | -0.00068 | -0.00115 | notebook stdout (exact) |
| 7 | `stop_codon` | 0.01588 | 0.02287 | +0.00698 | 0.01505 | 0.02027 | +0.00522 | -0.00083 | -0.00260 | notebook stdout (exact) |
| 8 | `promoter_Tissue_invariant` | 0.00285 | 0.00479 | +0.00194 | 0.00196 | 0.00363 | +0.00166 | -0.00089 | -0.00116 | notebook stdout (exact) |
| 9 | `promoter_Tissue_specific` | 0.00567 | 0.01010 | +0.00443 | 0.00424 | 0.00807 | +0.00384 | -0.00143 | -0.00202 | notebook stdout (exact) |
| 10 | `enhancer_Tissue_invariant` | 0.00508 | 0.00952 | +0.00445 | 0.00357 | 0.00640 | +0.00283 | -0.00151 | -0.00313 | notebook stdout (exact) |
| 11 | `splice_acceptor` | 0.00502 | 0.01198 | +0.00696 | 0.00313 | 0.00954 | +0.00641 | -0.00189 | -0.00244 | digitised ±4e-5 |
| 12 | `5UTR+` | 0.00889 | 0.01157 | +0.00268 | 0.00542 | 0.00832 | +0.00290 | -0.00347 | -0.00326 | digitised ±4e-5 |
| 13 | `5UTR-` | 0.00929 | 0.01074 | +0.00145 | 0.00573 | 0.00798 | +0.00225 | -0.00356 | -0.00276 | digitised ±4e-5 |
| 14 | `enhancer_Tissue_specific` | 0.01924 | 0.02031 | +0.00107 | 0.01474 | 0.01660 | +0.00186 | -0.00450 | -0.00371 | digitised ±4e-5 |
| 15 | `exon` | 0.01908 | 0.01706 | -0.00202 | 0.01150 | 0.01054 | -0.00096 | -0.00758 | -0.00653 | digitised ±4e-5 |
| 16 | `skipped_exon` | 0.02413 | 0.01862 | -0.00551 | 0.01560 | 0.01229 | -0.00331 | -0.00853 | -0.00633 | digitised ±4e-5 |
| 17 | `splice_donor` | 0.01489 | 0.02105 | +0.00616 | 0.00526 | 0.01432 | +0.00906 | -0.00963 | -0.00673 | digitised ±4e-5 |
| 18 | `always_on_exon` | 0.02729 | 0.01991 | -0.00738 | 0.01592 | 0.01234 | -0.00358 | -0.01137 | -0.00757 | digitised ±4e-5 |
| 19 | `ORF` | 0.02871 | 0.02424 | -0.00447 | 0.01616 | 0.01491 | -0.00125 | -0.01255 | -0.00932 | digitised ±4e-5 |
| 20 | `start_codon` | 0.02713 | 0.02870 | +0.00157 | 0.01458 | 0.02073 | +0.00615 | -0.01255 | -0.00798 | digitised ±4e-5 |
| 21 | `intron` | 0.02879 | 0.02309 | -0.00570 | 0.01434 | 0.01282 | -0.00152 | -0.01445 | -0.01027 | digitised ±4e-5 |

```
Mean signed shift  Pathogenic +0.00193, Benign +0.00264
Mean |shift|       Pathogenic 0.00432, Benign 0.00365
Max  |shift|       Pathogenic 0.01331, Benign 0.01432
Ratio regen/pub    Pathogenic median 1.41x, Benign median 1.42x
```

The regenerated bars are systematically larger for most elements (the 650M model
at 64 kb context produces larger BED-probability perturbations) but *smaller*
for the strongest coding elements (`always_on_exon` −0.0074, `intron` −0.0057,
`skipped_exon` −0.0055, `ORF` −0.0045, `exon` −0.0020), i.e.
the change is element-specific, not a global rescale — exactly what a different
model checkpoint produces and what no row filter could produce.

## Element ordering comparison

Ordering rule (unchanged): `sort_values("Difference_Benign_minus_Path",
ascending=False)`.

| position | published | regenerated | shift |
|---|---|---|---|
| 1 | `polyA_signal` | `polyA_signal` | — |
| 2 | `protein_coding_gene` | `protein_coding_gene` | — |
| 3 | `CTCF-bound` | `CTCF-bound` | — |
| 4 | `3UTR+` | `3UTR+` | — |
| 5 | `3UTR-` | `lncRNA` | ≠ |
| 6 | `lncRNA` | `promoter_Tissue_invariant` | ≠ |
| 7 | `stop_codon` | `3UTR-` | ≠ |
| 8 | `promoter_Tissue_invariant` | `promoter_Tissue_specific` | ≠ |
| 9 | `promoter_Tissue_specific` | `splice_acceptor` | ≠ |
| 10 | `enhancer_Tissue_invariant` | `stop_codon` | ≠ |
| 11 | `splice_acceptor` | `5UTR-` | ≠ |
| 12 | `5UTR+` | `enhancer_Tissue_invariant` | ≠ |
| 13 | `5UTR-` | `5UTR+` | ≠ |
| 14 | `enhancer_Tissue_specific` | `enhancer_Tissue_specific` | — |
| 15 | `exon` | `skipped_exon` | ≠ |
| 16 | `skipped_exon` | `exon` | ≠ |
| 17 | `splice_donor` | `splice_donor` | — |
| 18 | `always_on_exon` | `always_on_exon` | — |
| 19 | `ORF` | `start_codon` | ≠ |
| 20 | `start_codon` | `ORF` | ≠ |
| 21 | `intron` | `intron` | — |

Kendall τ between the published and regenerated orderings = **0.905**; total
rank displacement = 20 positions; 8 of 21 elements sit in a different slot. The
head (`polyA_signal`, `protein_coding_gene`, `CTCF-bound`, `3UTR+`) and the tail
(`always_on_exon`, `ORF`/`start_codon`, `intron`) are preserved; the middle
reshuffles.

For completeness, the same statistics run on the quality-filtered candidates
give *worse* agreement with the published ordering, which is further evidence
against hypothesis 1:

| candidate subset | n | P / B | significant | Kendall τ vs published order |
|---|---|---|---|---|
| **all labelled rows, 4 NaN rows dropped (used)** | 40,972 | 22,252 / 18,720 | **21/21** | **0.905** |
| all labelled rows, NaN kept | 40,976 | 22,252 / 18,724 | 0/21 (NaN p-values) | 0.905 |
| ≥2 gold stars **and** ClinVar rationale (Fig 2b's 9,786) | 9,786 | 6,932 / 2,854 | 21/21 | 0.819 |
| ≥3 gold stars | 1,133 | 1,060 / 73 | 5/21 | 0.562 |

## Significance count

* Published: **20/21** significant at Bonferroni α = 0.05/21 = 0.00238.
  `CTCF-bound` is the single non-significant element (p = 0.4616, difference
  +0.000093). Pixel-scanning the notebook's stored panel-e PNG confirms this:
  only **20** green (`#008000` @ α 0.75 → RGB 64,160,64) bar groups are present,
  and the missing one is at the x-slot of the 3rd category (`CTCF-bound`, pixels
  x = 272–315), where the fill is light-grey (RGB ≈ 222,222,222) and only ~1 px
  tall. So "all 21 bars green" in the brief is a mis-read of a sub-pixel
  light-grey bar; the true published count is 20/21.
* Regenerated: **21/21** significant (`CTCF-bound` p = 7.5e-11 in the 650M run).
  Panel e is genuinely all-green.

## Axis ranges

| | published | regenerated |
|---|---|---|
| panel d y-range | 0 → ~0.0305 | 0 → ~0.0305 |
| panel e y-range | +0.0025 → −0.0150 | +0.0040 → −0.0110 |

Panel d happens to land on the same range because the tallest published bar
(`intron`, 0.0288) and the tallest regenerated bar (`start_codon`, 0.0287) are
nearly equal. Panel e is compressed because the largest regenerated
Benign−Pathogenic gap (`intron`, −0.0103) is smaller than the published one
(−0.0145), while `polyA_signal` is larger (+0.0032 vs +0.0022).

## Honest statement of residual differences

1. **The bar heights do not match the published figure and cannot be made to
   match from any file currently on disk.** Mean |shift| is 0.0043 (Pathogenic)
   and 0.0037 (Benign) on values of order 0.005–0.030 — i.e. tens of percent.
   The regenerated panels are a *re-analysis with a newer model*, not a
   re-rendering of the published one.
2. **The element ordering differs in 8 of 21 slots** (τ = 0.905). The
   qualitative story is unchanged: `polyA_signal` and `protein_coding_gene` are
   the only elements more disrupted in Benign variants; coding/splicing elements
   (`intron`, `ORF`, `start_codon`, `always_on_exon`) are the most
   Pathogenic-skewed.
3. **Significance count differs by one** (21/21 vs 20/21), because `CTCF-bound`
   becomes significant under the 650M model.
4. Two published values used in the table (11 of 42 numbers) are digitised
   rather than exact, with a validated accuracy of ±4·10⁻⁵ — negligible relative
   to the discrepancies above, but they are estimates, not authoritative
   figures.
5. Everything that is *not* data-dependent is faithful: the statistics code is
   cell 24 verbatim, the figure size (14 × 6 in), colours (`#d62728` /
   `#1f77b4`, α = 0.75), 95 % CI construction, Bonferroni threshold, sort key,
   titles, axis labels and tick rotation are all unchanged.

## What is needed to reproduce the published panels exactly

Exactly one of:

* **`clinvar_full_deltas.parquet`** — the 40,976 × 14,771 NTv3-100M / 2 kb
  results table (or just its 21 `D_BED_*` columns + `label`); drop it into
  `/mnt/d/Research/NT_genomics/MAGI/` and re-run `regen_fig2_de_panels.py`
  unchanged — the script already prefers that filename; **or**
* **`clinvar_full_deltas_with_snp_ids.parquet`**, referenced at
  `analyses/stav-analysis.ipynb` cell 5, which is the same table plus ClinVar
  IDs and is likewise absent; **or**
* the saved `results_df` from cell 24 (all 21 rows with `Pathogenic_Mean`,
  `Pathogenic_SEM`, `Benign_Mean`, `Benign_SEM`, `P_Value`) as a CSV — enough to
  redraw both panels exactly without the raw deltas; **or**
* a re-run of `inference.py` with `MODEL_NAME = "InstaDeepAI/NTv3_100M_post"` and
  `CONTEXT_LEN = 2048` over `clinvar_input.parquet` (bit-exact reproduction is
  not guaranteed, but the model/context would then match).

If the authors instead wish the manuscript to carry the **650M** numbers, the
files in this directory are ready to use as-is; the panel-d/e captions and any
in-text numbers quoted from them would need updating.

## Output inventory (this directory)

| file | notes |
|---|---|
| `fig2d_per_element_by_variant_class_{faithful,print90mm,print180mm}.{pdf,svg,png}` | 9 files |
| `fig2e_per_element_benign_vs_pathogenic_{faithful,print90mm,print180mm}.{pdf,svg,png}` | 9 files |
| `fig2de_per_element_statistics.csv` | all 21 elements × 18 statistics |
| `fig2de_provenance.json` | source file, model, filters, counts, ordering |

All 6 PDFs verified true vector: `pdfimages -list` reports **0 embedded
rasters** for each. No SVG contains a `data:image` payload. All PNGs written at
**600 dpi** (faithful 8338 × 3533 px; print180mm 4190 × 1904; print90mm
2029 × 1193).
