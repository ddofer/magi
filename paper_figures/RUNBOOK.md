# RUNBOOK — regenerate the NAR figure panels on the remote machine

Copy-paste, top to bottom. Two commands do the actual work; everything else is
checking. Should take about 5 minutes plus model load.

**What gets regenerated here:** Fig 5d, Fig 6c (need GPU + NTv3-650M) and Fig 2d, Fig 2e
(need `clinvar_new_deltas.parquet`).

**What does NOT get regenerated, and must not be deleted:** Fig 5b, 5c, 7a, 7b, 8a are
already-finished vector files committed to this repo. They were built from Ensembl /
UniProt / EBI / UCSC and **some can no longer be rebuilt** — the Ensembl 115 archive now
returns HTTP 403 to automated clients, and the UCSC panel depends on a specific cart
configuration. Leave `fig5/fig5b*`, `fig5/fig5c*`, `fig7/`, `fig8/` alone.

---

## Step 0 — get the code

```bash
git clone https://github.com/ddofer/magi.git      # or: cd <existing clone> && git pull
cd magi/paper_figures
```

Everything below is run from `magi/paper_figures`. The scripts locate the repo from
their own path, so there is nothing to edit. If your MAGI data lives somewhere other
than the repo root, set `export MAGI_ROOT=/path/to/magi` first.

## Step 1 — check the prerequisites

```bash
# poppler-utils is needed for the verification step
command -v pdffonts >/dev/null && echo "poppler OK" || echo "MISSING: sudo apt-get install -y poppler-utils"

# these must all exist
ls -la ../gradio_app/data/hg38.2bit ../gradio_app/data/functional_tracks_metadata_human.csv
ls -la ../clinvar_new_deltas.parquet
ls -la ../analyses/paper_case_studies/mat1a_gly336arg/ranked_tracks.csv
ls -la ../analyses/paper_case_studies/alox15b_rs9895916/ranked_tracks.csv
```

`hg38.2bit` (835 MB) or `hg38.fa` + `.fai` must be present — the script needs reference
sequence. The two `ranked_tracks.csv` files are the published rankings that Fig 5d/6c are
validated against; without them the script still runs but cannot self-validate.

Pick your python (must have torch+CUDA, transformers, pandas, scipy, matplotlib, pyfaidx,
py2bit, huggingface_hub):

```bash
PY=python                      # <-- set this to the interpreter/venv you use
$PY - <<'EOF'
import torch, transformers, pandas, scipy, matplotlib, pyfaidx, py2bit
print("torch", torch.__version__, "| cuda:", torch.cuda.is_available())
print("transformers", transformers.__version__, "| matplotlib", matplotlib.__version__)
EOF
```

`cuda: True` is required for a sensible runtime. NTv3-650M is gated on HuggingFace — make
sure you are logged in (`huggingface-cli login`) or `HF_TOKEN` is set.

## Step 2 — regenerate Fig 5d and Fig 6c (MAGI web-app region tracks)

```bash
cd ../gradio_app && $PY ../paper_figures/regen_magi_app_panels.py --device cuda ; cd ../paper_figures
```

Run it from `gradio_app/` — the app modules resolve some data paths relative to the
working directory. ~40 s per variant after the model loads.

**The script pins `NTV3_CONTEXT_LEN=32768` itself.** The app's own default is now 16 kb
(commit `fe87c8f`), which shifts the BED deltas ~15% and changes *which CAGE tracks appear
in the panel*. Do not override it.

### What a correct run prints

```
Vector output: Type-42 fonts, live SVG text, family=Arial
...
  top-10 track set matches the published panel; pinned to the published row order (bf16 ties)
  track set identical: True
  max |Δ difference| over the 10 rendered tracks: 0.0117
  => MATCHES published panel
```

**`=> MATCHES published panel` must appear twice** (once for MAT1A, once for ALOX15B).
Expected `max |Δ difference|`: **0.0117** for Fig 5d, **0.0078** for Fig 6c. Small
variation is fine — NTv3 runs in bf16 and the grid spacing near Δ≈0.5 is 0.0039, so
differences of one or two steps are the floor of reproducibility.

**If you see `=> DIFFERS from published panel`, stop.** It means the track *set* changed,
which is a real difference, not rounding. Read `fig5/mat1a_gly336arg_validation.txt` (and
the fig6 equivalent) to see which track moved.

## Step 3 — regenerate Fig 2d and Fig 2e

```bash
$PY regen_fig2_de_panels.py
```

### What a correct run prints

```
SOURCE: clinvar_new_deltas.parquet
        NTv3-650M @ 64 kb context - matches the rest of the paper.
        Bar heights will NOT match the printed Fig 2d/2e (median ~1.4x larger)...
  dropped 4 inference-failure row(s) (all 21 D_BED_* NaN)
  40,972 variants (of 40,976) | labels: {'Pathogenic': 22252, 'Benign': 18720}
  N pathogenic=22,252  N benign=18,720  Bonferroni alpha=0.00238
  significant elements: 21/21
```

Confirm `SOURCE:` says **`clinvar_new_deltas.parquet`** and **NTv3-650M**. If your machine
also has `clinvar_full_deltas.parquet`, the script prints a NOTE saying it found it but is
*not* using it — that is correct and intended for the 650M decision.

Expect **21/21 significant** and 40,972 variants. The 4 dropped rows (chr9:14889,
chr9:15126, chr9:15978, chr9:16020) have all-NaN BED deltas from a failed inference batch;
leaving them in makes `scipy.ttest_ind` return NaN for every element and silently yields
0/21 significant.

> Fig 2d/2e will **not** match the printed figure numerically — the published panels came
> from NTv3-100M at 2 kb context via `clinvar_full_deltas.parquet`, which is missing here.
> Bars are a median ~1.4× taller and the element order shifts from position 5 onward.
> **Any Results text quoting Fig 2d/2e values must be re-checked.** Details:
> `fig2/fig2de_VALIDATION.md`.

## Step 4 — verify everything is print-ready

```bash
bash verify_outputs.sh
```

Must end with:

```
PASS - every PDF is pure vector with embedded non-Type-3 fonts; no SVG embeds bitmaps.
```

This is the check that matters. It confirms no PDF is secretly a screenshot, and that no
Type 3 fonts crept back in — matplotlib emits Type 3 by default and OUP production rejects
them.

## Step 5 — (optional) 600 dpi PNGs

Only needed for tools that cannot place vector art. The `.pdf`/`.svg` are the real
deliverables. The two regeneration scripts already write PNGs; for the pre-built panels
(5b, 5c, 7, 8) that are committed as vector only:

```bash
for f in fig5/fig5b*.pdf fig5/fig5c_MAT1A_transcript_alignment.pdf fig7/*.pdf fig8/*.pdf; do
  [ -e "$f" ] && pdftoppm -r 600 -png -singlefile "$f" "${f%.pdf}"
done
```

---

## Which file goes into the figure

Each regenerated panel comes in three sizes with identical data:

| suffix | use when |
|---|---|
| `_faithful` | keeping the published proportions **and** enlarging the panel |
| `_print90mm` | the panel stays ~90 mm wide — type set to land at ~5 pt on the page |
| `_print180mm` | the panel gets the full 180 mm text width — type at ~7.5 pt |

This matters because Fig 5d/6c are not merely low-resolution. The app draws them on a
14-inch canvas with 8.5 pt labels; at their current ~2.3-inch placement that type lands at
about **1.4 pt**. A sharper raster would still be unreadable — the layout has to change.
Use `_print90mm` if the panel keeps its current footprint.

## Before you drop panels in

1. **PowerPoint overlays are not in these files.** Arrows, red boxes and dashed lines,
   "SNP: rs118204006", "Pr: H488Q", exon numbers, `RI`/`NMD` tags, the aa-length labels,
   the "Representative UCSC Tracks" side label — all native PPT objects layered on top.
   Re-apply them. Each `fig*/​*_PROVENANCE.md` lists which overlays its panel expects, and
   gives marker positions as fractions of panel width.
2. **Figures 5 and 6 are single flattened TIFFs in the deck** (slide 6, slide 7), so a
   panel cannot simply be swapped — those figures must be re-composited. Panels 5a, 6a,
   6b and 6d have no regenerable source and stay at their current resolution.
3. **Composite from the PPTX, not the DOCX** — the DOCX stores Figures 7 and 8 as
   lower-resolution flattened copies.

## Troubleshooting

| symptom | cause / fix |
|---|---|
| `Unable to download ... NTv3_650M_post` | model is gated: `huggingface-cli login`, or `export HF_TOKEN=...` |
| `Sequence fetch failed` | `gradio_app/data/hg38.2bit` (or `hg38.fa` + `.fai`) missing |
| `=> DIFFERS from published panel` | usually `NTV3_CONTEXT_LEN` overridden in the environment — unset it; the script sets 32768 |
| `delta parquet not found` | pass it: `--deltas /path/to/clinvar_new_deltas.parquet` |
| `0/21 significant` in step 3 | the NaN-row drop was disabled; `DROP_ALL_NAN_BED` must stay `True` |
| CUDA OOM | `--device cpu` (slow but correct), or free GPU memory |
| `verify_outputs.sh`: `MISSING: pdffonts` | `sudo apt-get install -y poppler-utils` |
| fonts come out DejaVu not Arial | cosmetic only; Arial not installed. Type 42 embedding still applies, so output is still valid |
