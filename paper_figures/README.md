# MAGI — NAR revision: print-quality figure panels

Replacement artwork for the panels flagged as blurry in
`MAGI-NAR Fig1-8-23072026.pptx` / `NAR-MAGI VER8-Formated-13072026.docx`.

> ### ▶ To regenerate on another machine, follow **[RUNBOOK.md](RUNBOOK.md)**.
> It is a copy-paste sequence with the expected output of every step.

Everything here is **true vector** (PDF + SVG, verified 0 embedded rasters) with a
600 dpi PNG alongside for tools that cannot place vector art. Fonts are Arial,
subset-embedded as TrueType (Type 42) — *not* Type 3, which journal production systems
routinely reject.

> **Source projects.** The manuscript is MAGI, whose code lives in
> `/mnt/d/Research/NT_genomics/MAGI/` (not `OpenTargetsTransfer/OTRec/gradio/`, which is
> the unrelated Open Targets recommender). The web-app panels were regenerated from
> `MAGI/gradio_app/`.

---

## What is here

| Panel | Content matches published? | Files |
|---|---|---|
| **Fig 5d** MAGI app region tracks | ✅ **yes** — same 10 tracks, Δ agree to ≤0.012 | `fig5/mat1a_gly336arg_region_tracks_*` |
| **Fig 6c** MAGI app region tracks | ✅ **yes** — same 10 tracks, Δ agree to ≤0.008 | `fig6/alox15b_rs9895916_region_tracks_*` |
| **Fig 5b** MAT1A transcript models | ✅ yes — authentic Ensembl vector export | `fig5/fig5b_*` |
| **Fig 5c** MAT1A transcript alignment | ✅ yes — recomputed from Ensembl REST | `fig5/fig5c_*` |
| **Fig 7a** GLA transcript schematic | ✅ yes — 16 transcripts, exon edges within 0.3% | `fig7/fig7a_*` |
| **Fig 7b** GLA protein MSA | ✅ yes — Clustal Omega reproduces it column-for-column | `fig7/fig7b_*` |
| **Fig 8a** UCSC track view | ✅ yes — same window and 20 tracks, ≤3 bp drift | `fig8/fig8a_*` |
| **Fig 2d** per-element disruption | ❌ **no — different model, numbers differ** | `fig2/fig2d_*` |
| **Fig 2e** benign vs pathogenic | ❌ **no — different model, numbers differ** | `fig2/fig2e_*` |

**Figures 2d and 2e are the exception and need a decision from you.** The published
versions were produced with `NTv3_100M_post` at 2 kb context from a results file that no
longer exists; everything else in the paper uses `NTv3_650M_post` at 32–64 kb. What is in
`fig2/` is the same analysis code on the current 650M data, so bar heights differ by a
median factor of ~1.4×. See `FINDINGS_FOR_AUTHORS.md` §1 and `fig2/fig2de_VALIDATION.md`.

Read next:

- **`FINDINGS_FOR_AUTHORS.md`** — nine things found along the way that need your
  decision, including two factual errors in the current figures (Fig 5b "385 aa" should
  be 395 aa; Fig 7a "363 aa" should be 362 aa) and a typo ("Trsnscript strsnd").
- `_validation/` — current vs regenerated shown side by side at matched on-page size,
  plus the numeric write-ups.
- `RESOLUTION_AUDIT.md` — measured effective DPI of every figure in the manuscript.
- `_reference_current/` — the panels as they exist in the paper today, extracted from the
  PPTX at native resolution. Reference only; not deliverables.

---

## Which file do I place?

Each app/plot panel comes in three sizes. They contain the same data.

- **`*_faithful.*`** — identical content *and proportions* to the published panel.
  Use if you are keeping the current layout and will enlarge the panel.
- **`*_print90mm.*`** — re-laid out to be legible at **90 mm** wide (≈5 pt type).
- **`*_print180mm.*`** — re-laid out to be legible at **180 mm**, full text width (≈7.5 pt).

The distinction matters because Fig 5d/6c are not merely low-resolution: the app draws
them on a 14 in canvas with 8.5 pt labels, so at their current ~2.3 in placement the type
lands at roughly **1.4 pt**. A sharper raster would still be unreadable. The print-sized
variants are rendered *at final physical size* with type set in absolute points.

Externally sourced panels (5b, 5c, 7a, 7b, 8a) come as a single vector file each; scale
to taste.

---

## Read this before dropping panels in

**1. Manual PowerPoint overlays are not included.** Arrows, red boxes and dashed lines,
callout labels ("SNP: rs118204006", "Pr: H488Q", exon numbers, `RI`/`NMD` tags, the
"Representative UCSC Tracks" side label) are native PPT objects layered over the images,
not part of the underlying artwork. They must be re-applied. Each panel's
`*_PROVENANCE.md` lists exactly which overlays it expects and, where relevant, the
coordinate the marker should sit at.

**2. Figures 5 and 6 are flattened single images in the deck.** Slide 6 is one TIFF and
slide 7 is another, so panel d / panel c cannot simply be swapped — those figures have to
be re-composited from their parts. Panels not covered by this work (5a, 6a, 6b, 6d) are
only available at their current resolution.

**3. Composite from the PPTX, not the DOCX.** The DOCX stores Figures 7 and 8 as
flattened composites at lower resolution than the layered PPTX versions.

**4. The whole figure set is below NAR's 300 dpi floor**, not just the flagged panels —
see `RESOLUTION_AUDIT.md`. Figure 6 is the worst at 193 dpi, and its panels a, b and d
are not addressed here.

---

## Reproducing

See **[RUNBOOK.md](RUNBOOK.md)** for the full procedure with expected output. Short form,
from this directory:

```bash
# Fig 5d + 6c — MAGI web-app region tracks (needs GPU; ~40 s per variant)
cd ../gradio_app && python ../paper_figures/regen_magi_app_panels.py --device cuda; cd ../paper_figures

# Fig 2d + 2e — per-element disruption (needs ../clinvar_new_deltas.parquet)
python regen_fig2_de_panels.py

# confirm everything is genuinely vector, no Type 3 fonts
bash verify_outputs.sh
```

Paths are derived from the script location, so nothing needs editing; set `MAGI_ROOT` only
if the data lives outside the repo.

Fig 5b/5c, 7a/7b and 8a are **not** regenerated by these scripts — they are finished vector
files committed here, and some are no longer rebuildable (the Ensembl 115 archive returns
HTTP 403 to automated clients). Their generators and exact source URLs are preserved in
`fig5/scripts_fig5bc/`, `fig7/make_fig7*.py` and `fig8/fig8a_hgTracks_URL.txt`.

**Two traps if you re-run the app panels yourself:** the Gradio app now defaults to a
16 kb sequence window (the paper used 32 kb) and a 128 bp zoom (the paper used a 1000 bp
radius). Both change the output materially — the window change alters which CAGE tracks
appear. `regen_magi_app_panels.py` pins the paper's values; see
`_validation/VALIDATION_5d_6c.md`.
