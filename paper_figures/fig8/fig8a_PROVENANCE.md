# Figure 8a — provenance and regeneration recipe

Print-quality replacement for the low-resolution UCSC Genome Browser screenshot
used in MAGI Fig. 8a (NAR submission).

Generated: **2026-07-25** · Source: **UCSC Genome Browser v500**, `genome.ucsc.edu`
(rev. 2 — widened label column, tight page box; see §9 changelog)

---

## 1. Deliverables

| File | Type | Size | Notes |
|---|---|---|---|
| `fig8a_UCSC_hg38_chr19_rs897804.pdf` | **True vector PDF** | **540 × 344.03 pt** (7.500 × 4.778 in / 190.5 × 121.4 mm) | Server-side PostScript/PDF export from UCSC. No embedded raster images; all text is live text in an embedded Helvetica (Type1/CFF) subset. **/MediaBox and /CropBox are both set to the tight content box** `[0 58.97 540 403]`, so tools that ignore /CropBox (Ghostscript, LaTeX `\includegraphics`, Illustrator place, Word) still get the cropped panel. **Use this for the manuscript.** |
| `fig8a_UCSC_hg38_chr19_rs897804.png` | High-res raster | **4500 × 2867 px** (600 dpi at the PDF's native size) | Rasterised from the PDF above. ≈1286 dpi if placed at a 3.5 in single-column width; 600 dpi at full 7.5 in width. Replaces the old 1292 px (~200 dpi) screenshot. |
| `fig8a_hgTracks_URL.txt` | Text | — | The exact `curl`-ready UCSC configuration URL (§6). |
| `fig8a_PROVENANCE.md` | This file | — | — |

---

## 2. Genomic coordinates

| Item | Value |
|---|---|
| Assembly | **hg38 / GRCh38** (UCSC `db=hg38`) |
| Displayed range | **chr19:12,765,443–12,766,947** (1-based inclusive; **1,505 bp**) |
| Variant | **rs897804**, hg38 **chr19:12,766,150** (verified via `api.genome.ucsc.edu/search?search=rs897804&genome=hg38`) |
| Gene / protein change | **HOOK2**, p.His488Gln (exon 15 of 22) |

**Why this exact range?** It was reverse-engineered from the original low-res
panel rather than guessed. Feature edges of known genomic coordinate were
measured in `_reference_current/CURRENT_fig8a_UCSC.png` and solved against the
UCSC REST API:

* GeneHancer `GH19J012765` = chr19:12,765,626–12,766,626 → drawn at x = 301…1050 px (750 px / 1000 bp = **0.750 px/bp** in the original)
* CpG island `CpG: 22` = chr19:12,766,075–12,766,294 → drawn at x = 638…802 px
* ORegAnno `OREG0025273` starts left of the window and is clipped at x = 163 px → that is the left edge of the original's data area

Solving gives start = 12,765,442 (0-based) and a data width of 1,505 bp, i.e.
**chr19:12,765,443-12,766,947**. This is the ~1,500 nt window quoted in the
published legend and places rs897804 at 47.0 % of the data width — identical to
the original panel.

---

## 3. Tracks, in the exact top-to-bottom order of the figure

| # | UCSC track (table) name | Visibility | Label shown in figure |
|---|---|---|---|
| 0 | `ruler` | dense | Scale bar + chr19 position ruler + "hg38" |
| 1 | `mane` | pack | MANE Select Plus Clinical: Representative transcript from RefSeq & GENCODE (**HOOK2**) |
| 2 | `cCREregistry` (superTrack `cCREs`) | dense | ENCODE4 Registry of candidate Cis-Regulatory Elements (cCREs) |
| 3 | `rmsk` | dense | Repeating Elements by RepeatMasker |
| 4 | `cpgIslandExt` (superTrack `cpgIslandSuper`) | pack | CpG Islands (Islands < 300 Bases are Light Green) — item `CpG: 22` |
| 5 | `hmaSummaryUnmethylated` (composite `humanMethylationAtlasSummary`, superTrack `dnaMethylation`) | pack | All unmethylated regions — Adipocytes, Gastric-Ep, Head-Neck-Ep, Oligodend (+ 2nd Head-Neck-Ep row) |
| 6 | `oreganno` | pack | Regulatory elements from ORegAnno — OREG0025273, OREG1386284, OREG1362738, OREG1790119, OREG0677122 |
| 7–16 | `humanMethylationAtlasSignals` composite, 10 default-on subtracks, each **dense**: `neuronMerged`, `heartCardioMerged`, `endothelMerged`, `bloodTMerged`, `bloodBMerged`, `bloodGranulMerged`, `lungAlveoEpMerged`, `pancBetaMerged`, `liverHepMerged`, `colonEpMerged` | dense | Human Methylation Atlas WGBS cell type signals — Neurons Merged, Heart Cardiomyocytes Merged, Endothelial Merged, Blood T Cells Merged, Blood B Cells Merged, Blood Granulocytes Merged, Lung Alveolar Epithelium Merged, Pancreas Beta Cells Merged, Liver Hepatocytes Merged, Colon Epithelium Merged |
| 17 | `refSeqFuncElems` | pack | NCBI RefSeq Functional Elements (no items in this window — centre label only, as in the original) |
| 18 | `ReMapDensity` (composite `ReMap`; `ReMapTFs` hidden) | full | ReMap density, y-axis 1 … 1426 |
| 19 | `geneHancerRegElementsDoubleElite` (composite `geneHancer` / view `ghGeneHancer`) | pack | Enhancers and promoters from GeneHancer (Double Elite) — `GH19J012765` |

Explicitly hidden so they do not appear: `cpgIslandExtUnmasked`,
`hmaSummaryPutEnhancers`, `hmaSummaryU250`, `ReMapTFs`, `geneHancerRegElements`,
`ghInteraction` / `geneHancerInteractions*`, `ghGeneTss` / `geneHancerGenes*`,
`ghClusteredInteraction` / `geneHancerClusteredInteractions*`.

### Image parameters

| Parameter | Value | Rationale |
|---|---|---|
| `pix` | **1360** | Widest setting that still reproduces the original's **500-base ruler tick spacing**; at `pix ≥ 1380` UCSC switches to 100-base ticks. It also keeps the printed text at ~5.6 pt when the PDF is placed at its native 7.5 in width — the same physical text size as the original panel. |
| `textSize` | **14** | Matches the original. |
| `hgt.labelWidth` | **31** | Character width of the left label column. UCSC's default (20) truncates the longest Methylation-Atlas row labels; 31 is the length of the longest label present, `Lung Alveolar Epithelium Merged`. See §8. |
| `guidelines` | on | Light-blue vertical guides, as in the original. |

### Track data versions (as served on 2026-07-25)

* MANE — see `https://genome.ucsc.edu/gbdb/hg38/mane/README_versions.txt`
* ENCODE cCREs — **ENCODE Registry version 4, 2024** (includes ENCODE2, ENCODE3 and Roadmap Epigenomics)
* Human Methylation Atlas (Loyfer et al.) — **Data release version 2**
* ReMap density — `reMapDensity2022.bw` (ReMap 2022)
* GeneHancer — **January 2019 (V2: corrections to Experiment field)**
* RepeatMasker, CpG Islands, ORegAnno, NCBI RefSeq Functional Elements — current UCSC hg38 releases

---

## 4. Method that worked (vector, first choice)

UCSC's **server-side PostScript/PDF export** (`hgt.psOutput=on`) was used, driven
by a cookie-backed cart built with `curl`. This yields a genuine vector PDF
(no embedded bitmap), which is why the output is print-quality rather than a
scaled-up screenshot.

Three facts were essential and are worth recording:

1. **Left-hand track labels are a separate image in the HTML browser.** In
   drag-scroll mode `hgt_genome_*.png` contains only the data area; the side
   labels live in `hgtSide_genome_*.png`. The PDF/PS export renders a single
   complete image with the side labels included, so it avoids that problem
   entirely. (A naive `hgRenderTracks`/`hgt_genome_*.png` grab loses all the
   left labels.)
2. **Cross-group vertical ordering** (RepeatMasker sits between the cCREs and the
   CpG Islands, i.e. out of its normal "Repeats" group position) is *not*
   controllable via `<track>.priority`, which only sorts within a track group.
   It is controlled by the drag-reorder cart variables
   **`<trackName>_imgOrd=<row index>`** (see `dragReorder.setOrder` in
   `https://genome.ucsc.edu/js/utils.js`).
3. **`hgt.reset=1` is applied after the other parameters**, so it must be sent as
   its own separate request before the configuration request, not merged into it.

## 5. How to regenerate

```bash
CJ=$(mktemp)                 # cookie jar = the UCSC cart
UA="Mozilla/5.0"

# Step 1 — start from a clean cart (MUST be a separate request; see §4.3)
curl -s -L -c "$CJ" -b "$CJ" -A "$UA" \
  "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.reset=1" -o /dev/null

# Step 2 — apply position, tracks, visibilities, label width and vertical order
curl -s -L -c "$CJ" -b "$CJ" -A "$UA" "$(cat fig8a_hgTracks_URL.txt)" -o /dev/null

# Step 3 — ask UCSC for the vector PDF, then download it
curl -s -L -c "$CJ" -b "$CJ" -A "$UA" \
  "https://genome.ucsc.edu/cgi-bin/hgTracks?hgt.psOutput=on" -o ps.html
grep -o '\.\./trash/hgt/hgt_genome_[^"]*\.pdf' ps.html | head -1 \
  | sed 's|^\.\.|https://genome.ucsc.edu|' \
  | xargs -I{} curl -s -b "$CJ" -A "$UA" {} -o raw.pdf
```

`raw.pdf` is 540 × 403 pt with ~59 pt of empty band below the GeneHancer row.
Post-processing (PyMuPDF; **stays vector — nothing is rasterised to crop**):

1. **Tight page box.** The drawn content ends 334.83 pt from the top; adding the
   same 9.20 pt bottom margin the original panel had gives a page height of
   **344.03 pt**. Because PDF user space has its origin at the *bottom* left and
   the content sits at the *top* of the sheet, the tight box is
   `[0, 403−344.03, 540, 403]` = `[0 58.97 540 403]`. Set **both** `/MediaBox`
   and `/CropBox` to it (setting only `/CropBox` leaves Ghostscript, LaTeX and
   other MediaBox-based consumers showing the full 403 pt sheet with ~15 % white
   space at the bottom).
2. **Rasterise** the cropped PDF at 600 dpi → the 4500 × 2867 px PNG.

```python
import fitz
W, H, h = 540.0, 403.0, 344.03
d = fitz.open("raw.pdf"); pg = d[0]
box = f"[0 {H-h:.2f} {W:.2f} {H:.2f}]"
for k in ("MediaBox", "CropBox"):
    d.xref_set_key(pg.xref, k, box)
d.save("fig8a_UCSC_hg38_chr19_rs897804.pdf", garbage=4, deflate=True, clean=True)

d = fitz.open("fig8a_UCSC_hg38_chr19_rs897804.pdf")
d[0].get_pixmap(matrix=fitz.Matrix(600/72, 600/72),
                colorspace=fitz.csRGB, alpha=False
      ).save("fig8a_UCSC_hg38_chr19_rs897804.png")
```

Verified afterwards with MediaBox-based tools:
`pdfinfo` → `Page size: 540 x 344.03 pts`;
`gs -sDEVICE=bbox` → `%%BoundingBox: 1 0 540 344` (content fills the page).

## 6. Exact configuration URL (step 2)

Stored ready to `curl` in **`fig8a_hgTracks_URL.txt`**.

```
https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr19%3A12765443-12766947&hideTracks=1&pix=1360&textSize=14&hgt.labelWidth=31&guidelines=on&ruler=dense&mane=pack&cCREs=show&cCREregistry=dense&rmsk=dense&cpgIslandSuper=show&cpgIslandExt=pack&cpgIslandExtUnmasked=hide&dnaMethylation=show&humanMethylationAtlasSummary=pack&hmaSummaryUnmethylated=pack&hmaSummaryUnmethylated_sel=1&hmaSummaryPutEnhancers=hide&hmaSummaryPutEnhancers_sel=0&hmaSummaryU250=hide&hmaSummaryU250_sel=0&oreganno=pack&humanMethylationAtlasSignals=dense&neuronMerged=dense&heartCardioMerged=dense&endothelMerged=dense&bloodTMerged=dense&bloodBMerged=dense&bloodGranulMerged=dense&lungAlveoEpMerged=dense&pancBetaMerged=dense&liverHepMerged=dense&colonEpMerged=dense&neuronMerged_sel=1&heartCardioMerged_sel=1&endothelMerged_sel=1&bloodTMerged_sel=1&bloodBMerged_sel=1&bloodGranulMerged_sel=1&lungAlveoEpMerged_sel=1&pancBetaMerged_sel=1&liverHepMerged_sel=1&colonEpMerged_sel=1&refSeqFuncElems=pack&ReMap=full&ReMapDensity=full&ReMapDensity_sel=1&ReMapTFs=hide&ReMapTFs_sel=0&geneHancer=pack&ghGeneHancer=pack&geneHancerRegElementsDoubleElite=pack&geneHancerRegElementsDoubleElite_sel=1&geneHancerRegElements=hide&geneHancerRegElements_sel=0&ghInteraction=hide&ghGeneTss=hide&ghClusteredInteraction=hide&ruler_imgOrd=0&mane_imgOrd=1&cCREregistry_imgOrd=2&rmsk_imgOrd=3&cpgIslandExt_imgOrd=4&hmaSummaryUnmethylated_imgOrd=5&oreganno_imgOrd=6&neuronMerged_imgOrd=7&heartCardioMerged_imgOrd=8&endothelMerged_imgOrd=9&bloodTMerged_imgOrd=10&bloodBMerged_imgOrd=11&bloodGranulMerged_imgOrd=12&lungAlveoEpMerged_imgOrd=13&pancBetaMerged_imgOrd=14&liverHepMerged_imgOrd=15&colonEpMerged_imgOrd=16&refSeqFuncElems_imgOrd=17&ReMapDensity_imgOrd=18&geneHancerRegElementsDoubleElite_imgOrd=19
```

Opening that URL in a browser (after a `hgt.reset=1`) reproduces the same view
interactively; the PDF can then be obtained from **View → PDF/PS**.

---

## 7. ⚠ Manual annotations that MUST be re-applied by the authors

The files delivered here are the **clean browser image only**. Everything below
was added on top of the original screenshot in PowerPoint and is *not* present
and *not* reproducible from UCSC — it has to be re-drawn over the new asset:

1. The **red vertical line** marking the rs897804 position.
2. The green **"rs897804"** callout box (top right).
3. The **"Pr: H488Q"** box (top right).
4. The four red-outlined annotation boxes down the right-hand side:
   * "Unmethylated Cell and tissues"
   * "Open regulatory annotations"
   * "Regulatory (Experimental) density: TF, ChIp"
   * "GeneCard: Elite Promoter/Enhancer"
5. The rotated **"Representative UCSC Tracks"** label.
6. The thin **black frame** drawn around the whole panel.

### Where to put the rs897804 line — measured on the delivered files

The label column occupies the leftmost **19.27 %** of the image; the data area is
the remaining **80.73 %**. rs897804 (chr19:12,766,150) sits at **47.01 %** across
the data area, i.e. at **57.22 %** of the total image width:

| Asset | x position of the rs897804 line |
|---|---|
| `…rs897804.png` (4500 px wide) | **x ≈ 2575 px** from the left edge |
| `…rs897804.pdf` (540 pt wide) | **x ≈ 309.0 pt** from the left edge |
| any scaled placement | **x = 0.5722 × (placed width)** |

(These values changed from rev. 1 because the label column was widened; see §9.)

Note that in the original panel the right-hand annotation boxes *overlap and
obscure* the right-hand end of several tracks (the ENCODE4 cCRE and RepeatMasker
blocks and the right end of the methylation-atlas rows). Those data are fully
visible in the new asset; place the boxes to preserve the original appearance, or
move them into the surrounding white space now that the panel is sharper.

---

## 8. Validation against the original panel

The new PDF was rasterised to 1292 px (the original's width) and compared
pixel-wise with `_reference_current/CURRENT_fig8a_UCSC.png`.

**Matches:**

* Same assembly, same genomic window, same "500 bases" scale bar, and the same
  two ruler tick labels (12,766,000 / 12,766,500).
* Same 20 track rows in the same vertical order, same left labels
  (`HOOK2`, `ENCODE4 cCREs`, `RepeatMasker`, `CpG: 22`,
  `Adipocytes`/`Gastric-Ep`/`Head-Neck-Ep`/`Oligodend`, `OREG0025273`, the ten
  `… Merged` methylation rows, `1426`/`1`/`ReMap density`, `GH19J012765`).
* Identical HOOK2 exon/intron structure and strand arrows.
* Methylation-atlas rows are greyscale-shaded (dense mode) in both — confirmed by
  pixel sampling of the original, which contains no saturated colour in that
  block apart from the manually drawn red line.
* **Feature positions, expressed as a fraction of the data area, agree to
  ≤ 3 bp** (the label column is wider than the original's, so absolute pixel
  columns differ; the genomic mapping does not):

  | feature | original (fraction of data area) | new (fraction of data area) | max Δ |
  |---|---|---|---|
  | CpG island `CpG: 22` | 0.4207 … 0.5660 | 0.4219 … 0.5657 | 1.7 bp |
  | GeneHancer `GH19J012765` | 0.1222 … 0.7857 | 0.1237 … 0.7852 | 2.2 bp |
  | ORegAnno `OREG0025273` (clipped left) | 0.0000 … 0.6200 | 0.0019 … 0.6203 | 2.9 bp |
  | 2nd `Head-Neck-Ep` unmethylated block | 0.6191 … 0.6909 | 0.6194 … 0.6903 | 0.8 bp |

**Deliberate improvement over the original (documented, not a regression):**

* **Methylation-atlas row labels are no longer truncated.** The original panel
  clips them on the left — it reads "od Granulocytes Merged", "eolar Epithelium
  Merged", "reas Beta Cells Merged", "ver Hepatocytes Merged", "olon Epithelium
  Merged" and "Cardiomyocytes Merged" (verified by zooming into the reference
  PNG). `hgt.labelWidth=31` widens the label column so all ten labels render in
  full: *Neurons Merged, Heart Cardiomyocytes Merged, Endothelial Merged, Blood
  T Cells Merged, Blood B Cells Merged, Blood Granulocytes Merged, Lung Alveolar
  Epithelium Merged, Pancreas Beta Cells Merged, Liver Hepatocytes Merged, Colon
  Epithelium Merged.*
  **Cost:** the label gutter grows from 12.6 % to 19.3 % of the panel width, so
  the data area is compressed from 87.4 % to 80.7 % (horizontal scale 0.693 px/bp
  instead of 0.750 px/bp at 1292 px). The 500-base ruler ticks and the printed
  text size are unaffected. Raising `pix` instead of the label width would have
  preserved the proportions but shrunk the printed text below ~4.5 pt, which was
  judged worse for a print figure.

**Remaining small differences (all cosmetic):**

1. **ReMap density track is ~20 % shorter.** In the original it occupies ~149 px
   of 834; here ~122 px of 824. UCSC hard-caps this bigWig at `heightPer=128`
   (larger values are silently clamped — verified by request), so the original's
   taller box is not reproducible from the current trackDb. The density curve's
   shape, peak position and 1…1426 axis are identical. Overall aspect ratio is
   therefore 1.570 vs the original's 1.549 (1.3 % difference).
2. **Centre-label wording.** UCSC now writes *"Methylation Atlas: All unmethylated
   regions"*; the original screenshot (older UCSC release) read *"All unmethylated
   regions"*. A UCSC trackDb label change, not a data change.
3. At low zoom the ReMap density fill can show faint vertical seams — an
   anti-aliasing artefact of adjacent abutting vector rectangles (measured gap
   = −0.0004 pt, i.e. the rectangles actually overlap). At the delivered 600 dpi
   and at any realistic print size the fill is solid (worst-case seam luminance
   25/255 when downsampled to 300 dpi at 7 in).

---

## 9. Changelog

**rev. 2 (2026-07-25)**

* Added `hgt.labelWidth=31` — fixes truncation of the longest Methylation-Atlas
  row labels (a defect inherited from the original panel, see §8). Label column
  12.8 % → 19.3 % of width; ruler tick spacing and text size unchanged.
* Page box is now tight in **/MediaBox as well as /CropBox** (`[0 58.97 540 403]`).
  rev. 1 set only /CropBox, so MediaBox-based tools still rendered the full
  540 × 403 pt sheet with ~15 % empty space at the bottom.
* PDF page size 540 × 344.03 pt and PNG 4500 × 2867 px are unchanged from rev. 1
  (label width does not affect vertical layout).
* rs897804 red-line x-position updated: 2575 px (PNG) / 309.0 pt (PDF); it was
  2421 px / 290.5 pt in rev. 1.

**rev. 1 (2026-07-25)** — initial vector export.
