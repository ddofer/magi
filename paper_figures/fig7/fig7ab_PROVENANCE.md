# Figure 7a / 7b — print-quality regeneration

Prepared 2026-07-25. Replaces the low-resolution screenshots currently embedded in
`MAGI-NAR Fig1-8-23072026.pptx` (slide 8, `ppt/media/image21.png` = panel a,
`ppt/media/image23.png` = panel b; ~276 and ~199 effective DPI respectively).

Both panels are now **true vector** (SVG + PDF, live text, no rasterised elements —
verified with `pdfimages -list`: 0 embedded images in every PDF), with 600 dpi PNG
fallbacks.

---

## 0. Discrepancies requiring an author decision

Two factual mismatches between the published figure and the underlying Ensembl
annotation were found while regenerating panel a. **Neither has been silently
"corrected" in the artwork** — the regenerated panel reproduces exactly what was
published, so the authors must decide.

### 0.1 The panel shows 16 transcripts; the manuscript text says 17

> Supporting text: *"Fig. 7a analyzed all 17 transcripts reported for GLA gene
> (ENSG00000102393)."*

**Evidence.** The published panel image contains exactly **16** transcript rows
(counted directly off `CURRENT_fig7a_transcripts.png`; 16 row bands detected, 16 ID
labels). Those 16 IDs *with their exact version suffixes* match **Ensembl release
109** (feb2023 archive), which lists **16** GLA transcripts and no more:

| Ensembl release | GLA transcripts | Matches published panel? |
|---|---|---|
| 107 – 109 (to Feb 2023) | **16** | **yes — set and versions identical** |
| 110 – 115 | 17 | no (adds ENST00000710365.1; two versions bumped) |
| 116 (current) | 21 | no |

The 17th transcript, present from release 110 onwards, is
**ENST00000710365.1 (GLA-217, protein_coding, 8 exons, chrX:101,397,453–101,407,925)**.
It is **not** drawn in the published figure.

**Decision needed:** either change the text to "16 transcripts" (and cite Ensembl 109),
or add ENST00000710365.1 to the panel and re-state the release. Adding it would also
require re-checking the sentence *"in 4 of 5 transcripts intron retention (IR) between
Exon 5 to 6 were reported"* — release 109 annotates **5** retained-intron transcripts
(ENST00000466414.2, ENST00000468823.2, ENST00000675968.1, ENST00000479445.2,
ENST00000674142.1).

### 0.2 Row 12 is labelled `363 aa` but Ensembl gives it 362 aa

The right-margin PowerPoint label on row 12 (**ENST00000674634.2**, GLA-211) reads
`363 aa`. Ensembl 109 gives that transcript a **362 aa** translation; the 363 aa
protein-coding transcript is **ENST00000675592.1** (row 14), which is also labelled
`363 aa`. Panel b corroborates this independently: the corresponding TrEMBL entries
are **A0A6Q8PHD1 = 362 aa** and **A0A6Q8PGG0 = 363 aa** — two distinct isoforms, not
one length repeated.

```
ENST00000674634.2  GLA-211  protein_coding  Translation length = 362
ENST00000675592.1  GLA-212  protein_coding  Translation length = 363
```

**Decision needed:** change the row-12 label to `362 aa`, or document why 363 is
intended.

---

## 1. Deliverables

| File | Type | Nominal size | Notes |
|---|---|---|---|
| `fig7a_GLA_transcripts.svg` / `.pdf` | vector | 180.0 × 90.2 mm (510.24 × 255.75 pt) | **drop-in replacement** — reproduces the published panel's geometry |
| `fig7a_GLA_transcripts.png` | raster | 4251 × 2131 px @ 600 dpi | ≥ 710 dpi at the current 5.97-in PowerPoint placement |
| `fig7a_GLA_transcripts_trueScale.svg` / `.pdf` / `.png` | vector / raster | 180.0 × 83.9 mm; PNG 4251 × 1982 | strictly linear genomic x-axis, wider label margin (overlays must be repositioned) |
| `fig7b_GLA_MSA.svg` / `.pdf` | vector | 180.0 × 306.5 mm (510.24 × 868.75 pt) | **complete** alignment, all 470 columns, 9 blocks |
| `fig7b_GLA_MSA.png` | raster | 4251 × 7239 px @ 600 dpi | |
| `fig7b_GLA_MSA_asPublished.svg` / `.pdf` / `.png` | vector / raster | 180.0 × 238.1 mm (510.24 × 674.99 pt); PNG 4251 × 5624 | truncated at alignment column 385 (7 blocks), i.e. the same extent as the published screenshot |
| `make_fig7a.py`, `make_fig7b.py` | scripts | — | regeneration recipes (see §6) |
| `_data_fig7a_ensembl109.json` | data | — | cached Ensembl 109 REST payload |
| `_data_fig7b_input.fasta`, `_data_fig7b_clustalo.aln` | data | — | the exact 10 input sequences and the Clustal Omega alignment used |

Vector output is resolution-independent: the "nominal size" above is only the
`MediaBox`/`viewBox`; scaling the placed picture in PowerPoint or InDesign does not
degrade it.

**Font handling.** The **PDFs are self-contained** — Arial is subset-embedded as
TrueType (`ELUGPM+ArialMT`, `ELUGPM+Arial-BoldMT`), so they are the safest artwork to
send to the journal. The **SVGs keep text as live `<text>` elements referencing
`Arial, Liberation Sans, DejaVu Sans` by name** (editable, but the renderer must have
one of those fonts). Prefer the PDF for submission; use the SVG if the production
team needs to edit labels.

---

## 2. Panel a — data sources

**Gene** GLA (alpha-galactosidase A), `ENSG00000102393`, GRCh38,
**chrX:101,393,273–101,408,012, reverse strand (−1)**.
Because the gene is on the minus strand and the panel is drawn in ascending genomic
coordinate, exon numbering runs **7 → 1 left to right**.

**Source** Ensembl REST API, **release 109** (Feb 2023 archive), GRCh38.p13:

```
https://feb2023.rest.ensembl.org/lookup/id/ENSG00000102393?expand=1;content-type=application/json
```

Release 109 was chosen because it is the release whose GLA transcript set and
**version suffixes match the published panel exactly** (16 transcripts, incl.
`ENST00000674127.1` and `ENST00000486121.6`). From release 110 onwards Ensembl adds a
17th transcript and bumps two versions; release 116 (current) lists 21 transcripts.
See §5 for the implications.

**Transcript rows (top → bottom, exactly as published), Ensembl 109:**

| # | Transcript | Name | Biotype | Protein | PPT label |
|---|---|---|---|---|---|
| 1 | ENST00000218516.4 | GLA-201 | protein_coding | 429 aa | `429 aa` |
| 2 | ENST00000466414.2 | GLA-202 | retained_intron | — | `RI` |
| 3 | ENST00000480513.6 | GLA-205 | nonsense_mediated_decay | 199 aa | `NMD` |
| 4 | ENST00000468823.2 | GLA-203 | retained_intron | — | `RI` |
| 5 | ENST00000675968.1 | GLA-214 | retained_intron | — | `RI` |
| 6 | ENST00000479445.2 | GLA-204 | retained_intron | — | `RI` |
| 7 | ENST00000493905.6 | GLA-207 | nonsense_mediated_decay | 222 aa | `NMD` |
| 8 | ENST00000675799.1 | GLA-213 | nonsense_mediated_decay | 199 aa | `NMD` |
| 9 | ENST00000676156.1 | GLA-215 | protein_coding | 417 aa | `417 aa` |
| 10 | ENST00000674127.1 | GLA-209 | nonsense_mediated_decay | 198 aa | `NMD` |
| 11 | ENST00000674142.1 | GLA-210 | retained_intron | — | `RI` |
| 12 | ENST00000674634.2 | GLA-211 | protein_coding | **362 aa** | `363 aa` ⚠ |
| 13 | ENST00000676372.1 | GLA-216 | nonsense_mediated_decay | 270 aa | `NMD` |
| 14 | ENST00000675592.1 | GLA-212 | protein_coding | 363 aa | `363 aa` |
| 15 | ENST00000649178.1 | GLA-208 | protein_coding | 470 aa | `470 aa` |
| 16 | ENST00000486121.6 | GLA-206 | nonsense_mediated_decay | 199 aa | `NMD` |

Counts: 5 × retained_intron, 6 × NMD, 5 × protein_coding — matching the 5 `RI`,
6 `NMD` and 5 aa-labels present as PowerPoint text on slide 8.

**Exon coordinates** (GRCh38, ascending; canonical ENST00000218516.4):
101,397,803–101,398,099 (exon 7) · 101,398,370–101,398,567 (6) ·
101,398,785–101,398,946 (5) · 101,400,666–101,400,757 (4) ·
101,401,632–101,401,809 (3) · 101,403,811–101,403,985 (2) ·
101,407,710–101,407,925 (1). All other transcripts' exons are in the cached JSON.

**Export method** Not an Ensembl image export. The panel is re-drawn from the REST
exon coordinates with matplotlib (v3.10.9) → SVG/PDF/PNG. Style matched to the
original: exon boxes `#EDEDED` fill with `#565F64` 2-px stroke, introns as
`#565F64` 1.8-px lines, transcript IDs in Arial in a left-hand column.

**Coordinate mapping (default variant).** The published screenshot is *not* on a
linear genomic axis: it is linear at 0.147853 px/bp up to ~chrX:101,402,800 and then
compressed three times inside the large 3′ introns. To keep the figure a drop-in
replacement, this mapping was recovered by least-squares fitting all 210 exon edges
measured off the original PNG:

```
x(g) = 0.14785300894216793 * (g - 101397803) + 216.43247735032767
       - 148.38998773698228   if g > 101,402,800
       - 18.864733392321828   if g > 101,404,500
       - 218.41863848567635   if g > 101,406,500
```
(units = pixels of the original 1648 × 688 PNG; fit rms 0.33 px, max 2.58 px.)

Verified: the regenerated panel has the same number of exon boxes on every row and
every box edge lands within **4 px of 1358** (< 0.3 % of panel width) of the original.

The `_trueScale` variant discards this mapping and uses a single linear genomic axis
over chrX:101,393,273–101,407,990.

---

## 3. Panel b — data sources

**Sequences** UniProtKB, retrieved 2026-07-25 from
`https://rest.uniprot.org/uniprotkb/{accession}.fasta`
(cached verbatim in `_data_fig7b_input.fasta`). Row order as published:

| # | Accession | Status | SV | Length | Description |
|---|---|---|---|---|---|
| 1 | A0A669KB83 | TrEMBL | 2 | 222 | Alpha-galactosidase |
| 2 | A0A3B3IRU3 | TrEMBL | 1 | 199 | Alpha-galactosidase |
| 3 | V9GYN5 | TrEMBL | 1 | 222 | Alpha-galactosidase |
| 4 | A0A6Q8PGG0 | TrEMBL | 1 | 363 | Alpha-galactosidase |
| 5 | A0A6Q8PFA9 | TrEMBL | 1 | 417 | Alpha-galactosidase |
| 6 | **P06280** | **SwissProt** | 1 | **429** | Alpha-galactosidase A (AGAL_HUMAN) — canonical |
| 7 | A0AA34QW02 | TrEMBL | 1 | 454 | Alpha-galactosidase |
| 8 | A0A3B3IUC4 | TrEMBL | 1 | 470 | Alpha-galactosidase |
| 9 | A0A6Q8PHD1 | TrEMBL | 2 | 362 | Alpha-galactosidase |
| 10 | A0A6Q8PHM8 | TrEMBL | 1 | 270 | Alpha-galactosidase |

All ten are *Homo sapiens* (OX=9606), GN=GLA. The lengths in this column are exactly
the ten numbers currently overlaid as PowerPoint text in the right margin
(222, 199, 222, 363, 417, 429, 454, 470, 362, 270) — confirming the accession set is
the published one.

**Alignment** Clustal Omega **1.2.4**, run through the EBI Job Dispatcher REST API
(the same engine behind the UniProt "Align" tool that produced the original figure):

```
POST https://www.ebi.ac.uk/Tools/services/rest/clustalo/run
     email=<user>  sequence=@_data_fig7b_input.fasta  outfmt=clustal_num  order=input
job id: clustalo-R20260725-165508-0610-24646492-p1m
GET  https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/<jobid>/aln-clustal_num
```
All other parameters at their service defaults (guide tree from mBed, 0 iterations,
BLOSUM/Gonnet default matrix). Result cached in `_data_fig7b_clustalo.aln`.

* alignment length **470 columns**; **170 columns fully identical** across all 10 rows.
* This alignment is **column-for-column identical to the published one** — every
  residue counter visible in the original screenshot is reproduced
  (55 / 69·69·69·69·69·69·94·110·69·69 / 124…149·165 / 179…170…204·220 /
  222·199·222·234·234·234·259·275·234·234 / …), including the two TrEMBL-specific
  insertions in A0AA34QW02 (25 aa) and A0A3B3IUC4 (41 aa) in block 2.

**Row labels** are the **bare UniProt accession** (`A0A669KB83`, `P06280`, …), not the
original viewer's `db|ACC|ENTRY_HUMAN` triplet. The prefix and the `_HUMAN` entry name
are fully recoverable from the table above and from `_data_fig7b_input.fasta`; the
full-form strings as they appear in the published screenshot are:

```
tr|A0A669KB83|A0A669KB83_HUMAN    tr|A0AA34QW02|A0AA34QW02_HUMAN
tr|A0A3B3IRU3|A0A3B3IRU3_HUMAN    tr|A0A3B3IUC4|A0A3B3IUC4_HUMAN
tr|V9GYN5|V9GYN5_HUMAN            tr|A0A6Q8PHD1|A0A6Q8PHD1_HUMAN
tr|A0A6Q8PGG0|A0A6Q8PGG0_HUMAN    tr|A0A6Q8PHM8|A0A6Q8PHM8_HUMAN
tr|A0A6Q8PFA9|A0A6Q8PFA9_HUMAN    sp|P06280|AGAL_HUMAN
```

Dropping the redundant prefix/suffix shrinks the label gutter from 149 to 73.4
reference px and hands that width to the alignment: cell pitch rises from 11.7636 to
**13.230** reference px (**+12.5 %**) and the residue type with it. The SwissProt
canonical entry `P06280` is set **bold**, and its residue counter is set in **red**
(`#BB2739`), preserving the emphasis the authors added by hand.

**Rendering** matplotlib 3.10.9, Arial, 55 alignment columns per block, cell pitch
13.230 reference px, row pitch 13.3 reference px. Column shading (identity fraction
computed over **all 10 rows**, gaps counting as mismatches) uses the exact colours
sampled from the published figure:

| Identity | Colour |
|---|---|
| 100 % | `#9293F9` |
| ≥ 80 % | `#A6A7FA` |
| ≥ 60 % | `#D5D8FC` |
| ≥ 40 % | `#EEF0F1` |
| < 40 % | none (`#FCFEFF`) |

Verified cell-by-cell against the original: all five tiers reproduce the reference
RGB values exactly.

**Signal peptide** P06280 SIGNAL 1–31 (UniProt). Because the alignment is gap-free in
that region, this is alignment columns **1–31**, drawn as a `#BB2739` box around all
ten rows of block 1 plus a `#C23D4D` bar on the `A0A669KB83:Signal` feature track —
as in the published panel.

The feature track is emitted **only for blocks whose column range intersects the
annotated feature**, i.e. block 1 only; blocks 2–9 carry no feature row and reserve
no vertical space for one. (The UniProt viewer repeats the track label under every
block even where the feature is absent; that produces orphan labels with nothing
beside them, so it is not reproduced.) Because the label no longer has a 149 px
gutter to sit in, it is placed immediately to the **right of the red bar** rather than
in the left margin — at the same type size as the row labels, instead of the ~4.6 pt
it would need to fit the narrowed gutter.

**Met267Ile position** — see §4.

---

## 4. Manual PowerPoint overlays the authors must re-apply

Everything below is native PowerPoint text/shape on slide 8, **not** part of the
image, and is therefore *not* reproduced in these files.

### Panel a
* exon-number boxes `7 6 5 4 3 2 1`
* the **red vertical line** for c.801G>A
* right-margin labels: `429 aa`, `RI` ×5, `NMD` ×6, `417 aa`, `363 aa` ×2, `470 aa`
* panel letter `a`

**Where the red line goes.** NM_000169.3:c.801 is the last nucleotide of exon 5,
**chrX:101,398,785 (GRCh38)** — the low-coordinate (3′-in-transcript) edge of the
exon-5 box. In the regenerated panel that is:

* x = **361.6** in original-screenshot pixels
* = **26.63 %** of the panel width, measured from the left edge of the image
* i.e. flush with the **left edge of the third exon box from the left**

Because the geometry was matched to < 0.3 % of panel width, the existing PowerPoint
line, exon numbers and right-margin labels should still line up after replacing the
picture — **but the picture's crop must be cleared** (see §6).

### Panel b
* right-margin totals `222 199 222 363 417 429 454 470 362 270` + `aa`
  → **now redundant, delete them**: the final block of `fig7b_GLA_MSA.*` carries
  exactly these numbers as real, in-figure text, with `429` already set in red on the
  bold `P06280` row. (In `fig7b_GLA_MSA_asPublished.*` the alignment stops at column
  385, so the last visible counters are the running totals 222/199/222/278/332/**344**
  /369/385/344/270, not the full lengths — keep the overlay if you use that variant.)
* the `Met` / `Ile` text labels
* the **blue up-arrow** marking Met267Ile
* panel letter `b`

**Where the blue arrow goes.** P06280 Met267 is **alignment column 308**, i.e.
**block 6 (of 9), position 33 within the block**, row 6 (`sp|P06280|AGAL_HUMAN`).
The residue sits at the end of the conserved motif `…AGPGGWNDP**M**` — the last
residue encoded by in-frame exon 5 (A0A6Q8PGG0 terminates immediately after it and
A0A6Q8PHM8 diverges immediately after it, which is what makes the exon boundary
visible in the alignment).

Coordinates as fractions of the placed picture (0 = left/top edge):

| Variant | x (column centre) | y (top / bottom of the P06280 row) |
|---|---|---|
| `fig7b_GLA_MSA` (9 blocks) | 0.5953 | 0.6137 / 0.6227 |
| `fig7b_GLA_MSA_asPublished` (7 blocks) | 0.5953 | 0.7899 / 0.8015 |

Place the arrow tip at (x, y_bottom) pointing up.
(The P06280 row is the 6th of 10 in every block; its counter is the red one.)

---

## 5. Differences from the published screenshots — read this

**Panel a**

1. **16 transcripts, not 17.** The published panel shows 16; the manuscript text says
   "all 17 transcripts reported for GLA". Ensembl 109 (which matches the panel
   exactly) lists 16. From release 110 a 17th appears —
   **ENST00000710365.1 (GLA-217, protein_coding, 8 exons, chrX:101,397,453–101,407,925)**.
   Either the text should say 16, or that transcript should be added (it is *not*
   in the published image). Current Ensembl (116) lists 21 transcripts.
   *Nothing was silently added or removed: the regenerated panel contains the same
   16 transcripts, in the same order, as the published one.*
2. **Version drift.** If the authors prefer current Ensembl, note
   `ENST00000674127.1 → .2`, `ENST00000486121.6 → .7`, `ENST00000649178.1 → .2`.
   The regenerated panel keeps the published `.1`/`.6`/`.1` labels (release 109).
3. **`363 aa` label on row 12.** Ensembl gives ENST00000674634.2 a **362 aa**
   translation (ENST00000675592.1 is the 363 aa one). The corresponding TrEMBL
   entries in panel b are likewise 362 (A0A6Q8PHD1) and 363 (A0A6Q8PGG0) aa.
   Consider changing the row-12 label to `362 aa`.
4. **ENST00000674142.1** has a distal exon at chrX:101,393,273–101,393,829 that falls
   outside the displayed window. As in the original, it is represented only by the
   intron line running off the left edge of the plot area. It is drawn in full in the
   `_trueScale` variant.
5. **Non-linear axis.** The default variant deliberately preserves the original's
   piecewise x-axis (§2) so the overlays still fit. If a strictly proportional axis is
   preferred, use `fig7a_GLA_transcripts_trueScale.*` and reposition the overlays.
6. The duplicate right-hand column of transcript IDs present in the original PNG is
   omitted — the PowerPoint frame already cropped it away (`srcRect r="17602"`).

**Panel b**

1. **The published screenshot is truncated at alignment column 385** (7 of 9 blocks);
   the C-terminal ~85 columns were cut off. `fig7b_GLA_MSA.*` shows the **complete**
   470-column alignment (9 blocks), which is why it is taller than the original.
   `fig7b_GLA_MSA_asPublished.*` is provided if the original extent must be kept.
2. The UniProt-viewer chrome (row check-boxes, entry icons, row rules, the grey
   "selected row" highlight on A0A669KB83) is not reproduced — it is user-interface,
   not data. The check-box/icon column was already cropped off in PowerPoint
   (`srcRect l="4706"`).
3. **Row labels are bare accessions** (`A0A669KB83`) rather than the viewer's
   `tr|A0A669KB83|A0A669KB83_HUMAN`. The prefix and entry name carry no information
   the reader needs and cost 12.5 % of the alignment width; the full forms are
   tabulated in §3 and preserved in `_data_fig7b_input.fasta`. `P06280` is **bold**
   with a **red** residue counter, marking the SwissProt canonical entry (this
   replaces the viewer's gold icon / ticked check-box, and makes the authors'
   "429 in red" PowerPoint overlay redundant).
4. **The `A0A669KB83:Signal` feature row appears once**, under block 1 — the only
   block whose columns intersect the signal peptide. The original viewer repeated the
   label under every block even where the red bar was absent. The label sits to the
   right of the bar rather than in the left gutter (see §3).
5. **Type size.** The layout is still type-hungry (55 columns per row). Sizes in the
   delivered vector files, expressed at the 180 mm nominal width:

   | Element | reference px | pt @ 180 mm | pt @ 111 mm (current placement) |
   |---|---|---|---|
   | residues | 11.81 | 6.9 | 4.3 |
   | row labels / feature label | 11.0 | 6.5 | 4.0 |
   | residue counters | 11.0 | 6.5 | 4.0 |

   The bare-accession change already bought +12.5 % over the published layout, but the
   *current* PowerPoint placement (4.36 in ≈ 111 mm) still puts the residues at
   ≈ 4.3 pt, below most journal minimums. **Recommendation: place panel b at
   ≥ 160 mm wide** (≈ 6.1 pt residues) or full page width 180 mm (≈ 6.9 pt). The
   vector files scale losslessly, so this is purely a layout decision. Further gains
   would require fewer columns per block (a taller figure) or dropping the feature
   track.

**Nothing in the figure legend needs to change.** Every statement in the published
legend remains accurate for these files.

---

## 6. How to regenerate

```bash
cd NAR_figures_highres/fig7

# Panel a  (uses cached _data_fig7a_ensembl109.json; delete it to re-fetch)
python3 make_fig7a.py                 # drop-in replacement
python3 make_fig7a.py --truescale     # linear genomic axis variant

# Panel b  (uses cached _data_fig7b_input.fasta + _data_fig7b_clustalo.aln)
python3 make_fig7b.py                 # complete alignment
python3 make_fig7b.py --as-published  # truncated at column 385
```

Requirements: Python 3.10, matplotlib ≥ 3.8, network access only if the cached data
files are deleted. Arial is picked up from `/mnt/c/Windows/Fonts/arial.ttf` under WSL;
otherwise matplotlib falls back to Liberation Sans / DejaVu Sans.

To re-run the alignment from scratch:

```bash
for a in A0A669KB83 A0A3B3IRU3 V9GYN5 A0A6Q8PGG0 A0A6Q8PFA9 \
         P06280 A0AA34QW02 A0A3B3IUC4 A0A6Q8PHD1 A0A6Q8PHM8; do
  curl -s "https://rest.uniprot.org/uniprotkb/$a.fasta"
done > _data_fig7b_input.fasta

JOB=$(curl -s -X POST \
  --form email=YOU@EXAMPLE.ORG \
  --form sequence="$(cat _data_fig7b_input.fasta)" \
  --form outfmt=clustal_num --form order=input \
  https://www.ebi.ac.uk/Tools/services/rest/clustalo/run)
curl -s "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/$JOB/aln-clustal_num" \
  -o _data_fig7b_clustalo.aln
```

### Swapping the images into the PowerPoint deck

Slide 8 of `MAGI-NAR Fig1-8-23072026.pptx`:

* panel a = `rId3` → `ppt/media/image21.png`, currently cropped
  `srcRect t="1070" r="17602"`, placed 5.97 × 2.99 in at (0.97, 0.75) in.
  **Insert the new file and reset the crop to zero** — the new image already *is*
  the cropped region, so the overlays keep their positions.
* panel b = `rId5` → `ppt/media/image23.png`, currently cropped
  `srcRect l="4706" t="723"`, placed 4.36 × 6.11 in at (8.41, 0.90) in.
  Reset the crop to zero, re-proportion the frame and delete the right-margin count
  overlay. Aspect ratios (width : height): published crop 1 : 1.403 →
  `fig7b_GLA_MSA_asPublished` **1 : 1.323**, `fig7b_GLA_MSA` (full alignment)
  **1 : 1.703**.
* PowerPoint cannot place SVG/PDF reliably for print; use the 600-dpi PNGs in the
  deck and supply the **PDF or SVG** to the journal as the production artwork.
