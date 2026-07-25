# Figure 5, panels b and c — provenance and regeneration recipe

Print-quality replacements for the low-resolution screenshots currently used in
Figure 5b and 5c of the MAGI manuscript (NAR submission).

Generated: 2026-07-25. All source data from **Ensembl release 116 (June 2026,
GENCODE 50)**, assembly **GRCh38.p14**, via `https://www.ensembl.org` and
`https://rest.ensembl.org`.

---

## 1. Identifiers, validated against the Ensembl REST API

| Item | Value | Check |
|---|---|---|
| Gene | **MAT1A**, `ENSG00000151224` | `rest.ensembl.org/lookup/symbol/homo_sapiens/MAT1A` |
| Locus (release 116) | chr10:**80,269,780-80,290,150**, **reverse (−) strand** | same |
| Locus (release 115, used for the published panel) | chr10:80,269,780-**80,290,126** | see §5 |
| Variant | **rs118204006**, chr10:**80,274,599 C>T** (GRCh38) | `rest.ensembl.org/variation/human/rs118204006` |
| Variant synonyms | VCV000001211, RCV000001270, VAR_031244, OMIM 610550.0010 | same |
| Clinical significance | pathogenic / likely pathogenic / uncertain significance | same |
| HGVS (MANE Select) | `NM_000429.3:c.1006G>A`, `p.Gly336Arg`; codon `Gga>Aga` | VEP REST |

### Transcripts shown

| Row | Name | Transcript ID | Protein | Length | Exons | Consequence of rs118204006 |
|---|---|---|---|---|---|---|
| 1 | MAT1A-205 | `ENST00000871619` | ENSP00000541678 | **370 aa** | 10 | missense, p.Gly311Arg (CDS 931) |
| 2 | MAT1A-208 | `ENST00000871622` | ENSP00000541681 | **329 aa** | 8 | missense, p.Gly270Arg (CDS 808) |
| 3 | MAT1A-207 | `ENST00000871621` | ENSP00000541680 | **373 aa** | 9 | **splice_polypyrimidine_tract_variant + intron_variant** |
| 4 | MAT1A-201 (MANE Select, gold) | `ENST00000372213` | ENSP00000361287 | **395 aa** | 9 | missense, **p.Gly336Arg** (CDS 1006) |
| — | MAT1A-206 (panel c only) | `ENST00000871620` | ENSP00000541679 | 393 aa | 9 | missense, p.Gly334Arg (CDS 1000) |

> **ERROR IN THE CURRENT FIGURE — please correct.** Panel b currently labels
> MAT1A-201 as "**385 aa**". Ensembl 116 (`ENSP00000361287`), UniProt **Q00266**
> and RefSeq `NP_000420.1` all give **395 aa**. The other three labels
> (370 / 329 / 373 aa) are correct.

MAT1A-207 is intronic at this position because its exon 8 is truncated
(chr10:80,274,520-80,274,**587**) relative to the other transcripts
(chr10:80,274,520-80,274,**653**); the variant at 80,274,599 therefore falls in
intron 7/8 of MAT1A-207. This is the observation the legend refers to
("MAT1A-207 may influence splicing regulatory elements with splicing donors used
differently than the canonical one"), and it is confirmed by VEP on release 116.

Exon coordinates used for the zoom (MAT1A-201 numbering; = exons 7-9 of the
10-exon MAT1A-205, which is the numbering used in the published legend):

| MAT1A-205 exon | MAT1A-201 exon | GRCh38 coordinates |
|---|---|---|
| 9 (contains variant) | 8 | 80,274,520-80,274,653 (MAT1A-207: 80,274,520-80,274,587) |
| 8 | 7 | 80,275,017-80,275,199 |
| 7 | 6 | 80,276,376-80,276,594 (MAT1A-205/-208: 80,276,376-80,276,519) |

---

## 2. Panel b — transcript models

**Source page:** Ensembl *Region in detail* (`Location/View`), bottom panel
(`ViewBottom`), track **"Genes (Comprehensive set from GENCODE 50)"**
(config key `transcript_core_ensembl`) with renderer
**"Expanded with labels"** (`transcript_label`).

URLs:

* full view — `https://www.ensembl.org/Homo_sapiens/Location/View?r=10:80269780-80290126;g=ENSG00000151224`
* zoom (exons 7-9) — `https://www.ensembl.org/Homo_sapiens/Location/View?r=10:80273870-80276800;g=ENSG00000151224`

The full view deliberately uses the **release-115 gene span**
(`80,269,780-80,290,126`) so that the horizontal scale is identical to the
published panel.

**Export:** Ensembl's own *Export image* dialog, "Custom image" → **SVG**
(true vector). Programmatically this is a `POST` to
`https://www.ensembl.org/Homo_sapiens/ImageExport/ImageOutput` with

```
filename=MAT1A_full.svg  format=custom  image_format=svg  resize=  scale=
r=10:80269780-80290126   data_type=Location  component=ViewBottom
data_action=View         db=core  strain=0  decodeURL=1
extra={"highlightedTracks":[]}   submit=Download
```

(The same dialog also offers PDF and PNG; `image_format=png` plus `scale=5`
gives a 5× raster.) Native SVG canvas: **1600 px wide**.

**Post-processing** (`extract_rows.py`, kept with the scratch files, reproduced
in §6): Ensembl release 116 lists **27** MAT1A transcripts, so the raw export is
1600 × 1914 px with 27 transcript rows. The script

1. finds the `<text>` label of each requested transcript and derives its row band
   (rows are on a strict 44 px pitch, glyph 9 px high, label baseline +14 px);
2. copies every `<rect>` / `<path>` / `<text>` whose vertical centre falls in that
   band into a `<g transform="translate(0,dy)">` (path data is **never** rewritten —
   Ensembl mixes absolute `M` with relative `l` commands);
3. re-draws the full-height vertical grid lines at the new height;
4. drops the pale gene-highlight band (added because `g=` is in the URL);
5. splits Ensembl's two-line `<text>` nodes (they contain a literal newline that
   SVG would otherwise collapse onto one line) into two `<text>` elements 13 px apart;
6. crops the track-name gutter with `viewBox="148 0 1452 194"`.

Result: exactly the four transcripts of the published panel, in the published
order (MAT1A-205, -208, -207, -201), rendered as authentic Ensembl vector art.

### Files

| File | Type | Size |
|---|---|---|
| `fig5b_MAT1A_transcripts.svg` | **vector** | 1452 × 194 pt |
| `fig5b_MAT1A_transcripts.pdf` | **vector**, embedded fonts | 1452 × 194 pt |
| `fig5b_MAT1A_transcripts.png` | raster | 6000 × 802 px |
| `fig5b_MAT1A_transcripts_zoom_ex7-9.svg` | **vector** | 1452 × 194 pt |
| `fig5b_MAT1A_transcripts_zoom_ex7-9.pdf` | **vector**, embedded fonts | 1452 × 194 pt |
| `fig5b_MAT1A_transcripts_zoom_ex7-9.png` | raster | 6000 × 802 px |

At 180 mm print width the vector files are resolution-independent; the PNGs are
≈ 850 dpi at that width.

---

## 3. Panel c — transcript sequence alignment

**Source page:** Ensembl *Gene → Transcript comparison*.

```
https://www.ensembl.org/Homo_sapiens/Gene/TranscriptComparison?db=core;g=ENSG00000151224
  ;r=10:80269780-80290150
  ;t1=ENST00000372213;t2=ENST00000871619;t3=ENST00000871620
  ;t4=ENST00000871621;t5=ENST00000871622
```

The transcript set is chosen with the **"Select transcripts"** button in the
left-hand menu (`MultiSelector/Gene/TranscriptComparisonSelector`); the selection
is then encoded as the ordered `t1…t5` URL parameters above, so the URL alone
reproduces the view. (Repeating a plain `t=` parameter does **not** work.)

Display configuration — Ensembl defaults, unchanged:
`display_width=120` bp/row, `snp_display=on` (show variants),
`hide_long_snps=on` (hide variants > 10 bp), `hide_rare_snps=off`,
`consequence_filter=off`, `line_numbering=sequence`, `exons_only=off`,
`flank5_display=0`, `flank3_display=0`.

**The row shown is gene position 15,481-15,600.** With `line_numbering=sequence`
the coordinate is the offset along the gene's genomic sequence read 5'→3' on the
transcript (−) strand, i.e. position 1 = the highest genomic coordinate of the
gene. In release 116 (gene 5' end = chr10:80,290,150):

* row 15,481-15,600 = chr10:**80,274,670** down to chr10:**80,274,551**
* the variant, chr10:80,274,599, is gene position **15,552** = **column 72** of the row
* row sequence:
  `CTGCACCTTGCTGGAAGGTTTCCTATGCCATTGGTGTGGCCGAGCCGCTGTCCATTTCCATCTTCACCTACGGAACCTCTCAGAAGACAGAGCGAGAGCTGCTGGATGTGGTGCATAAGA`

The published panel shows the same row number, 15481, but its sequence starts
24 nt further 3' (see §5).

The colour key reproduced above the alignment is Ensembl's own
variant-consequence key (the "Variants" block of `_adornment_key`), colours
unchanged:

| Consequence | Colour | | Consequence | Colour |
|---|---|---|---|---|
| 3 prime UTR | `#7ac5cd` | | Splice acceptor | `#FF581A` |
| 5 prime UTR | `#7ac5cd` | | Splice donor | `#FF581A` |
| Coding sequence | `#458b00` | | Splice region | `#ff7f50` |
| Frameshift | `#9400D3` | | Start lost | `#ffd700` |
| Inframe deletion | `#ff69b4` | | Stop gained | `#ff0000` |
| Intronic | `#02599c` | | Stop lost | `#ff0000` |
| Missense | `#ffd700` | | Synonymous | `#76ee00` |

**Export:** the Transcript comparison page is HTML text, not an Ensembl image,
so it has no "Export image" button. The page was rendered in headless Chromium
(Playwright), a DOM `Range` was placed over the five transcript lines of the
15,481 block and cloned into a standalone `<pre class="text_sequence">`, the
"Variants" section of the colour key was cloned above it, everything else was
hidden, and the result was exported with Chromium **print-to-PDF** (true vector,
live text) plus a `device_scale_factor=6` element screenshot.

### IMPORTANT — the live page can no longer produce this panel as published

On **Ensembl release 116 the Transcript-comparison page renders per-base variant
markup for MAT1A-201 only**. The other four rows come back as plain
"translated sequence" (blue) / "intron" (grey) text with no consequence
backgrounds. This was verified exhaustively: the page was left open for 10 min,
the target block was scrolled into view, and the lazily-applied "adornment"
markup was polled every 15 s — the count of marked bases stayed at
`MAT1A-201: 61, MAT1A-205/206/207/208: 0` and the global styled-element count
froze at 6627, i.e. adornment had *completed*, not stalled. The release-115
figure did have all five rows marked, and Ensembl's release-115 archive
(`sep2025.archive.ensembl.org`) returns HTTP 403 to automated clients, so it
could not be re-rendered there.

Two versions are therefore supplied:

**(a) Primary — recomputed from Ensembl data (`build_fig5c.py`).**
Same window, same five transcripts, same Ensembl colour key, but the per-base
consequence for *every* transcript is taken from the Ensembl REST API instead of
from the web page:

* variants in the window: `rest.ensembl.org/overlap/region/human/10:80274551-80274670?feature=variation`
  → **65** features;
* per-transcript consequences: `POST rest.ensembl.org/vep/human/id` for all 65 IDs
  → 59 return `transcript_consequences`; the 6 HGMD features
  (CD962068, CI158698, CM001220, CM001221, CM106330, CM950788) return none and
  fall back to their region-level `consequence_type`, resolved per transcript by
  exon/intron state (which is what Ensembl's own page does with them);
* per base, the **most severe** consequence wins, using Ensembl's standard
  severity order; the colour and text colour are Ensembl's own;
* unmarked bases are coloured by Ensembl's exon/intron key —
  translated `#1044ee`, UTR `#cd6839`, non-coding exon `#333333`, intron `#aaaaaa`;
* variant bases are underlined, exactly as Ensembl does.

The window sequence was checked byte-for-byte against the live Ensembl page row
(identical), and the reading frame of MAT1A-201/-205/-206/-208 is the same at
this locus (CDS 1006 / 931 / 1000 / 808, all ≡ 1 mod 3), which is why those four
rows are identically coloured — exactly as in the published panel.

**(b) Secondary — the authentic release-116 live page** capture, supplied
unedited so the limitation above is documented and auditable.

### Files

| File | Type | Size |
|---|---|---|
| `fig5c_MAT1A_transcript_alignment.svg` | **vector** (a) | 1237 × 193 pt |
| `fig5c_MAT1A_transcript_alignment.pdf` | **vector**, embedded fonts, no rasters (a) | 1234 × 197 pt |
| `fig5c_MAT1A_transcript_alignment.png` | raster, 6× (a) | 7404 × 1188 px |
| `fig5c_MAT1A_transcript_alignment_fromSVG.png` | raster rendered from the SVG (a) | 7400 × 1155 px |
| `fig5c_MAT1A_transcript_alignment.html` | the HTML source of (a), for restyling | — |
| `fig5c_MAT1A_transcript_alignment_ensembl116_livepage.pdf` | **vector** (b), Ensembl live page | 1157 × 145 pt |
| `fig5c_MAT1A_transcript_alignment_ensembl116_livepage.png` | raster, 6× (b) | 6948 × 870 px |

---

## 4. Manual PowerPoint overlays the authors must RE-APPLY

None of the following are produced by Ensembl; they are annotations on the
original slide and must be redrawn on top of the new assets:

**Panel b (full view)**
* the four amino-acid labels at the left — `370 aa`, `329 aa`, `373 aa`, and
  **`395 aa`** (currently, incorrectly, `385 aa`);
* the **red dashed vertical line** at the variant. Its position is
  chr10:80,274,599, i.e. at fraction
  `(80274599 − 80269780) / (80290126 − 80269780) = 0.2369` of the *plotted*
  width (the plot area is the whole width of the supplied SVG/PNG);
* the **blue zoom box** and the two blue connector lines down to the inset.

**Panel b (zoom inset)**
* the outer blue frame;
* the red dotted vertical line at the variant — fraction
  `(80274599 − 80273870) / (80276800 − 80273870) = 0.2488` of the plotted width;
* the `SNP: rs118204006` callout box.

**Panel c**
* the **red rectangle** around the variant base — it is **column 72** of the
  120-column row (the base `G`, in `...ACCTAC[G]GAACC...`);
* the **red arrow** above it and the `SNP: rs118204006` label;
* the horizontal **"Transcript strand (−)"** arrow beneath the alignment.
  *Note: the current figure reads "**Trsnscript strsnd (-)**" — a typo that
  should be fixed to "Transcript strand (−)".*

---

## 5. Differences from the currently-published panels — read this

1. **Ensembl release.** The published panels were made on **Ensembl 115**
   (Sep 2025). This is provable from the panel-c coordinate: the row labelled
   15481 in the published figure begins at chr10:80,274,646, which implies a gene
   5' end of chr10:**80,290,126** — the release-115 span. In release 116 the gene
   was extended 24 bp by a new NMD transcript (MAT1A-226 / `ENST00001130619`,
   5' end 80,290,150), so the same row number now starts 24 nt earlier and the
   variant sits in column 72 instead of column 48. Ensembl **archive sites
   (`e115.ensembl.org` / `sep2025.archive.ensembl.org`) return HTTP 403 to
   automated clients**, so the panels could not be regenerated on release 115;
   everything here is release 116.
2. **Transcript count.** Release 115 annotated 8 MAT1A transcripts; release 116
   annotates **27**. The panel-b export therefore contains 27 rows and the four
   published ones are extracted programmatically (§2). In release 116 a new
   transcript (MAT1A-227 / `ENST00001140301`) sits between MAT1A-207 and
   MAT1A-201 in Ensembl's row order; the extraction removes it so the four rows
   remain adjacent, exactly as published.
3. **`395 aa`, not `385 aa`** — see §1.
4. **Zoom window.** The published inset framing was measured off the screenshot
   (exon-8 left edge at ~22 % of the panel width, exon-6 right edge at ~93 %) and
   reproduced as **chr10:80,273,890-80,276,800**. The start is deliberately 7 bp
   3' of the end of exon 9 (80,273,883) so that, as in the published inset, the
   transcripts enter the frame as dashed intron lines rather than as a clipped
   exon block.
5. **Panel-c variant markup** — Ensembl 116 no longer renders it for four of the
   five rows; see §3 and §7.
6. **Panel c 5'/3' coordinate labels.** Ensembl prints both the start (15481) and
   end (15600) coordinate of each row. The published panel is cropped after
   ~82 nt so only 15481 is visible. Crop as you prefer.

---

## 6. How to regenerate

Requirements: `python3` with `playwright` (chromium installed:
`python3 -m playwright install chromium`), and `rsvg-convert` for
SVG → PDF/PNG conversion.

### Panel b

```python
# 1. open Ensembl once and accept the cookie banner (sets the session cookie)
# 2. load the region WITH g=ENSG00000151224 — this is what makes Ensembl add the
#    "Genes (Comprehensive set from GENCODE)" track to the session
#      https://www.ensembl.org/Homo_sapiens/Location/View?r=<REGION>;g=ENSG00000151224
# 3. apply the track configuration.  Config URLs only take effect when fetched
#    as XHR, and the path is Config/Location/ViewBottom (NOT .../View):
#      fetch(url, {headers:{'X-Requested-With':'XMLHttpRequest'}, credentials:'same-origin'})
#    base = https://www.ensembl.org/Homo_sapiens/Config/Location/ViewBottom
#             ?db=core;g=ENSG00000151224;r=<REGION>;submit=1;
#    then, in order:  transcript_core_ensembl=transcript_label
#                     contig=off  seq=off  chr_band_core=off  age_of_base=off
#                     alignment_compara_9593_constrained=off
#                     dna_align_core_alt_seq_mapping=off  gencode_primary=off
#                     mane_select=off  nstd166=off  regulatory_build=off
#                     variation_set_ph_variants=off
#    (transcript_core_ensembl returns HTTP 500 if step 2 has not run first)
# 4. reload the region page, wait for img.imagemap
# 5. POST the ImageExport form shown in §2 -> SVG bytes
# 6. python3 extract_rows.py <exported>.svg <out>.svg
# 7. rsvg-convert -f pdf out.svg -o out.pdf ; rsvg-convert -w 6000 out.svg -o out.png
```

`REGION` = `10:80269780-80290126` (full) or `10:80273870-80276800` (zoom).

### Panel c

```
open https://www.ensembl.org/Homo_sapiens/Gene/TranscriptComparison?db=core;g=ENSG00000151224
   ;t1=ENST00000372213;t2=ENST00000871619;t3=ENST00000871620
   ;t4=ENST00000871621;t5=ENST00000871622
wait for pre.text_sequence, then scroll down to the 15481 block and wait for the
lazily-applied variant "adornment" markup, then clone the five rows + the
Variants key into a standalone element and print-to-PDF / screenshot it.
```

The working scripts are kept next to the assets in **`scripts_fig5bc/`**:

| Script | Purpose |
|---|---|
| `gen_b4.py` | drives Ensembl and exports the panel-b full-view SVG/PDF/PNG |
| `gen_b6.py` | same for the exon 7-9 zoom window (`10:80273890-80276800`) |
| `extract_rows.py` | pulls the four transcript rows out of the exported SVG and restacks them |
| `gen_c5.py` | captures the authentic Ensembl-116 Transcript-comparison block (version b) |
| `build_fig5c.py` | rebuilds panel c from Ensembl REST data (version a); writes HTML + SVG |
| `render_html.py` | renders an HTML file to vector PDF + 6× PNG with Playwright |

Paths inside the scripts point at the scratch directory they were written in;
edit the `SD` constant before re-running.

---

## 7. Verification against the current figure

Every generated file was rendered back and compared, element by element, with
`_reference_current/CURRENT_fig5_ALL.png`.

### Panel b — full view: **matches**

| Check | Result |
|---|---|
| Four transcripts, in the published top-to-bottom order 205 / 208 / 207 / 201 | ✔ |
| Transcript IDs `ENST00000871619`, `…871622`, `…871621`, `ENST00000372213` | ✔ identical |
| Label format `< MAT1A-2xx - ENSTxxxxxxxxxxx` / `protein coding` | ✔ identical |
| Colours: three brick-red (protein coding) + MAT1A-201 gold (MANE Select) | ✔ |
| MAT1A-205 extends further 3' (left) than the others; extra 20 bp exon 1 at the far right | ✔ |
| MAT1A-208 lacks the exon present in 205/207/201 at ~78 % of the width (8 vs 9/10 exons) | ✔ |
| Pale vertical grid lines | ✔ |
| Horizontal scale (gene span 80,269,780-80,290,126) | ✔ |

### Panel b — zoom: **matches**

| Check | Result |
|---|---|
| Three exons per transcript (exons 7-9 of MAT1A-205 = 6-8 of MAT1A-201) | ✔ |
| MAT1A-207's variant-containing exon visibly **shorter** (68 bp vs 134 bp) | ✔ |
| MAT1A-207 and MAT1A-201's 5'-most exon visibly **longer** (219 bp vs 144 bp) | ✔ |
| Dashed intron continuation at both frame edges | ✔ |
| Same relative exon positions as the published inset (±1 % of panel width) | ✔ |

### Panel c — primary (recomputed): **matches, with the caveats below**

| Check | Result |
|---|---|
| Five rows, order MAT1A-201 / -205 / -206 / -207 / -208 | ✔ |
| Row coordinate `15481` | ✔ identical |
| 120 nt per row (Ensembl default) | ✔ |
| Colour key: 14 chips, Ensembl's exact hex colours, wrapped 7 + 7 as published | ✔ |
| MAT1A-201/-205/-206/-208 identically coloured; MAT1A-207 distinct | ✔ |
| MAT1A-207 shown as intronic (dark blue `#02599c`) over the variant region | ✔ |
| Underlining of variant bases; grey introns; blue translated bases | ✔ |
| Window sequence identical to the live Ensembl page row | ✔ verified byte-for-byte |

Differences, all traceable to the release change:

1. The published row starts 24 nt further 3' (release-115 gene span), so the
   variant sits in **column 72** here versus column 48 in the published panel.
   The published crop corresponds to **columns 25-106** of the new row.
2. A handful of individual bases differ in colour because release 116 contains
   dbSNP variants that release 115 did not, and because VEP now calls
   `splice_polypyrimidine_tract_variant` for MAT1A-207 across the 3' end of the
   window (coral) where the published panel showed plain intron grey.
3. The 2 bp deletion `rs1397155707` is drawn over both bases here; Ensembl draws
   it over one.

### Panel c — secondary (Ensembl 116 live page): **does not match, by design**

Only MAT1A-201 carries variant markup (see §3). Supplied for auditability only;
do not use it as the figure.

### Not reproduced (must be re-added by hand)

Everything listed in §4 — the amino-acid labels, red dashed/dotted variant lines,
blue zoom box and connectors, red variant box and arrow, "SNP: rs118204006"
callouts and the "Transcript strand (−)" arrow. None of these come from Ensembl.
