#!/usr/bin/env python3
"""
Fig. 7a - GLA (ENSG00000102393) transcript-structure schematic.

Vector re-draw of the low-resolution Ensembl screenshot used in the MAGI/NAR
manuscript.  Exon coordinates come from the Ensembl REST API, release 109
(feb2023 archive), GRCh38 - the release whose transcript set and version
numbers match the published panel exactly (16 transcripts).

Two variants are produced:
  * default      - reproduces the *geometry* of the published screenshot
                   (piecewise-linear x-axis with the two/three intron
                   compressions present in the original), so the native
                   PowerPoint overlays (exon numbers, red variant line,
                   429 aa / RI / NMD labels) still line up.
  * --truescale  - strictly linear genomic x-axis, wider left margin,
                   larger labels.  Scientifically cleaner, but the
                   PowerPoint overlays must be repositioned.

Usage:  python3 make_fig7a.py [--truescale]
"""
import argparse
import json
import os
import sys
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_data_fig7a_ensembl109.json")

ENSEMBL_URL = ("https://feb2023.rest.ensembl.org/lookup/id/ENSG00000102393"
               "?expand=1;content-type=application/json")

# transcript rows, top-to-bottom, exactly as in the published panel
ORDER = [
    "ENST00000218516.4", "ENST00000466414.2", "ENST00000480513.6",
    "ENST00000468823.2", "ENST00000675968.1", "ENST00000479445.2",
    "ENST00000493905.6", "ENST00000675799.1", "ENST00000676156.1",
    "ENST00000674127.1", "ENST00000674142.1", "ENST00000674634.2",
    "ENST00000676372.1", "ENST00000675592.1", "ENST00000649178.1",
    "ENST00000486121.6",
]

# ---------------------------------------------------------------- style ----
EXON_FILL = "#EDEDED"
STROKE = "#565F64"
TEXT = "#000000"

# Reference-screenshot geometry (units = pixels of the original 1648x688 PNG).
# The PowerPoint picture frame crops it to x in [0, 1357.9], y in [7.36, 688].
REF_W, REF_H = 1357.9, 688.0
REF_Y0 = 7.36
BAND_TOP0, BAND_PITCH, BAND_H = 35.0, 39.8667, 33.0
LABEL_RIGHT = 202.0            # right edge of the transcript-ID label column
BOX_PAD = 0.0                   # boxes span the exon exactly (see PROVENANCE)

# Piecewise-linear genomic -> reference-pixel map, least-squares fitted to the
# 210 exon edges measured off the original screenshot.  Anchored at G0 to keep
# full numerical precision.  rms residual 0.33 px, max 2.58 px.
FIT = dict(G0=101397803,
           slope=0.14785300894216793,      # reference px per bp
           x0=216.43247735032767,          # reference x of G0
           breaks=[101402800, 101404500, 101406500],
           jumps=[-148.38998773698228, -18.864733392321828,
                  -218.41863848567635])


def fetch():
    if not os.path.exists(CACHE):
        sys.stderr.write("fetching %s\n" % ENSEMBL_URL)
        with urllib.request.urlopen(ENSEMBL_URL, timeout=120) as fh:
            data = json.load(fh)
        with open(CACHE, "w") as fh:
            json.dump(data, fh)
    with open(CACHE) as fh:
        return json.load(fh)


def use_arial():
    for path in ("/mnt/c/Windows/Fonts/arial.ttf", "/mnt/c/Windows/Fonts/arialbd.ttf"):
        if os.path.exists(path):
            fm.fontManager.addfont(path)
    names = {f.name for f in fm.fontManager.ttflist}
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = (["Arial"] if "Arial" in names else []) + \
        ["Liberation Sans", "DejaVu Sans"]
    plt.rcParams["svg.fonttype"] = "none"   # keep text as text in the SVG
    plt.rcParams["pdf.fonttype"] = 42       # embed TrueType (editable) in PDF


def xmap_ref(g):
    x = FIT["slope"] * (g - FIT["G0"]) + FIT["x0"]
    for b, j in zip(FIT["breaks"], FIT["jumps"]):
        if g > b:
            x += j
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truescale", action="store_true")
    args = ap.parse_args()

    gene = fetch()
    tx = {t["id"] + "." + str(t["version"]): t for t in gene["Transcript"]}
    missing = [t for t in ORDER if t not in tx]
    if missing:
        raise SystemExit("transcripts missing from Ensembl payload: %s" % missing)

    use_arial()

    rows = []
    for tid in ORDER:
        t = tx[tid]
        exons = sorted((e["start"], e["end"]) for e in t["Exon"])
        rows.append((tid, t["biotype"], exons))

    if args.truescale:
        gmin = min(e[0] for _, _, ex in rows for e in ex)
        gmax = max(e[1] for _, _, ex in rows for e in ex) + 1
        span = gmax - gmin
        label_w, plot_w, right_pad = 260.0, 1180.0, 20.0
        canvas_w = label_w + plot_w + right_pad
        pad = 0.0

        def X(g):
            return label_w + (g - gmin) / span * plot_w
        fig_w_mm = 180.0
    else:
        canvas_w = REF_W
        pad = BOX_PAD
        X = xmap_ref
        fig_w_mm = 180.0

    canvas_h = REF_H - REF_Y0
    fig_w_in = fig_w_mm / 25.4
    scale = fig_w_in / canvas_w                       # inches per reference px
    fig_h_in = canvas_h * scale
    px2pt = scale * 72.0                              # points per reference px

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, canvas_w)
    ax.set_ylim(REF_H, REF_Y0)
    ax.axis("off")
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    lw_box = 2.0 * px2pt
    lw_line = 1.8 * px2pt
    fs_label = 18.0 * px2pt

    for i, (tid, biotype, exons) in enumerate(rows):
        top = BAND_TOP0 + i * BAND_PITCH
        mid = top + BAND_H / 2.0

        # intron line: full transcript span, clipped to the canvas.  A
        # transcript whose 5'/3' end lies outside the window (only
        # ENST00000674142.1, whose distal exon is at chrX:101,393,273-101,393,829)
        # gets a stub running to the edge of the plotting area, exactly as in
        # the original panel.
        left_limit = 0.0 if args.truescale else LABEL_RIGHT + 6.0
        x0 = max(left_limit, X(exons[0][0]))
        x1 = min(canvas_w, X(exons[-1][1] + 1))
        ax.plot([x0, x1], [mid, mid], color=STROKE, lw=lw_line,
                solid_capstyle="butt", zorder=2)

        for (gs, ge) in exons:
            bx0, bx1 = X(gs) - pad, X(ge + 1) + pad
            if bx1 < 0 or bx0 > canvas_w:
                continue                              # off-canvas exon
            bx0, bx1 = max(bx0, 0.0), min(bx1, canvas_w)
            ax.add_patch(Rectangle((bx0, top), bx1 - bx0, BAND_H,
                                   facecolor=EXON_FILL, edgecolor=STROKE,
                                   lw=lw_box, joinstyle="miter", zorder=3))

        lx = LABEL_RIGHT if not args.truescale else label_w - 14.0
        ax.text(lx, mid, tid, ha="right", va="center",
                fontsize=fs_label, color=TEXT, zorder=5)

    stem = "fig7a_GLA_transcripts" + ("_trueScale" if args.truescale else "")
    for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 600})):
        out = os.path.join(HERE, stem + "." + ext)
        fig.savefig(out, facecolor="white", **kw)
        sys.stderr.write("wrote %s\n" % out)
    plt.close(fig)


if __name__ == "__main__":
    main()
