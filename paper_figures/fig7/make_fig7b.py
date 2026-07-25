#!/usr/bin/env python3
"""
Fig. 7b - multiple sequence alignment of human GLA / alpha-galactosidase A
protein isoforms (UniProtKB SwissProt P06280 + 9 TrEMBL entries).

Vector re-draw of the low-resolution UniProt "Align" screenshot used in the
MAGI/NAR manuscript.  The alignment itself is *not* re-derived by eye: the ten
FASTA sequences are pulled from the UniProt REST API and aligned with Clustal
Omega (EBI Job Dispatcher), which reproduces the published alignment
column-for-column (verified against every visible residue-count in the
original screenshot).

Layout, cell size, colour scheme and the red signal-peptide box are matched to
the original.  The blue Met267Ile arrow of the published figure is a native
PowerPoint overlay and is deliberately NOT drawn here; see PROVENANCE for the
exact coordinates at which to re-apply it.

Usage:
    python3 make_fig7b.py                 # complete alignment (9 blocks)
    python3 make_fig7b.py --as-published  # truncated at column 385 (7 blocks)
"""
import argparse
import os
import sys
import urllib.request

import matplotlib
import matplotlib.textpath
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA = os.path.join(HERE, "_data_fig7b_input.fasta")
ALN = os.path.join(HERE, "_data_fig7b_clustalo.aln")

# row order exactly as in the published panel
ACCS = ["A0A669KB83", "A0A3B3IRU3", "V9GYN5", "A0A6Q8PGG0", "A0A6Q8PFA9",
        "P06280", "A0AA34QW02", "A0A3B3IUC4", "A0A6Q8PHD1", "A0A6Q8PHM8"]
CANONICAL = "P06280"
SIGNAL_FROM, SIGNAL_TO = 1, 31          # P06280 signal peptide (UniProt)
SIGNAL_TRACK_LABEL = "A0A669KB83:Signal"
COLS_PER_BLOCK = 55

# --------------------------------------------------------------- colours ---
BG = "#FCFEFF"
SHADE = [(1.00, "#9293F9"),             # column fully identical in all 10 rows
         (0.80, "#A6A7FA"),
         (0.60, "#D5D8FC"),
         (0.40, "#EEF0F1")]
RESIDUE_TXT = "#101018"
GAP_TXT = "#3A3A3A"
NAME_TXT = "#3C3C3C"
ACC_TXT = "#31547D"
COUNT_TXT = "#444444"
RED = "#BB2739"
RED_BAR = "#C23D4D"

# ------------------------------------------- geometry (reference pixels) ---
# units are pixels of the original 912 x 1228 screenshot; the PowerPoint frame
# crops it to x >= 42.9, y >= 8.9
X_LEFT = 43.0
NAME_X = 47.0
NAME_GAP = 10.0                          # gap between name column and alignment
SEQ_X1 = 858.0                           # right edge of the alignment area
COUNT_X = 898.0
CANVAS_W = 912.0 - X_LEFT
Y_TOP = 8.9
BLOCK1_TOP = 12.0
ROW_H = 13.3
ROWS_H = 10 * ROW_H                      # 133.0
SIGNAL_GAP = 13.0                        # rows bottom -> feature bar top
SIGNAL_BAR_H = 10.0
GAP_AFTER_FEATURE = 22.5                 # inter-block gap, block with feature row
GAP_PLAIN = 32.0                         # inter-block gap, block without one
BOTTOM_MARGIN = 10.0

FS_RES_PER_CELL = 0.8926                 # residue type size as a fraction of cell width
FS_NAME = 11.0                           # font sizes, reference px
FS_COUNT = 11.0
FS_SIGNAL = 11.0

FIG_W_MM = 180.0                         # nominal size; vector output scales freely


def use_arial():
    for p in ("/mnt/c/Windows/Fonts/arial.ttf", "/mnt/c/Windows/Fonts/arialbd.ttf"):
        if os.path.exists(p):
            fm.fontManager.addfont(p)
    names = {f.name for f in fm.fontManager.ttflist}
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = (["Arial"] if "Arial" in names else []) + \
        ["Liberation Sans", "DejaVu Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42


def fetch_sequences():
    if os.path.exists(FASTA):
        return
    chunks = []
    for acc in ACCS:
        url = "https://rest.uniprot.org/uniprotkb/%s.fasta" % acc
        sys.stderr.write("fetching %s\n" % url)
        with urllib.request.urlopen(url, timeout=120) as fh:
            txt = fh.read().decode()
        lines = txt.strip().split("\n")
        chunks.append(lines[0].split()[0] + "\n" + "".join(l.strip() for l in lines[1:]))
    with open(FASTA, "w") as fh:
        fh.write("\n".join(chunks) + "\n")


def read_alignment():
    """Parse the cached Clustal Omega output (clustal_num format)."""
    if not os.path.exists(ALN):
        raise SystemExit(
            "%s not found.  Regenerate it with:\n"
            "  curl -X POST --form email=YOU@EXAMPLE.ORG "
            "--form sequence=\"$(cat %s)\" --form outfmt=clustal_num "
            "--form order=input "
            "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run\n"
            "  curl .../result/<jobid>/aln-clustal_num -o %s" % (ALN, FASTA, ALN))
    seqs, order = {}, []
    for line in open(ALN):
        if line.startswith("CLUSTAL") or not line.strip() or line.startswith(" "):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name, block = parts[0], parts[1]
        if set(block) <= set("*:. "):
            continue
        if name not in seqs:
            seqs[name] = []
            order.append(name)
        seqs[name].append(block)
    return order, {k: "".join(v) for k, v in seqs.items()}


def shade_for(column):
    """Background colour for an alignment column (identity over all rows)."""
    top = max(column.count(c) for c in set(column) if c != "-") if any(
        c != "-" for c in column) else 0
    frac = top / float(len(column))
    for thr, col in SHADE:
        if frac >= thr - 1e-9:
            return col
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--as-published", action="store_true",
                    help="truncate at alignment column 385 like the original crop")
    args = ap.parse_args()

    fetch_sequences()
    order, seqs = read_alignment()
    if len(order) != len(ACCS):
        raise SystemExit("expected %d sequences, got %d" % (len(ACCS), len(order)))
    # re-order rows to the published order, keyed on accession
    by_acc = {n.split("|")[1]: n for n in order}
    order = [by_acc[a] for a in ACCS]
    aln_len = len(seqs[order[0]])
    total_cols = 385 if args.as_published else aln_len
    n_blocks = (total_cols + COLS_PER_BLOCK - 1) // COLS_PER_BLOCK

    columns = ["".join(seqs[n][i] for n in order) for i in range(aln_len)]
    shades = [shade_for(c) for c in columns]

    use_arial()
    fig_w_in = FIG_W_MM / 25.4
    scale = fig_w_in / CANVAS_W
    px2pt = scale * 72.0

    # renderer-free text advance, in reference pixels (cached)
    _t2p = matplotlib.textpath.TextToPath()
    _cache = {}

    def adv(txt, weight="normal", size=None):
        size = FS_NAME if size is None else size
        key = (txt, weight, size)
        if key not in _cache:
            prop = fm.FontProperties(
                family=plt.rcParams["font.sans-serif"][0],
                size=size * px2pt, weight=weight)
            w = _t2p.get_text_width_height_descent(txt, prop, False)[0]
            _cache[key] = w / px2pt          # points -> reference px
        return _cache[key]

    # Row labels are bare accessions, so the name column is sized to the widest
    # accession rather than to a fixed 145 px; the width this frees goes to the
    # alignment cells, which enlarges the residue type by the same proportion.
    name_w = max(adv(a, "bold", FS_NAME) for a in ACCS)
    seq_x0 = NAME_X + name_w + NAME_GAP
    cell_w = (SEQ_X1 - seq_x0) / COLS_PER_BLOCK
    fs_res = cell_w * FS_RES_PER_CELL

    # Vertical layout: the feature row is emitted only for blocks whose column
    # range actually intersects the annotated signal peptide, and blocks that do
    # not carry one reserve no space for it.
    def has_feature(c0, c1):
        return c1 > SIGNAL_FROM - 1 and c0 < SIGNAL_TO

    layout, y = [], BLOCK1_TOP
    for b in range(n_blocks):
        c0 = b * COLS_PER_BLOCK
        c1 = min(c0 + COLS_PER_BLOCK, total_cols)
        feat = has_feature(c0, c1)
        layout.append((c0, c1, y, feat))
        y += ROWS_H + (SIGNAL_GAP + SIGNAL_BAR_H if feat else 0.0)
        if b < n_blocks - 1:
            y += GAP_AFTER_FEATURE if feat else GAP_PLAIN
    canvas_h = y + BOTTOM_MARGIN - Y_TOP

    fig = plt.figure(figsize=(fig_w_in, canvas_h * scale))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(X_LEFT, 912.0)
    ax.set_ylim(Y_TOP + canvas_h, Y_TOP)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.add_patch(Rectangle((X_LEFT, Y_TOP), CANVAS_W, canvas_h,
                           facecolor=BG, edgecolor="none", zorder=0))

    for (c0, c1, btop, feat) in layout:
        for r, name in enumerate(order):
            rtop = btop + r * ROW_H
            seq = seqs[name]
            # merged shading runs (avoids hairline seams in PDF/SVG)
            run_start, run_col = None, None
            for c in range(c0, c1 + 1):
                col = shades[c] if (c < c1 and seq[c] != "-") else None
                if col != run_col:
                    if run_col is not None:
                        x = seq_x0 + (run_start - c0) * cell_w
                        ax.add_patch(Rectangle(
                            (x, rtop), (c - run_start) * cell_w, ROW_H,
                            facecolor=run_col, edgecolor="none", zorder=1))
                    run_start, run_col = c, col
            for c in range(c0, c1):
                ch = seq[c]
                ax.text(seq_x0 + (c - c0 + 0.5) * cell_w, rtop + ROW_H / 2.0, ch,
                        ha="center", va="center_baseline",
                        fontsize=fs_res * px2pt,
                        color=RESIDUE_TXT if ch != "-" else GAP_TXT, zorder=3)

            # row label: bare UniProt accession; the SwissProt canonical entry is
            # set bold and its residue counter in red
            acc = name.split("|")[1]
            canon = acc == CANONICAL
            y = rtop + ROW_H / 2.0
            ax.text(NAME_X, y, acc, ha="left", va="center_baseline",
                    fontsize=FS_NAME * px2pt, color=ACC_TXT,
                    fontweight="bold" if canon else "normal", zorder=3)

            cnt = len(seq[:c1].replace("-", ""))
            ax.text(COUNT_X, y, str(cnt), ha="right", va="center_baseline",
                    fontsize=FS_COUNT * px2pt,
                    color=RED if canon else COUNT_TXT,
                    fontweight="bold" if canon else "normal", zorder=3)

        if not feat:
            continue
        # signal-peptide feature track (only under blocks that contain it)
        sy = btop + ROWS_H + SIGNAL_GAP
        s0, s1 = max(SIGNAL_FROM - 1, c0), min(SIGNAL_TO, c1)
        bar_x0 = seq_x0 + (s0 - c0) * cell_w
        bar_w = (s1 - s0) * cell_w
        ax.add_patch(Rectangle((bar_x0, sy), bar_w, SIGNAL_BAR_H,
                               facecolor=RED_BAR, edgecolor="none", zorder=3))
        # red box around the signal peptide across all rows
        ax.add_patch(Rectangle((bar_x0, btop), bar_w, ROWS_H,
                               facecolor="none", edgecolor=RED,
                               lw=2.0 * px2pt, joinstyle="miter", zorder=4))
        # label sits beside the bar: the narrower accession-only gutter can no
        # longer hold it at a legible size
        ax.text(bar_x0 + bar_w + 8.0, sy + SIGNAL_BAR_H / 2.0,
                SIGNAL_TRACK_LABEL, ha="left", va="center_baseline",
                fontsize=FS_SIGNAL * px2pt, color="#111111",
                fontweight="bold", zorder=3)

    stem = "fig7b_GLA_MSA" + ("_asPublished" if args.as_published else "")
    for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 600})):
        out = os.path.join(HERE, stem + "." + ext)
        fig.savefig(out, facecolor="white", **kw)
        sys.stderr.write("wrote %s\n" % out)
    plt.close(fig)


if __name__ == "__main__":
    main()
