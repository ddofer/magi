#!/usr/bin/env python3
"""
Regenerate Figure 2 panels d and e (per-element disruption) at print quality.

Panel d: "Per-Element Disruption Magnitude by Variant Class"
         mean |ΔPred| per BED element, Pathogenic vs Benign, 95% CI
Panel e: "Per-Element Disruption: Benign vs Pathogenic (95% CI)"
         mean |ΔPred| difference (Benign - Pathogenic), Bonferroni-significant in green

Statistics are copied verbatim from the cell that produced the published
figure (analyses/genomic_featImp_experiment_v3.ipynb, cell 24) so the numbers
are identical; only the output device changes (vector instead of a 245 dpi
screen-resolution PNG crop).

------------------------------------------------------------------------------
INPUT SUBSET - data archaeology note (2026-07, NAR revision)
------------------------------------------------------------------------------
The notebook that produced the published panels sets, in cell 5:

    MODEL_NAME         = "InstaDeepAI/NTv3_100M_post"
    CONTEXT_LEN        = 2048
    OUTPUT_RESULTS_FILE = "clinvar_full_deltas.parquet"

and cell 8's stored output confirms it loaded 40,976 variants / 14,771 columns
(7,362 D_BW_* + 21 D_BED_*) with labels {Pathogenic: 22252, Benign: 18724}.
Cell 24 then runs on ALL of those rows - there is NO quality/gold-star/SNP
subsetting anywhere between cell 8 and cell 24 (cell 11 only copies final_df
and appends an Impact_Score column).

`clinvar_full_deltas.parquet` no longer exists anywhere on disk.  The file that
does exist, `clinvar_new_deltas.parquet`, is the SAME 40,976 input variants
(same 22,252 / 18,724 label split, produced from the same clinvar_input.parquet)
re-scored with a DIFFERENT model: inference.py now uses
`InstaDeepAI/NTv3_650M_post` with CONTEXT_LEN = 64*1024, and its output carries
only 3,602 D_BW_* tracks instead of 7,362.  The two runs are therefore not
related by row selection:

    old (published)  D_BED_* range [-0.961818, +0.946396], mean|d| 0.010791
    new (on disk)    D_BED_* range [-0.774623, +0.714373], mean|d| 0.013053

A subset can only shrink a range, never widen it, so no filter over
clinvar_new_deltas.parquet can reproduce the published numbers.

AUTHOR DECISION (2026-07): the paper standardises on NTv3-650M, so this script
uses `clinvar_new_deltas.parquet` (650M / 64 kb) BY DEFAULT and never silently
switches to `clinvar_full_deltas.parquet` even if that file is present. The
resulting bars are a median ~1.4x taller than the printed figure, and the
element ordering shifts from position 5 onward, so any Results text quoting
Fig 2d/2e values must be re-checked. The analysis population is the SAME
unfiltered 40,976 variants the notebook used; the only row filter is the drop
of inference-failure rows (see DROP_ALL_NAN_BED below). Full audit in
fig2/fig2de_VALIDATION.md.

To generate the ORIGINAL published (100M / 2 kb) numbers instead, pass that file
explicitly:
    python regen_fig2_de_panels.py --deltas /path/to/clinvar_full_deltas.parquet

Usage:
    python regen_fig2_de_panels.py [--outdir DIR] [--deltas FILE]
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats


def configure_print_fonts() -> str:
    """Journal-production settings for vector output.

    matplotlib defaults to Type 3 fonts in PDF and outlines text in SVG. Type 3
    is routinely rejected or queried by journal production systems, and outlined
    text cannot be edited by the typesetter. Force TrueType (type 42) embedding
    and live SVG text, and prefer Arial to match NAR figure style.

    Arial is looked up in the usual Linux locations plus the Windows font
    directory (this repo is often used from WSL). Falls back to DejaVu Sans,
    which is still embedded as TrueType - the Type 3 problem is fixed either way.
    """
    from matplotlib import font_manager

    for path in (
        "/mnt/c/Windows/Fonts/arial.ttf", "/mnt/c/Windows/Fonts/arialbd.ttf",
        "/mnt/c/Windows/Fonts/ariali.ttf", "/mnt/c/Windows/Fonts/arialbi.ttf",
        "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
        "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf",
        "/Library/Fonts/Arial.ttf",
    ):
        if Path(path).exists():
            try:
                font_manager.fontManager.addfont(path)
            except Exception:
                pass

    available = {f.name for f in font_manager.fontManager.ttflist}
    family = next((f for f in ("Arial", "Liberation Sans", "Helvetica")
                   if f in available), "DejaVu Sans")

    matplotlib.rcParams.update({
        "pdf.fonttype": 42,      # TrueType, not Type 3
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # keep <text> live instead of outlining it
        "font.family": "sans-serif",
        "font.sans-serif": [family, "DejaVu Sans"],
    })
    return family


MAGI = Path(os.environ.get("MAGI_ROOT", Path(__file__).resolve().parents[1]))

# Delta table used by default.
#
# AUTHOR DECISION (2026-07): the paper standardises on NTv3-650M, so these
# panels are generated from the 650M / 64 kb re-run, NOT from the NTv3-100M /
# 2 kb file the published panels came from.  Bar heights therefore differ from
# the printed figure by a median factor of ~1.4x - see fig2de_VALIDATION.md.
#
# This is deliberately NOT a "use the published file if you find it" fallback:
# on a machine that still has clinvar_full_deltas.parquet, silently preferring
# it would produce 100M figures when 650M was asked for.  If you ever do want
# the published 100M numbers, pass them explicitly:
#     --deltas /path/to/clinvar_full_deltas.parquet
DEFAULT_DELTAS = MAGI / "clinvar_new_deltas.parquet"   # NTv3-650M, 64 kb
PUBLISHED_SOURCE = MAGI / "clinvar_full_deltas.parquet"  # NTv3-100M, 2 kb

# Row filter #1 (and the only one).  `clinvar_new_deltas.parquet` contains 4
# rows - chr9:14889 A>G, chr9:15126 G>C, chr9:15978 C>T, chr9:16020 A>G, all
# Benign - whose 21 D_BED_* values are ALL NaN, i.e. the inference call failed
# for that batch.  The published run had no NaNs (its per-element Welch t-tests
# all returned finite p-values), so these rows must be removed: scipy's
# ttest_ind propagates NaN and would otherwise return NaN for every element,
# making 0/21 rather than 21/21 elements Bonferroni-significant.
DROP_ALL_NAN_BED = True

# Published panel geometry, read off the embedded PPTX image (1248x1118 px,
# placed at 5.10 in wide => ~245 dpi effective). Panels are stacked in one
# source image; the deck crops the top half for (d) and the bottom for (e).
FAITHFUL_FIGSIZE = (14, 6)          # exactly what the notebook used
PRINT_TARGETS = [                    # (name, width_in, height_in, base_pt)
    ("print90mm", 3.543, 2.10, 5.0),
    ("print180mm", 7.087, 3.30, 7.0),
]


def compute_per_element_stats(df: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """Verbatim reimplementation of the notebook's per-element comparison."""
    delta_bed_cols = [c for c in df.columns if c.startswith("D_BED_")]
    path_df = df[df["label"] == "Pathogenic"]
    ben_df = df[df["label"] == "Benign"]

    results = []
    for col in delta_bed_cols:
        element = col.replace("D_BED_", "")
        p_vals = path_df[col].abs()
        b_vals = ben_df[col].abs()

        p_mean, b_mean = p_vals.mean(), b_vals.mean()
        p_sem, b_sem = p_vals.sem(), b_vals.sem()
        p_std, b_std = p_vals.std(), b_vals.std()

        p_ci = stats.t.interval(0.95, len(p_vals) - 1, loc=p_mean, scale=p_sem)
        b_ci = stats.t.interval(0.95, len(b_vals) - 1, loc=b_mean, scale=b_sem)
        _, p_value = stats.ttest_ind(p_vals, b_vals, equal_var=False)

        pooled_std = np.sqrt((p_std**2 + b_std**2) / 2)
        cohens_d = (b_mean - p_mean) / pooled_std if pooled_std > 0 else 0
        diff = b_mean - p_mean
        diff_se = np.sqrt(p_sem**2 + b_sem**2)

        results.append({
            "Element": element,
            "Pathogenic_Mean": p_mean, "Pathogenic_SEM": p_sem,
            "Pathogenic_CI_Lower": p_ci[0], "Pathogenic_CI_Upper": p_ci[1],
            "Benign_Mean": b_mean, "Benign_SEM": b_sem,
            "Benign_CI_Lower": b_ci[0], "Benign_CI_Upper": b_ci[1],
            "Difference_Benign_minus_Path": diff,
            "Diff_CI_Lower": diff - 1.96 * diff_se,
            "Diff_CI_Upper": diff + 1.96 * diff_se,
            "Cohens_d": cohens_d, "P_Value": p_value,
            "N_Pathogenic": len(p_vals), "N_Benign": len(b_vals),
        })

    res = pd.DataFrame(results).sort_values(
        "Difference_Benign_minus_Path", ascending=False
    )
    bonf = 0.05 / len(res)
    res["Significant"] = res["P_Value"] < bonf
    return res, bonf


def plot_panel_d(res: pd.DataFrame, figsize, s: float):
    """Grouped Pathogenic/Benign bars with 95% CI. *s* scales type/strokes."""
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(res))
    w = 0.35
    for offset, key, sem_key, label, colour in (
        (-w / 2, "Pathogenic_Mean", "Pathogenic_SEM", "Pathogenic", "#d62728"),
        (+w / 2, "Benign_Mean", "Benign_SEM", "Benign", "#1f77b4"),
    ):
        ax.bar(x + offset, res[key], w, label=label, color=colour, alpha=0.75,
               edgecolor="black", linewidth=0.5 * s,
               yerr=1.96 * res[sem_key], capsize=4 * s,
               error_kw={"linewidth": 1.0 * s})

    ax.set_xticks(x)
    ax.set_xticklabels(res["Element"], rotation=45, ha="right", fontsize=10 * s)
    ax.set_ylabel("Mean |ΔPred| (with 95% CI)", fontsize=11 * s, fontweight="bold")
    ax.set_title("Per-Element Disruption Magnitude by Variant Class",
                 fontsize=12 * s, fontweight="bold")
    ax.legend(fontsize=10 * s, loc="upper right")
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.8 * s)
    ax.tick_params(axis="y", labelsize=10 * s)
    fig.tight_layout()
    return fig


def plot_panel_e(res: pd.DataFrame, n_tests: int, figsize, s: float):
    """Benign-minus-Pathogenic difference with 95% CI; green = Bonferroni-sig."""
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(res))
    colours = ["green" if sig else "lightgray" for sig in res["Significant"]]
    errors = np.array([
        res["Difference_Benign_minus_Path"].values - res["Diff_CI_Lower"].values,
        res["Diff_CI_Upper"].values - res["Difference_Benign_minus_Path"].values,
    ])
    ax.bar(x, res["Difference_Benign_minus_Path"], color=colours, alpha=0.75,
           edgecolor="black", linewidth=0.5 * s, yerr=errors,
           capsize=5 * s, error_kw={"linewidth": 1.5 * s})
    ax.axhline(0, color="black", linewidth=1.5 * s, linestyle="-", alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(res["Element"], rotation=45, ha="right", fontsize=10 * s)
    ax.set_ylabel("Mean |ΔPred| Difference (Benign - Pathogenic)",
                  fontsize=11 * s, fontweight="bold")
    ax.set_title(
        "Per-Element Disruption: Benign vs Pathogenic (95% CI)\n"
        f"Green = Significant (Bonferroni-corrected), n={n_tests} elements",
        fontsize=12 * s, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.8 * s)
    ax.tick_params(axis="y", labelsize=10 * s)
    fig.tight_layout()
    return fig


def save_all(fig, base: Path, dpi: int = 600) -> None:
    for ext in ("pdf", "svg", "png"):
        kw = {"bbox_inches": "tight", "facecolor": "white"}
        if ext == "png":
            kw["dpi"] = dpi
        fig.savefig(base.with_suffix(f".{ext}"), **kw)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path,
                    default=Path(__file__).resolve().parent / "fig2",
                    help="output directory (default: <script dir>/fig2)")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--deltas", type=Path, default=None,
                    help=f"delta parquet to use. Default: {DEFAULT_DELTAS.name} "
                         "(NTv3-650M / 64 kb). Pass clinvar_full_deltas.parquet "
                         "explicitly if you want the published NTv3-100M / 2 kb "
                         "numbers instead.")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    font_family = configure_print_fonts()
    print(f"Vector output: Type-42 fonts, live SVG text, family={font_family}")

    deltas = args.deltas if args.deltas is not None else DEFAULT_DELTAS
    if not deltas.exists():
        raise SystemExit(f"delta parquet not found: {deltas}")

    is_published_source = deltas.name == PUBLISHED_SOURCE.name
    print("=" * 72)
    if is_published_source:
        print(f"SOURCE: {deltas.name}")
        print("        NTv3-100M @ 2 kb context - the ORIGINAL published source.")
        print("        Output should reproduce the printed Fig 2d/2e numbers.")
    else:
        print(f"SOURCE: {deltas.name}")
        print("        NTv3-650M @ 64 kb context - matches the rest of the paper.")
        print("        Bar heights will NOT match the printed Fig 2d/2e (median")
        print("        ~1.4x larger). This is the intended behaviour; see")
        print("        fig2de_VALIDATION.md. Any Results text quoting Fig 2d/2e")
        print("        values must be re-checked against the new numbers.")
        if PUBLISHED_SOURCE.exists():
            print()
            print(f"NOTE:   {PUBLISHED_SOURCE.name} also exists on this machine but is")
            print("        NOT being used. To generate the published 100M numbers:")
            print(f"            --deltas {PUBLISHED_SOURCE}")
    print("=" * 72)

    import pyarrow.parquet as pq
    names = pq.ParquetFile(deltas).schema.names
    bed_cols = [c for c in names if c.startswith("D_BED_")]
    cols = bed_cols + ["label"]
    print(f"Reading {len(cols)} columns from {deltas.name} ...")
    df = pd.read_parquet(deltas, columns=cols)
    # The results parquet stores `label` as a boolean; the notebook maps it the
    # same way before plotting (True -> Pathogenic, False -> Benign).
    if df["label"].dtype == bool:
        df["label"] = df["label"].map({True: "Pathogenic", False: "Benign"})
    # Analysis population = every labelled variant, exactly as in the notebook
    # (no gold-star / SNP-only / rationale subsetting is applied there).
    df = df[df["label"].isin(["Pathogenic", "Benign"])].copy()
    n_before = len(df)

    n_dropped = 0
    if DROP_ALL_NAN_BED:
        failed = df[bed_cols].isna().all(axis=1)
        n_dropped = int(failed.sum())
        df = df[~failed].copy()
        print(f"  dropped {n_dropped} inference-failure row(s) "
              f"(all {len(bed_cols)} D_BED_* NaN)")
    assert not df[bed_cols].isna().any().any(), "residual NaNs in D_BED_* columns"

    print(f"  {len(df):,} variants (of {n_before:,}) | "
          f"labels: {df['label'].value_counts().to_dict()}")

    res, bonf = compute_per_element_stats(df)
    res.to_csv(args.outdir / "fig2de_per_element_statistics.csv", index=False)

    n_p = int(res["N_Pathogenic"].iloc[0])
    n_b = int(res["N_Benign"].iloc[0])
    print(f"  N pathogenic={n_p:,}  N benign={n_b:,}  Bonferroni alpha={bonf:.5f}")
    print(f"  significant elements: {int(res['Significant'].sum())}/{len(res)}")
    print("\n  element order (as plotted, left to right):")
    print("   " + ", ".join(res["Element"]))

    for name, fn in (("fig2d_per_element_by_variant_class", plot_panel_d),
                     ("fig2e_per_element_benign_vs_pathogenic",
                      lambda r, fs, s: plot_panel_e(r, len(res), fs, s))):
        fig = fn(res, FAITHFUL_FIGSIZE, 1.0)
        save_all(fig, args.outdir / f"{name}_faithful", args.dpi)
        plt.close(fig)
        for label, w, h, base_pt in PRINT_TARGETS:
            fig = fn(res, (w, h), base_pt / 10.0)
            save_all(fig, args.outdir / f"{name}_{label}", args.dpi)
            plt.close(fig)
        print(f"  wrote {name}_(faithful|print90mm|print180mm).(pdf|svg|png)")

    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_data": str(deltas),
        "source_data_is_published_source": is_published_source,
        "published_source_expected": str(PUBLISHED_SOURCE),
        "published_source_present": PUBLISHED_SOURCE.exists(),
        "published_source_model": "InstaDeepAI/NTv3_100M_post, CONTEXT_LEN=2048",
        "default_source_model": "InstaDeepAI/NTv3_650M_post, CONTEXT_LEN=65536",
        "font_family": font_family,
        "row_filter": ("all labelled variants (no gold-star / SNP-only / "
                       "rationale subsetting, matching the notebook); rows with "
                       "all-NaN D_BED_* dropped"),
        "n_rows_before_filter": int(n_before),
        "n_rows_dropped_all_nan_bed": n_dropped,
        "source_code": "MAGI/analyses/genomic_featImp_experiment_v3.ipynb cell 24",
        "n_variants_total": int(len(df)),
        "n_pathogenic": n_p,
        "n_benign": n_b,
        "n_elements": int(len(res)),
        "bonferroni_alpha": bonf,
        "n_significant": int(res["Significant"].sum()),
        "element_order_left_to_right": list(res["Element"]),
    }
    (args.outdir / "fig2de_provenance.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8")
    print(f"\nWrote {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
