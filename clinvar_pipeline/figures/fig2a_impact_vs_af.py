"""Fig 2a — MAGI impact score vs gnomAD AF (from fig2a_impact_vs_af.ipynb)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_auc_score

from config import PATHS
from figures.benchmark_data import load_impact_benchmark
from figures.common import DPI, ensure_dirs

LABEL_COLORS = {"Pathogenic": "#d73027", "Benign": "#4575b4"}
AF_RARE = 1e-4
AF_COMMON = 0.01
WIN_NAME = "Global_z_sum_log"
WIN_GROUP = "Global"


def _evaluate_score(values: np.ndarray, labels_bin: np.ndarray) -> dict:
    valid = np.isfinite(values)
    v, y = values[valid], labels_bin[valid]
    if valid.sum() < 30 or y.sum() < 5 or (1 - y).sum() < 5:
        return {"AUC": np.nan, "Cohen_d": np.nan}
    auc = max(roc_auc_score(y, v), 1.0 - roc_auc_score(y, v))
    p, b = v[y == 1], v[y == 0]
    pooled = np.sqrt(
        ((len(p) - 1) * p.std() ** 2 + (len(b) - 1) * b.std() ** 2)
        / max(len(p) + len(b) - 2, 1)
    )
    d = abs(p.mean() - b.mean()) / max(pooled, 1e-12)
    return {"AUC": float(auc), "Cohen_d": float(d)}


def _plot_impact_vs_af(df: pd.DataFrame, winner: dict, out_path: Path) -> None:
    score_col = WIN_NAME
    plot_df = df[
        df["gnomADe_AF"].notna()
        & (df["gnomADe_AF"] > 0)
        & df[score_col].notna()
        & np.isfinite(df[score_col].values)
    ].copy()
    plot_df["log10_AF"] = np.log10(plot_df["gnomADe_AF"])
    plot_df = plot_df[np.isfinite(plot_df["log10_AF"]) & np.isfinite(plot_df[score_col])]
    if plot_df.empty:
        raise ValueError("No variants with AF > 0 and finite impact score")

    impact_thresh = plot_df.loc[plot_df["label"] == "Pathogenic", score_col].median()

    fig, ax = plt.subplots(figsize=(10, 7.5))
    for lbl, zo, alpha in [("Benign", 2, 0.30), ("Pathogenic", 3, 0.45)]:
        sub = plot_df[plot_df["label"] == lbl]
        is_snp = sub["variant_type"].astype(str).str.upper() == "SNP"
        if is_snp.any():
            ax.scatter(
                sub.loc[is_snp, "gnomADe_AF"],
                sub.loc[is_snp, score_col],
                c=LABEL_COLORS[lbl],
                s=12,
                alpha=alpha,
                edgecolors="none",
                marker="o",
                zorder=zo,
                rasterized=True,
                label=f"{lbl} (SNP)",
            )
        if (~is_snp).any():
            ax.scatter(
                sub.loc[~is_snp, "gnomADe_AF"],
                sub.loc[~is_snp, score_col],
                c=LABEL_COLORS[lbl],
                s=16,
                alpha=alpha,
                edgecolors="none",
                marker="^",
                zorder=zo,
                rasterized=True,
                label=f"{lbl} (Indel)",
            )

    for lbl in ["Benign", "Pathogenic"]:
        sub = plot_df[plot_df["label"] == lbl]
        try:
            sns.kdeplot(
                x=sub["log10_AF"],
                y=sub[score_col],
                ax=ax,
                color=LABEL_COLORS[lbl],
                levels=4,
                linewidths=0.8,
                alpha=0.5,
                zorder=4,
            )
        except (ValueError, np.linalg.LinAlgError):
            pass

    ax.set_xscale("log")
    ax.set_xlabel("gnomAD Exome Allele Frequency", fontsize=13, fontweight="bold", labelpad=10)
    ax.set_ylabel(
        f"MAGI Continuous Impact Score ({WIN_GROUP})",
        fontsize=13,
        fontweight="bold",
        labelpad=10,
    )
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, _p: f"{x:.0e}" if x < 1e-3 else f"{x:.1%}")
    )
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=7, steps=[1, 2, 5, 10]))
    ax.grid(True, alpha=0.3, linestyle="-", linewidth=0.5)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.axhline(impact_thresh, color="#333333", ls=":", lw=1.2, alpha=0.6, zorder=1)
    ax.axvline(AF_COMMON, color="#333333", ls="--", lw=1.2, alpha=0.6, zorder=1)
    ax.fill_between([AF_COMMON, xlim[1]], impact_thresh, ylim[1], alpha=0.08, color="#ff8c00", zorder=0)
    ax.fill_between([xlim[0], AF_RARE], impact_thresh, ylim[1], alpha=0.08, color="#d73027", zorder=0)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.text(0.02, 0.97, "True\nPathogenic", transform=ax.transAxes, fontsize=11, fontweight="bold",
            fontstyle="italic", color="#d73027", ha="left", va="top", zorder=20)
    ax.text(0.98, 0.97, '"Exonerated"', transform=ax.transAxes, fontsize=11, fontweight="bold",
            fontstyle="italic", color="#cc6600", ha="right", va="top", zorder=20)
    ax.text(AF_COMMON * 1.3, ylim[1] * 0.60, "AF = 1%", fontsize=8, color="#777777",
            va="center", ha="left", fontstyle="italic", rotation=90, zorder=20)

    leg = ax.legend(loc="lower left", fontsize=9, framealpha=0.92, edgecolor="#cccccc",
                    markerscale=1.8, ncol=2, fancybox=True)
    leg.get_frame().set_linewidth(0.6)

    info = (
        f"n = {len(plot_df):,}  |  Score: {WIN_NAME}\n"
        f"AUC = {winner['AUC']:.3f}  |  Cohen's d = {winner['Cohen_d']:.2f}"
    )
    ax.text(0.98, 0.02, info, transform=ax.transAxes, fontsize=7.5, color="#888888",
            va="bottom", ha="right", family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#dddddd", alpha=0.90))
    ax.set_title(
        "MAGI Continuous Impact Score vs. Population Allele Frequency",
        fontsize=15,
        fontweight="bold",
        pad=16,
    )
    plt.subplots_adjust(left=0.1, right=0.95, top=0.93, bottom=0.11)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def run() -> None:
    ensure_dirs()
    df = load_impact_benchmark()
    if WIN_NAME not in df.columns:
        raise KeyError(f"{WIN_NAME} missing — run scripts/06_compute_impact_scores.py")

    labels_bin = (df["label"] == "Pathogenic").astype(int).values
    winner = _evaluate_score(df[WIN_NAME].values.astype(float), labels_bin)
    winner["Group"] = WIN_GROUP

    print(f"  Benchmark: {len(df):,} labelled variants")
    print(f"  Score: {WIN_NAME}  AUC={winner['AUC']:.3f}  d={winner['Cohen_d']:.3f}")

    out = Path(PATHS["fig2a"])
    _plot_impact_vs_af(df, winner, out)
    print(f"  Fig2a → {out}")
