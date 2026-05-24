"""Plot helpers ported from parent fig2/fig3 notebooks."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde
from sklearn.metrics import roc_auc_score

from config import IMPACT_SCORE_COLS
from figures.common import DPI, normalize_label_series

SCORE_COLS = list(IMPACT_SCORE_COLS)
DEFAULT_QUARTILE_ORDER = ["Q1", "Q2", "Q3", "Q4"]

CONC_COLORS = {
    "CONCORDANT": "forestgreen",
    "PARTIAL": "goldenrod",
    "DISCORDANT": "firebrick",
}

MECH_COLOR_MAP = {
    "DISCORDANT": "#e25b45",
    "PARTIAL": "#f0a830",
    "CONCORDANT": "#6abf69",
    "NOT_APPLICABLE": "#bbbbbb",
}


def _safe_crosstab(row_vals, col_vals, row_order, col_order) -> pd.DataFrame:
    ct = pd.crosstab(row_vals, col_vals)
    return ct.reindex(index=row_order, columns=col_order, fill_value=0)


def _make_quartile_groups(series: pd.Series, labels=DEFAULT_QUARTILE_ORDER) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    out = pd.Series(index=s.index, dtype="object")
    valid = s.notna()
    if valid.sum() == 0:
        return out
    try:
        out.loc[valid] = pd.qcut(s.loc[valid], q=4, labels=labels, duplicates="drop").astype(str)
    except ValueError:
        ranked = s.loc[valid].rank(method="average", pct=True)
        out.loc[valid] = pd.cut(
            ranked,
            bins=[0, 0.25, 0.5, 0.75, 1.0],
            labels=labels,
            include_lowest=True,
        ).astype(str)
    return out


def plot_overall_concordance_by_variant_type(
    merged: pd.DataFrame,
    out_path: Path,
    *,
    variant_order: list[str] | None = None,
) -> None:
    """Fig 3a — CONCORDANT (solid) + PARTIAL (hatched) by SNP vs INDEL."""
    variant_order = variant_order or ["snp", "indel"]
    d = merged.copy()
    d["concordance"] = d["concordance"].astype(str).str.strip().str.upper()
    d["variant_type"] = d["variant_type"].astype(str).str.strip().str.lower()

    counts = (
        d.groupby(["variant_type", "concordance"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=variant_order)
    )
    frac = counts.div(counts.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(5, 5))
    palette = sns.color_palette("deep")
    colors = {"snp": palette[0], "indel": palette[1]}

    for i, vt in enumerate(variant_order):
        if vt not in frac.index:
            continue
        c = float(frac.loc[vt].get("CONCORDANT", 0))
        p = float(frac.loc[vt].get("PARTIAL", 0))
        ax.bar(i, c, color=colors.get(vt, "grey"), edgecolor="black")
        ax.bar(i, p, bottom=c, color=colors.get(vt, "grey"), edgecolor="black", hatch="//")

    cp_counts = counts[["CONCORDANT", "PARTIAL"]].sum(axis=1)
    totals = counts.sum(axis=1)
    labels = [
        f"{vt.upper()}\n({int(cp_counts.get(vt, 0))}/{int(totals.get(vt, 0))})"
        for vt in variant_order
    ]
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1)
    ticks = np.arange(0, 1.01, 0.1)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{int(t * 100)}%" for t in ticks])
    ax.set_ylabel("Fraction of variants")
    ax.set_title("Concordance by Variant Type")
    ax.legend(
        handles=[
            Patch(facecolor="white", edgecolor="black", label="CONCORDANT"),
            Patch(facecolor="white", edgecolor="black", hatch="//", label="PARTIAL"),
        ]
    )
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_concordance_bars(
    df: pd.DataFrame,
    group_col: str = "signal_category",
    concord_col: str = "concordance",
    group_order=None,
    conc_order=("CONCORDANT", "PARTIAL", "DISCORDANT"),
    conc_colors=None,
    normalize: bool = True,
    figsize=(12, 5),
    title: str | None = None,
    ax=None,
    score_cols=SCORE_COLS,
    quartile_labels=DEFAULT_QUARTILE_ORDER,
    use_quantile_binning_for_scores: bool = True,
):
    """Fig 3b/c — dual stacked bars (concordant+partial | discordant) by quartile."""
    d = df.copy()
    d[concord_col] = d[concord_col].astype(str).str.strip().str.upper()
    if conc_colors is None:
        conc_colors = CONC_COLORS

    is_score_mode = group_col in score_cols and use_quantile_binning_for_scores
    if is_score_mode:
        group_label = f"{group_col}_quartile"
        d[group_label] = _make_quartile_groups(d[group_col], labels=quartile_labels)
        d[group_label] = d[group_label].astype(str).str.strip().str.upper()
        group_order_final = [q.upper() for q in quartile_labels] if group_order is None else group_order
        plot_group_col = group_label
        xlabel = f"{group_col} quartile"
    else:
        d[group_col] = d[group_col].astype(str).str.strip().str.upper()
        plot_group_col = group_col
        group_order_final = list(pd.Index(d[group_col].dropna().unique())) if group_order is None else group_order
        xlabel = group_col

    conc_order = [c.upper() for c in conc_order]
    ct = _safe_crosstab(d[plot_group_col], d[concord_col], group_order_final, conc_order)

    if normalize:
        plot_df = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        ylabel = "Fraction"
    else:
        plot_df = ct.copy()
        ylabel = "Count"

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    x = np.arange(len(group_order_final))
    width = 0.35
    stacked = [c for c in conc_order if c != "DISCORDANT"]
    separate = [c for c in conc_order if c == "DISCORDANT"]

    bottom = np.zeros(len(group_order_final))
    for conc in stacked:
        if conc not in plot_df.columns:
            continue
        vals = plot_df[conc].values
        ax.bar(
            x - width / 2 - 0.02,
            vals,
            width,
            bottom=bottom,
            label=conc,
            color=conc_colors.get(conc, "#999999"),
            edgecolor="white",
            linewidth=0.5,
        )
        bottom += vals

    for conc in separate:
        if conc not in plot_df.columns:
            continue
        vals = plot_df[conc].values
        ax.bar(
            x + width / 3 + 0.04,
            vals,
            width,
            label=conc,
            color=conc_colors.get(conc, "#999999"),
            edgecolor="white",
            linewidth=0.5,
        )

    row_totals = ct.sum(axis=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n(n={int(row_totals.get(g, 0))})" for g in group_order_final])
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.legend(title=concord_col, fontsize=8, loc="upper right")
    if title:
        ax.set_title(title)
    return fig, ax, plot_df


def plot_concordance_by_mechanism(
    df: pd.DataFrame,
    vt_type: str,
    out_png: Path,
    out_csv_prefix: Path,
    *,
    mechanism_col: str = "NT_Mechanism",
    concordance_col: str = "concordance",
    min_samples: int = 20,
) -> None:
    """Fig 3e — horizontal stacked bars per NT mechanism (parent notebook)."""
    d = df.copy()
    d[concordance_col] = d[concordance_col].astype(str).str.strip().str.upper()

    title = (
        "SNP: Concordance by signal-assigned mechanism"
        if vt_type == "snp"
        else "INDEL: Concordance by signal-assigned mechanism"
    )

    ct = pd.crosstab(d[mechanism_col], d[concordance_col])
    preferred = ["DISCORDANT", "PARTIAL", "NOT_APPLICABLE", "CONCORDANT"]
    concordance_order = [c for c in preferred if c in ct.columns]
    concordance_order += [c for c in ct.columns if c not in concordance_order]
    ct = ct.reindex(columns=concordance_order, fill_value=0)

    row_totals = ct.sum(axis=1)
    if min_samples > 0:
        keep = row_totals[row_totals >= min_samples].index
        ct = ct.loc[keep]
        row_totals = row_totals.loc[keep]

    ct_frac = ct.div(row_totals, axis=0)
    sort_idx = row_totals.sort_values(ascending=True).index
    ct_frac = ct_frac.loc[sort_idx]
    row_totals = row_totals.loc[sort_idx]

    ct_frac.to_csv(out_csv_prefix.with_name(out_csv_prefix.name + "_proportions.csv"))
    ct.to_csv(out_csv_prefix.with_name(out_csv_prefix.name + "_counts.csv"))

    n_mechs = len(ct_frac)
    figsize = (12, max(4, 0.55 * n_mechs + 1.5))
    fig, ax = plt.subplots(figsize=figsize)

    lefts = np.zeros(n_mechs)
    bars_by_cat = {}
    for cat in concordance_order:
        if cat not in ct_frac.columns:
            continue
        widths = ct_frac[cat].values
        b = ax.barh(
            range(n_mechs),
            widths,
            left=lefts,
            height=0.65,
            color=MECH_COLOR_MAP.get(cat, "#999999"),
            edgecolor="white",
            linewidth=0.5,
            label=cat,
        )
        bars_by_cat[cat] = b
        lefts += widths

    ax.set_yticks(range(n_mechs))
    ax.set_yticklabels(ct_frac.index, fontsize=9)
    for i, (_, n) in enumerate(row_totals.items()):
        ax.text(1.02, i, f"n={n}", va="center", ha="left", fontsize=8, color="#555555")

    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_xlabel("Fraction", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)

    legend_cats = [c for c in concordance_order if c in bars_by_cat and c != "NOT_APPLICABLE"]
    ax.legend(
        [bars_by_cat[c] for c in legend_cats],
        legend_cats,
        loc="lower right",
        fontsize=8,
        framealpha=0.9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(out_png, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_global_z_density_panel(
    snp_df: pd.DataFrame,
    indel_df: pd.DataFrame,
    out_path: Path,
) -> None:
    """Fig 2 — 1×3 KDE panel (SNP, Indel, Combined) matching fig2_density_label.ipynb."""
    neg_label = "Benign"
    pos_label = "Pathogenic"
    color_map = {neg_label: "#4C72B0", pos_label: "#DD8452"}

    concat = pd.concat([snp_df, indel_df], ignore_index=True)
    datasets = {"SNP": snp_df, "Indel": indel_df, "Combined": concat}
    panel_labels = ["a", "b", "c"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    for ax, (panel, (dataset_name, df)) in zip(axes, zip(panel_labels, datasets.items())):
        plot_df = df[["Global_z_sum_log", "label"]].dropna().copy()
        plot_df["label"] = normalize_label_series(plot_df["label"])
        plot_df = plot_df[plot_df["label"].isin([neg_label, pos_label])]

        if plot_df.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(dataset_name)
            continue

        y_true = (plot_df["label"] == pos_label).astype(int)
        y_score = plot_df["Global_z_sum_log"].values
        auc = roc_auc_score(y_true, y_score)
        if auc < 0.5:
            auc = roc_auc_score(y_true, -y_score)
            direction_note = f"Higher values → {neg_label}"
        else:
            direction_note = f"Higher values → {pos_label}"

        x = np.linspace(plot_df["Global_z_sum_log"].min(), plot_df["Global_z_sum_log"].max(), 500)
        for lbl in [neg_label, pos_label]:
            vals = plot_df.loc[plot_df["label"] == lbl, "Global_z_sum_log"].values
            if len(vals) > 1:
                kde = gaussian_kde(vals)
                y = kde(x)
                ax.plot(x, y, label=f"{lbl} (n={len(vals)})", color=color_map[lbl], linewidth=2)
                ax.fill_between(x, y, alpha=0.25, color=color_map[lbl])

        ax.set_title(dataset_name, fontsize=13, pad=12)
        ax.set_xlabel("Global_z_sum_log", fontsize=11)
        ax.text(
            0.03,
            0.95,
            f"ROC AUC = {auc:.3f}\n{direction_note}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )
        ax.text(-0.15, 1.08, panel, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top", ha="left")

    axes[0].set_ylabel("Density", fontsize=11)
    handles, labels = axes[0].get_legend_handles_labels()
    if not handles:
        for ax in axes:
            h, l = ax.get_legend_handles_labels()
            if h:
                handles, labels = h, l
                break
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.05))

    fig.suptitle(
        "Distribution of variants (SNPs and Indels) along MAGI's impact score (global)",
        fontsize=14,
        y=1.12,
    )
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
