"""Fig 3d — indel concordance by frameshift modulo (from fig3d_indel_frame.ipynb)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from config import PATHS
from clinvar.impact_scoring import attach_impact_to_signaled
from figures.common import DPI, concordance_color_list, ensure_dirs
from figures.merge_eval import load_merged_eval_frame


def _resolve_indel_size(df: pd.DataFrame) -> pd.Series:
    for col in ("size", "indel_size"):
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce")
    signaled = attach_impact_to_signaled(
        PATHS["indel_annotated_signaled"],
        PATHS["indel_deltas_impact"],
    )
    if "indel_size" in signaled.columns and "#VariationID" in df.columns:
        lookup = signaled.drop_duplicates("#VariationID").set_index("#VariationID")["indel_size"]
        return df["#VariationID"].map(lookup)
    raise KeyError("indel_size/size not found in LLM results or signaled parquet")


def run() -> None:
    ensure_dirs()
    merged = load_merged_eval_frame()
    df = merged[merged["variant_type"].astype(str).str.lower() == "indel"].copy()
    if df.empty:
        print("  Fig3d skipped: no indel rows in merged eval frame")
        return

    df["size"] = _resolve_indel_size(df)
    df = df[df["size"].notna()].copy()
    if df.empty:
        print("  Fig3d skipped: could not resolve indel sizes")
        return

    df["mod3_label"] = (df["size"].abs() % 3).map({0: "0 (in-frame)", 1: "1", 2: "2"})
    df["concordance"] = df["concordance"].astype(str).str.strip().str.upper()

    row_order = ["0 (in-frame)", "1", "2"]
    col_order = ["CONCORDANT", "PARTIAL", "DISCORDANT"]
    mod3_props = (
        pd.crosstab(df["mod3_label"], df["concordance"], normalize="index") * 100
    ).reindex(index=row_order, columns=col_order, fill_value=0)

    fig, ax = plt.subplots(figsize=(8, 5))
    mod3_props.plot(
        kind="bar",
        stacked=True,
        color=concordance_color_list(col_order),
        edgecolor="white",
        ax=ax,
        legend=True,
    )
    counts_per_mod3 = df["mod3_label"].value_counts().reindex(row_order)
    ax.set_xticklabels(
        [f"{label}\n(n={int(counts_per_mod3[label])})" for label in row_order],
        rotation=0,
    )
    ax.set_ylabel("Percentage")
    ax.set_xlabel("|size| mod 3")
    ax.set_title("Concordance distribution by indel frame shift")
    ax.legend(title="Concordance", bbox_to_anchor=(1.05, 1), loc="upper left")
    out = Path(PATHS["fig3d"])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Fig3d → {out}")
