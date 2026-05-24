"""Fig 1 — label and rationale coverage panel (from fig1.ipynb) → fig1c.png."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from config import PATHS
from figures.common import DPI, ensure_dirs, load_annotated_pair, tag_variant_type


def _load_cohort_funnel() -> dict[str, dict[str, int]]:
    path = Path(PATHS["cohort_funnel"])
    if path.exists():
        return json.loads(path.read_text())
    from clinvar.variant_prep import measure_cohort_funnel, write_cohort_funnel

    print("  Computing cohort funnel (quality vs rationaled) …")
    funnel = measure_cohort_funnel(PATHS)
    write_cohort_funnel(PATHS, funnel)
    return funnel


def _draw_label_pie(ax, df: pd.DataFrame, title: str, colors: dict[str, str]) -> None:
    counts = df["label"].value_counts()
    raw = [(lbl, int(counts.get(lbl, 0)), colors[lbl]) for lbl in colors]
    pie_labels = [l for l, c, _ in raw if c > 0]
    pie_data = [c for _, c, _ in raw if c > 0]
    pie_colors = [col for _, c, col in raw if c > 0]
    if not pie_data:
        ax.set_title(f"{title} (no data)")
        return
    explode = [0.10 if l == "VUS" else 0 for l in pie_labels]
    ax.pie(
        pie_data,
        labels=pie_labels,
        autopct="%1.1f%%",
        startangle=140,
        colors=pie_colors,
        explode=explode,
    )
    ax.set_title(title)


def _draw_rationale_pie(ax, total: int, with_rationale: int, variant_name: str) -> None:
    """Quality cohort vs variants with detailed ClinVar rationale (fig1.ipynb)."""
    missing = max(0, total - with_rationale)
    pie_data = [with_rationale, missing]
    pie_labels = ["Detailed\nRationale", "Missing Rationale\n(Knowledge Gap)"]
    colors = ["#7BAF2B", "#E25815"]
    explode = [0, 0.08]

    wedges, texts, autotexts = ax.pie(
        pie_data,
        labels=pie_labels,
        colors=colors,
        explode=explode,
        startangle=140,
        autopct="%1.0f%%",
        pctdistance=0.48,
        labeldistance=0.48,
        wedgeprops={"edgecolor": "black", "linewidth": 1.0},
        textprops={"fontsize": 11, "fontweight": "bold", "color": "black"},
    )
    if len(texts) > 1:
        texts[1].set_text("")
    if len(autotexts) > 1 and total > 0:
        pct_missing = 100 * missing / total
        autotexts[0].set_text("")
        autotexts[1].set_text(f"~{round(pct_missing):d}%")
    ax.set_title(f"{variant_name}: ClinVar rationale coverage")


def run() -> None:
    ensure_dirs()
    snp_ann, indel_ann = load_annotated_pair()
    all_labeled = tag_variant_type(snp_ann, indel_ann)
    funnel = _load_cohort_funnel()

    vus_path = Path(PATHS["uncertain_ids"])
    if vus_path.exists():
        vus_ids = pd.read_csv(vus_path)["#VariationID"].tolist()
        mask = all_labeled["#VariationID"].isin(vus_ids)
        all_labeled.loc[mask, "label"] = "VUS"

    colors = {"Pathogenic": "#C1602D", "Benign": "#89AD47", "VUS": "#E0B646"}
    subsets = {
        "snp": all_labeled[all_labeled["variant_class"] == "snp"],
        "indel": all_labeled[all_labeled["variant_class"] == "indel"],
    }

    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    _draw_label_pie(axes[0, 0], subsets["snp"], "SNP labels", colors)
    _draw_label_pie(axes[0, 1], subsets["indel"], "Indel labels", colors)
    for ax, (name, key) in zip(axes[1], [("SNP", "snp"), ("Indel", "indel")]):
        stats = funnel.get(key, {})
        total = int(stats.get("quality", 0))
        with_rationale = int(stats.get("rationaled", 0))
        if total > 0:
            _draw_rationale_pie(ax, total, with_rationale, name)
        else:
            ax.set_title(f"{name} rationale (no funnel data)")

    out = Path(PATHS["fig1c"])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"  Fig1 → {out}")
