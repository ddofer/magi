"""Fig 3a–c — concordance overview and Global_z quartile bars (from fig3_a_b_c.ipynb)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from config import PATHS
from figures.common import DPI, ensure_dirs
from figures.merge_eval import load_merged_eval_frame
from figures.plotting import plot_concordance_bars, plot_overall_concordance_by_variant_type


def run() -> None:
    ensure_dirs()
    merged = load_merged_eval_frame()

    plot_overall_concordance_by_variant_type(merged, Path(PATHS["fig3a"]))

    quartile_outputs = {
        "snp": PATHS["fig3b"],
        "indel": PATHS["fig3c"],
    }
    for vt_key, out_png in quartile_outputs.items():
        sub = merged[merged["variant_type"].astype(str).str.lower() == vt_key].copy()
        if sub.empty or "Global_z_sum_log" not in sub.columns:
            print(f"  Fig3 {vt_key} quartiles skipped: no Global_z_sum_log")
            continue
        vt_label = vt_key.upper()
        fig, ax, df_out = plot_concordance_bars(
            sub,
            group_col="Global_z_sum_log",
            normalize=True,
            title=f"{vt_label}: Concordance by Global_z_sum_log quartiles",
        )
        fig.savefig(out_png, dpi=DPI, bbox_inches="tight", facecolor="white")
        csv_path = Path(out_png).with_suffix(".csv")
        df_out.to_csv(csv_path)
        plt.close(fig)

    print(f"  Fig3abc → {PATHS['fig3a']}, {PATHS['fig3b']}, {PATHS['fig3c']}")
