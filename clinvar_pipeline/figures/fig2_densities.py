"""Fig 2 — Global_z_sum_log KDE density by label (from fig2_density_label.ipynb)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from config import PATHS
from figures.common import ensure_dirs, normalize_label_series
from figures.plotting import plot_global_z_density_panel


def run() -> None:
    ensure_dirs()
    snp = pd.read_parquet(PATHS["snp_deltas_impact"])
    indel = pd.read_parquet(PATHS["indel_deltas_impact"])
    for name, df in [("SNP", snp), ("Indel", indel)]:
        if "label" not in df.columns:
            raise KeyError(
                f"{name} impact parquet missing 'label' — re-run: "
                "python3 scripts/06_compute_impact_scores.py"
            )
        if "Global_z_sum_log" not in df.columns:
            raise KeyError(
                f"{name} impact parquet missing Global_z_sum_log — re-run: "
                "python3 scripts/06_compute_impact_scores.py"
            )
        df["label"] = normalize_label_series(df["label"])

    out_path = Path(PATHS["fig2b"])
    plot_global_z_density_panel(snp, indel, out_path)
    print(f"  Fig2 → {out_path}")
