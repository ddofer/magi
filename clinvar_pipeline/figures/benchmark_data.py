"""Load labelled benchmark frames for fig2a and related plots."""

from __future__ import annotations

import pandas as pd

from config import PATHS
from figures.common import normalize_label_series


def load_impact_benchmark() -> pd.DataFrame:
    """Impact-scored benchmark (Global_z columns + AF) from ``output/impact/``."""
    snp = pd.read_parquet(PATHS["snp_deltas_impact"])
    indel = pd.read_parquet(PATHS["indel_deltas_impact"])
    df = pd.concat([snp, indel], ignore_index=True)
    df["label"] = normalize_label_series(df["label"])
    return df[df["label"].isin(["Pathogenic", "Benign"])].copy()
