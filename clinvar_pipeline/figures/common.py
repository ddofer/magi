"""Shared paths and label helpers for figure scripts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from config import PATHS

FIGURES_DIR = Path(PATHS["figures_root"])
DPI = 300

# Shared concordance palette (green / orange / red) for fig3b–3e
CONCORDANCE_COLORS = {
    "CONCORDANT": "forestgreen",
    "PARTIAL": "goldenrod",
    "DISCORDANT": "firebrick",
    "NOT_APPLICABLE": "#bbbbbb",
}

DEFAULT_CONCORDANCE_ORDER = ("CONCORDANT", "PARTIAL", "DISCORDANT")


def concordance_color_list(categories: list[str] | tuple[str, ...]) -> list[str]:
    return [CONCORDANCE_COLORS.get(c, "#999999") for c in categories]


def ensure_dirs() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def normalize_label_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.map({True: "Pathogenic", False: "Benign"})
    if pd.api.types.is_numeric_dtype(s):
        return s.map({1: "Pathogenic", 0: "Benign", True: "Pathogenic", False: "Benign"})
    out = s.astype(str)
    out = out.replace({"True": "Pathogenic", "False": "Benign", "1": "Pathogenic", "0": "Benign"})
    return out


def load_annotated_pair() -> tuple[pd.DataFrame, pd.DataFrame]:
    snp = pd.read_parquet(PATHS["snp_annotated"])
    indel = pd.read_parquet(PATHS["indel_annotated"])
    for df in (snp, indel):
        if "label" in df.columns:
            df["label"] = normalize_label_series(df["label"])
    return snp, indel


def tag_variant_type(snp: pd.DataFrame, indel: pd.DataFrame) -> pd.DataFrame:
    snp = snp.copy()
    indel = indel.copy()
    snp["variant_class"] = "snp"
    indel["variant_class"] = "indel"
    return pd.concat([snp, indel], ignore_index=True)
