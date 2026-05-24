"""Load labelled benchmark frames for fig2a and related plots."""

from __future__ import annotations

import pandas as pd

from config import PATHS
from clinvar.impact_scoring import (
    _harmonize_labels,
    _merge_af,
    _read_delta_subset,
    prepare_combined_deltas,
)
from figures.common import normalize_label_series


def load_benchmark_deltas(*, bed_only: bool = False, snp_only: bool = False) -> pd.DataFrame:
    """Raw labelled delta benchmark (unfiltered cohort)."""
    if snp_only:
        df = _read_delta_subset(PATHS["snp_deltas"])
        df["variant_type"] = "SNP"
        df = _merge_af(df, [PATHS["af_snps"]])
    else:
        df = prepare_combined_deltas(
            PATHS["snp_deltas"],
            PATHS["indel_deltas"],
            [PATHS["af_snps"], PATHS["af_indels"]],
        )
    df = _harmonize_labels(df)
    df["label"] = normalize_label_series(df["label"])
    df = df[df["label"].isin(["Pathogenic", "Benign"])].copy()
    if "variant_type" not in df.columns:
        df["variant_type"] = "Indel"
    df.loc[df["variant_type"].isna(), "variant_type"] = "SNP"
    if bed_only:
        keep = [c for c in df.columns if c.startswith("D_BED_") or c in ("label", "variant_type")]
        keep += [c for c in ("chrom", "pos", "ref", "alt", "gnomADe_AF") if c in df.columns]
        df = df[list(dict.fromkeys(keep))]
    return df


def load_impact_benchmark() -> pd.DataFrame:
    """Impact-scored benchmark (Global_z columns + AF) from ``output/impact/``."""
    snp = pd.read_parquet(PATHS["snp_deltas_impact"])
    indel = pd.read_parquet(PATHS["indel_deltas_impact"])
    df = pd.concat([snp, indel], ignore_index=True)
    df["label"] = normalize_label_series(df["label"])
    return df[df["label"].isin(["Pathogenic", "Benign"])].copy()
