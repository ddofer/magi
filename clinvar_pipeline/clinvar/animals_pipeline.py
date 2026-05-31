"""OMIA animal variant signal extraction (deltas → signaled parquets)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from clinvar.signal_extraction import extract_signals_comprehensive_fast


def _zscore_column(df: pd.DataFrame, col: str) -> None:
    if col not in df.columns:
        return
    std = df[col].std()
    if std and std > 0:
        df[col] = (df[col] - df[col].mean()) / std


def extract_animals_signals(
    paths: dict,
    variant_type: str,
    *,
    k_signals: int = 10,
) -> pd.DataFrame:
    """Read animal delta parquet, extract top-k signals, write signaled parquet."""
    if variant_type == "snp":
        input_file = paths["animals_snp_deltas"]
        output_file = paths["animals_snp_signaled"]
    elif variant_type == "indel":
        input_file = paths["animals_indel_deltas"]
        output_file = paths["animals_indel_signaled"]
    else:
        raise ValueError(f"Unknown variant type: {variant_type}")

    print(f"Loading {variant_type} animal deltas from {input_file} ...")
    variant_df = pq.read_table(input_file).to_pandas()
    print(f"Loaded {len(variant_df):,} variants")

    metadata_df = None
    metadata_path = paths.get("tracks_metadata")
    if metadata_path and Path(metadata_path).exists():
        metadata_df = pd.read_csv(metadata_path)
        print(f"Loaded metadata for {len(metadata_df):,} tracks")

    delta_cols = [
        c for c in variant_df.columns if c.startswith("D_BED_") or c.startswith("D_BW_")
    ]
    print(f"Found {sum(c.startswith('D_BED_') for c in delta_cols)} BED and "
          f"{sum(c.startswith('D_BW_') for c in delta_cols)} BW delta columns")

    for col in delta_cols:
        variant_df[col] = pd.to_numeric(variant_df[col], errors="coerce")
    for col in ("MLM_logprob_delta", "MLM_logprob_ref"):
        _zscore_column(variant_df, col)

    print("Extracting signals (vectorized) ...")
    signal_df = extract_signals_comprehensive_fast(
        variant_df,
        delta_cols,
        metadata_df,
        k=k_signals,
        variant_type=variant_type,
    )
    for col in signal_df.columns:
        variant_df[col] = signal_df[col].values

    variant_df = variant_df.reset_index()
    rename_map = {"index": "#VariationID"}
    if "rationale" in variant_df.columns:
        rename_map["rationale"] = "FullRationale"
    variant_df = variant_df.rename(columns=rename_map)

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    variant_df.to_parquet(out_path, index=False)
    print(f"Saved → {out_path}  ({len(variant_df):,} variants)")
    return variant_df
