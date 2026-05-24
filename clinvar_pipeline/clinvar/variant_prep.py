"""Load NT delta parquets, map ClinVar IDs, z-score MLM columns, strict filter."""

from __future__ import annotations

import gc
import json
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from clinvar.clinvar_enrichment import fetch_rationales_and_filter, load_submission_summary
from clinvar.constants import (
    CLINVAR_SUMMARY_COLS,
    REVIEW_STATUS_TO_GOLD_STARS,
    Z_COLS_INDEL,
    Z_COLS_SNP,
)


def zscore_cols(df: pd.DataFrame, z_cols: list[str]) -> pd.DataFrame:
    for col in z_cols:
        if col in df.columns:
            df[col] = (df[col] - df[col].mean()) / df[col].std()
    return df


# BED/BW delta + reference columns are loaded after strict filtering (memory).
TRACK_COLUMN_PREFIXES = (
    "D_BED_",
    "D_BW_",
    "REF_BED_",
    "ALT_BED_",
    "REF_BW_",
    "ALT_BW_",
)
_ROW_COL = "_pq_row"
_DEFAULT_TRACK_COL_BATCH = 120


def is_track_column(name: str) -> bool:
    return name.startswith(TRACK_COLUMN_PREFIXES)


def _track_columns(names: list[str]) -> list[str]:
    return [c for c in names if is_track_column(c)]


def _slim_columns(names: list[str]) -> list[str]:
    track = set(_track_columns(names))
    return [c for c in names if c not in track]


def attach_track_columns(
    df: pd.DataFrame,
    path: str,
    *,
    row_col: str = _ROW_COL,
    column_batch_size: int = _DEFAULT_TRACK_COL_BATCH,
) -> pd.DataFrame:
    """Load deferred BED/BW columns for rows already present in *df*."""
    if row_col not in df.columns:
        raise KeyError(f"Missing row index column {row_col!r}")

    base = df.drop(columns=[row_col])
    if df.empty:
        return base

    pf = pq.ParquetFile(path)
    track_cols = _track_columns(pf.schema.names)
    if not track_cols:
        return base

    indices = pa.array(df[row_col].astype("int64").tolist())
    n_track = len(track_cols)
    bed_n = sum(c.startswith("D_BED_") for c in track_cols)
    print(
        f"  Attaching {n_track:,} track columns ({bed_n} BED deltas) "
        f"for {len(df):,} variants ..."
    )

    parts: list[pd.DataFrame] = []
    for start in range(0, n_track, column_batch_size):
        batch_cols = track_cols[start : start + column_batch_size]
        table = pq.read_table(path, columns=batch_cols)
        taken = table.take(indices)
        parts.append(taken.to_pandas())
        del table, taken
        gc.collect()

    track_df = pd.concat(parts, axis=1)
    track_df.index = base.index
    return pd.concat([base, track_df], axis=1)


def load_snp_deltas(path: str, *, defer_tracks: bool = True) -> pd.DataFrame:
    pf = pq.ParquetFile(path)
    all_cols = pf.schema.names
    if defer_tracks:
        slim_cols = _slim_columns(all_cols)
        track_cols = _track_columns(all_cols)
        bed_n = sum(c.startswith("D_BED_") for c in track_cols)
        print(
            f"  Loading {len(slim_cols)} slim columns "
            f"({len(track_cols)} track cols deferred, {bed_n} BED deltas) ..."
        )
        df = pf.read(columns=slim_cols).to_pandas()
        df[_ROW_COL] = np.arange(len(df), dtype=np.int64)
    else:
        delta_cols = [c for c in all_cols if c.startswith("D_BED_")]
        print(f"  Loading {len(all_cols)} columns ({len(delta_cols)} BED features) ...")
        df = pf.read(columns=all_cols).to_pandas()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    return df


def load_indel_deltas(path: str, *, defer_tracks: bool = True) -> pd.DataFrame:
    pf = pq.ParquetFile(path)
    all_cols = pf.schema.names
    if defer_tracks:
        slim_cols = _slim_columns(all_cols)
        track_cols = _track_columns(all_cols)
        print(
            f"  Loading {len(slim_cols)} slim columns "
            f"({len(track_cols)} track cols deferred) ..."
        )
        df = pf.read(columns=slim_cols).to_pandas()
        df[_ROW_COL] = np.arange(len(df), dtype=np.int64)
    else:
        print(f"  Loading {len(all_cols)} columns ...")
        df = pf.read(columns=all_cols).to_pandas()
    df.rename(columns={"variant_id": "#VariationID"}, inplace=True)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    return df


def _normalize_variant_summary(var_df: pd.DataFrame) -> pd.DataFrame:
    var_df.rename(
        columns={
            "Chromosome": "chrom",
            "PositionVCF": "pos",
            "ReferenceAlleleVCF": "ref",
            "AlternateAlleleVCF": "alt",
            "VariationID": "#VariationID",
            "ReviewStatus": "review_status",
        },
        inplace=True,
    )
    var_df["chrom"] = var_df["chrom"].astype(str)
    if len(var_df) and not str(var_df["chrom"].iloc[0]).startswith("chr"):
        var_df["chrom"] = "chr" + var_df["chrom"]
    var_df["pos"] = pd.to_numeric(var_df["pos"], errors="coerce").astype("Int64")
    var_df["gold_stars"] = var_df["review_status"].map(REVIEW_STATUS_TO_GOLD_STARS)
    return var_df


def _load_variant_summary_filtered(
    variant_summary_path: str,
    *,
    assembly: str | None = None,
    variant_types: set[str] | None = None,
) -> pd.DataFrame:
    """Stream variant_summary and keep only rows matching assembly/type filters."""
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        variant_summary_path,
        sep="\t",
        compression="gzip",
        usecols=CLINVAR_SUMMARY_COLS,
        chunksize=250_000,
        low_memory=False,
    ):
        if assembly is not None:
            chunk = chunk[chunk["Assembly"] == assembly]
        if variant_types is not None:
            chunk = chunk[chunk["Type"].isin(variant_types)]
        if not chunk.empty:
            chunks.append(chunk)

    if not chunks:
        return _normalize_variant_summary(
            pd.DataFrame(columns=[c for c in CLINVAR_SUMMARY_COLS if c != "Assembly"])
        )

    var_df = pd.concat(chunks, ignore_index=True)
    del chunks
    gc.collect()
    return _normalize_variant_summary(var_df)


def map_snp_variation_ids(delta_df: pd.DataFrame, variant_summary_path: str) -> pd.DataFrame:
    print("  Loading ClinVar variant_summary (GRCh38 SNVs only) ...")
    snp_df = _load_variant_summary_filtered(
        variant_summary_path,
        assembly="GRCh38",
        variant_types={"single nucleotide variant"},
    )

    print("  Merging VariationIDs ...")
    merged = delta_df.merge(
        snp_df[["chrom", "pos", "ref", "alt", "#VariationID", "GeneSymbol", "gold_stars"]],
        on=["chrom", "pos", "ref", "alt"],
        how="left",
    )
    del snp_df
    gc.collect()
    n_matched = merged["#VariationID"].notna().sum()
    print(f"  Matched {n_matched:,} / {len(merged):,} variants to ClinVar IDs")
    return merged


def map_indels_variation_ids(delta_df: pd.DataFrame, variant_summary_path: str) -> pd.DataFrame:
    print("  Loading ClinVar variant_summary (indels only) ...")
    var_df = _load_variant_summary_filtered(
        variant_summary_path,
        variant_types={"Deletion", "Insertion", "Indel"},
    )
    indel_df = var_df

    gene_symbol_per_id = (
        indel_df.dropna(subset=["GeneSymbol"])
        .drop_duplicates(subset=["#VariationID"])[["#VariationID", "GeneSymbol", "review_status"]]
    )

    print("  Merging GeneSymbols ...")
    merged = delta_df.merge(gene_symbol_per_id, on="#VariationID", how="left")
    merged["gold_stars"] = merged["review_status"].map(REVIEW_STATUS_TO_GOLD_STARS)

    n_matched = merged["GeneSymbol"].notna().sum()
    print(f"  Matched {n_matched:,} / {len(merged):,} variants to ClinVar IDs")
    return merged


def measure_cohort_funnel(
    paths: dict,
    *,
    run_snp: bool = True,
    run_indel: bool = True,
) -> dict[str, dict[str, int]]:
    """Count quality vs rationaled cohorts without writing parquets (for figures)."""
    funnel: dict[str, dict[str, int]] = {}
    sub_df = None

    if run_snp:
        sub_df = load_submission_summary(paths["submission_summary"])
        delta_df = load_snp_deltas(paths["snp_deltas"], defer_tracks=True)
        delta_df = zscore_cols(delta_df, Z_COLS_SNP)
        mapped_df = map_snp_variation_ids(delta_df, paths["variant_summary"])
        del delta_df
        gc.collect()
        _, stats = fetch_rationales_and_filter(
            mapped_df, paths["submission_summary"], sub_df=sub_df
        )
        funnel["snp"] = stats
        del mapped_df, sub_df
        sub_df = None
        gc.collect()

    if run_indel:
        sub_df = load_submission_summary(paths["submission_summary"])
        indel_df = load_indel_deltas(paths["indel_deltas"], defer_tracks=True)
        indel_df = zscore_cols(indel_df, Z_COLS_INDEL)
        mapped_df = map_indels_variation_ids(indel_df, paths["variant_summary"])
        del indel_df
        gc.collect()
        _, stats = fetch_rationales_and_filter(
            mapped_df, paths["submission_summary"], sub_df=sub_df
        )
        funnel["indel"] = stats
        del mapped_df, sub_df
        gc.collect()

    return funnel


def write_cohort_funnel(paths: dict, funnel: dict[str, dict[str, int]]) -> None:
    from pathlib import Path

    out = Path(paths["cohort_funnel"])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(funnel, indent=2) + "\n")
    print(f"  Cohort funnel → {out}")


def process_snps(paths: dict, sub_df: Optional[pd.DataFrame] = None) -> tuple[pd.DataFrame, dict[str, int]]:
    print("\n══ PART 1: SNP Variants ══")
    delta_df = load_snp_deltas(paths["snp_deltas"])
    delta_df = zscore_cols(delta_df, Z_COLS_SNP)
    mapped_df = map_snp_variation_ids(delta_df, paths["variant_summary"])
    del delta_df
    gc.collect()
    strict_df, funnel = fetch_rationales_and_filter(
        mapped_df, paths["submission_summary"], sub_df=sub_df
    )
    del mapped_df
    gc.collect()
    if _ROW_COL in strict_df.columns:
        strict_df = attach_track_columns(strict_df, paths["snp_deltas"])
    strict_df.to_parquet(paths["snp_strict"], index=False)
    print(f"  Saved → {paths['snp_strict']}")
    return strict_df, funnel


def process_indels(paths: dict, sub_df: Optional[pd.DataFrame] = None) -> tuple[pd.DataFrame, dict[str, int]]:
    print("\n══ PART 2: Indel Variants ══")
    indel_df = load_indel_deltas(paths["indel_deltas"])
    indel_df = zscore_cols(indel_df, Z_COLS_INDEL)
    mapped_df = map_indels_variation_ids(indel_df, paths["variant_summary"])
    del indel_df
    gc.collect()
    strict_df, funnel = fetch_rationales_and_filter(
        mapped_df, paths["submission_summary"], sub_df=sub_df
    )
    del mapped_df
    gc.collect()
    if _ROW_COL in strict_df.columns:
        strict_df = attach_track_columns(strict_df, paths["indel_deltas"])
    strict_df.to_parquet(paths["indel_strict"], index=False)
    print(f"  Saved → {paths['indel_strict']}")
    return strict_df, funnel
