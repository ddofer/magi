"""Global impact scores with benign-referenced z-scoring (per variant type).

Ported from ``impact_score.ipynb``. Z-score background: **Benign-labelled rows**
within each ``variant_type`` group (SNP vs indel types pooled for indels).
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from config import IMPACT_MERGE_KEYS, IMPACT_METADATA_COLS, IMPACT_SCORE_COLS

TRACK_PREFIXES = ("D_BED_", "D_BW_")


def _delta_columns(names: Iterable[str]) -> list[str]:
    return sorted(c for c in names if c.startswith(TRACK_PREFIXES))


def _read_delta_subset(path: str) -> pd.DataFrame:
    pf = pq.ParquetFile(path)
    optional = ["label", "variant_type", "variant_id", "#VariationID"]
    keep = list(dict.fromkeys([*IMPACT_MERGE_KEYS, *optional] + _delta_columns(pf.schema.names)))
    keep = [c for c in keep if c in pf.schema.names]
    return pf.read(columns=keep).to_pandas()


def _harmonize_labels(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out["label"].dtype == bool:
        out["label"] = out["label"].map({True: "Pathogenic", False: "Benign"})
    elif out["label"].dtype in (np.int64, np.int32, float):
        out["label"] = out["label"].map({1: "Pathogenic", 0: "Benign", True: "Pathogenic", False: "Benign"})
    return out


def z_score_matrix(
    df: pd.DataFrame,
    cols: list[str],
    benign_mask: np.ndarray,
    *,
    per_type: bool = True,
) -> np.ndarray:
    """Z-score columns against benign distribution (optionally per variant_type)."""
    matrix = df[cols].values.astype(np.float64)
    result = np.empty_like(matrix)

    groups = (
        df["variant_type"].unique()
        if per_type and "variant_type" in df.columns
        else [None]
    )
    for vt in groups:
        idx = (
            (df["variant_type"] == vt).values
            if vt is not None
            else np.ones(len(df), dtype=bool)
        )
        ben = idx & benign_mask
        if ben.sum() < 10:
            ben = benign_mask
        mu = np.nanmean(matrix[ben], axis=0)
        sigma = np.nanstd(matrix[ben], axis=0)
        sigma[sigma < 1e-8] = 1.0
        result[idx] = (matrix[idx] - mu) / sigma
    return result


def agg_scores(matrix: np.ndarray, prefix: str) -> dict[str, np.ndarray]:
    am = np.abs(matrix)
    return {
        f"{prefix}_mean_abs": np.nanmean(am, axis=1),
        f"{prefix}_max_abs": np.nanmax(am, axis=1),
        f"{prefix}_sum_log": np.nansum(np.log1p(am), axis=1),
    }


def top_k_mean(matrix: np.ndarray, k: int) -> np.ndarray:
    k = min(k, matrix.shape[1])
    if k == 0:
        return np.zeros(matrix.shape[0])
    return np.sort(np.abs(matrix), axis=1)[:, -k:].mean(axis=1)


def compute_impact_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Add BED/BW/Global_z impact columns to a labelled delta dataframe."""
    df = _harmonize_labels(df)
    d_bed = sorted(c for c in df.columns if c.startswith("D_BED_"))
    d_bw = sorted(c for c in df.columns if c.startswith("D_BW_"))
    benign_mask = (df["label"] == "Benign").values

    if d_bed:
        bed_raw = df[d_bed].values.astype(np.float64)
        for name, arr in agg_scores(np.abs(bed_raw), "BED").items():
            df[name] = arr
        df["BED_top3_mean"] = top_k_mean(bed_raw, 3)

    if d_bw:
        bw_z = z_score_matrix(df, d_bw, benign_mask, per_type=True)
        np.nan_to_num(bw_z, copy=False, nan=0.0)
        bw_abs = np.abs(bw_z)
        for name, arr in agg_scores(bw_abs, "BW_z").items():
            df[name] = arr
        df["BW_z_top10_mean"] = top_k_mean(bw_z, 10)

    if d_bed and d_bw:
        bed_z = z_score_matrix(df, d_bed, benign_mask, per_type=True)
        np.nan_to_num(bed_z, copy=False, nan=0.0)
        global_abs = np.hstack([np.abs(bed_z), np.abs(bw_z)])
        for name, arr in agg_scores(global_abs, "Global_z").items():
            df[name] = arr
        df["Composite_mean"] = 0.5 * df["BED_mean_abs"] + 0.5 * df["BW_z_mean_abs"]
        df["Composite_top"] = 0.5 * df["BED_top3_mean"] + 0.5 * df["BW_z_top10_mean"]
    return df


def _merge_af(df: pd.DataFrame, af_paths: list[str]) -> pd.DataFrame:
    frames = [pd.read_csv(p) for p in af_paths if p]
    if not frames:
        df["gnomADe_AF"] = np.nan
        return df
    af = pd.concat(frames, ignore_index=True)
    for c in [c for c in af.columns if c.startswith("gnomAD") or c == "UKBB_AF"]:
        af[c] = pd.to_numeric(af[c].replace("-", np.nan), errors="coerce")
    af = af.drop_duplicates(subset=IMPACT_MERGE_KEYS, keep="first")
    n = len(df)
    merged = df.merge(af[[*IMPACT_MERGE_KEYS, "gnomADe_AF"]], on=IMPACT_MERGE_KEYS, how="left")
    if len(merged) != n:
        raise RuntimeError(f"AF merge changed row count: {n} → {len(merged)}")
    return merged


def prepare_combined_deltas(
    snp_path: str,
    indel_path: str,
    af_paths: list[str],
) -> pd.DataFrame:
    snp = _read_delta_subset(snp_path)
    indel = _read_delta_subset(indel_path)
    snp["variant_type"] = "SNP"
    if "variant_type" not in indel.columns:
        indel["variant_type"] = "Indel"
    df = pd.concat([snp, indel], ignore_index=True)
    return _merge_af(df, af_paths)


def slim_impact_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in IMPACT_MERGE_KEYS if c in df.columns]
    meta = [c for c in IMPACT_METADATA_COLS if c in df.columns]
    extra = [c for c in IMPACT_SCORE_COLS if c in df.columns]
    return df[cols + meta + extra].copy()


def attach_impact_to_signaled(
    signaled_path: str,
    impact_path: str,
) -> pd.DataFrame:
    """Join impact score columns onto a signaled cohort parquet (in memory)."""
    signaled = pd.read_parquet(signaled_path)
    impact = pd.read_parquet(impact_path)
    merge_cols = [c for c in IMPACT_MERGE_KEYS if c in signaled.columns and c in impact.columns]
    score_cols = [c for c in IMPACT_SCORE_COLS if c in impact.columns]
    return signaled.merge(impact[merge_cols + score_cols], on=merge_cols, how="left")


def merge_impact_onto_signaled(
    signaled_path: str,
    impact_slim: pd.DataFrame,
) -> pd.DataFrame:
    """Deprecated alias — prefer ``attach_impact_to_signaled``."""
    signaled = pd.read_parquet(signaled_path)
    merge_cols = [c for c in IMPACT_MERGE_KEYS if c in signaled.columns]
    score_cols = [c for c in IMPACT_SCORE_COLS if c in impact_slim.columns]
    return signaled.merge(
        impact_slim[merge_cols + score_cols],
        on=merge_cols,
        how="left",
    )
