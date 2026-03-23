#!/usr/bin/env python3
"""
Analysis and visualization module for the MAGI Gradio app.

Adapted from elegant_solution.py for the app.

Provides:
- Impact score computation (top-K mean absolute delta)
- Feature summarization with direction labels
- Fingerprint-style plots for ranked signals
- BigWig track description mapping

Usage:
    from analysis import compute_impact_scores, make_fingerprint_plot
    df = compute_impact_scores(df)
    fig = make_fingerprint_plot(row, bed_names, bw_names, metadata_df)
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# Use non-interactive backend for server deployment
matplotlib.use("Agg")


BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
MAGI_BASELINE_FILE = DATA_DIR / "magi_baseline_stats.csv"
_MAGI_BASELINE_CACHE: Optional[Dict[str, Dict[str, pd.Series]]] = None

MLM_SIGNAL_SPECS = [
    ("LLR", "LLR", True),
    ("MLM_logprob_delta", "Log-prob Δ", True),
    ("MLM_Delta", "MLM Δ", True),
    ("MLM_KL_mean", "KL mean", False),
    ("MLM_KL_max", "KL max", False),
    ("EMB_cosine_dist", "Embedding cosine dist", False),
    ("EMB_l2_dist", "Embedding L2 dist", False),
    ("EMB_max_pos_dist", "Embedding max-pos dist", False),
    ("EMB_mean_pos_dist", "Embedding mean-pos dist", False),
]


# ============================================================================
# DATA STRUCTURES
# ============================================================================
@dataclass
class MechanisticSummary:
    """Concise summary of top disrupted features for a variant."""

    variant_id: str
    gene_name: str
    region_class: str
    top_features: List[Tuple[str, float, str]]  # (name, delta, GOF/LOF)
    impact_score_bed: float
    llr: float
    description: Optional[str] = None


# ============================================================================
# IMPACT SCORING
# ============================================================================
def identify_delta_columns(df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """
    Identify BED and BigWig delta columns in DataFrame.

    Returns:
        (bed_delta_cols, bw_delta_cols)
    """
    delta_bed_cols = [c for c in df.columns if c.startswith("D_BED_")]
    delta_bw_cols = [c for c in df.columns if c.startswith("D_BW_")]
    return delta_bed_cols, delta_bw_cols


def compute_impact_scores(
    df: pd.DataFrame,
    bed_k: int = 3,
    bw_k: int = 10,
    magi_baseline: Optional[Dict[str, Dict[str, pd.Series]]] = None,
) -> pd.DataFrame:
    """
    Compute impact scores from delta columns.

    Impact_Score_BED: mean of top-K largest absolute BED deltas
    Impact_Score_BW: mean of top-K largest absolute BigWig deltas

    Args:
        df: DataFrame with D_BED_* and D_BW_* columns
        bed_k: Number of top BED deltas to average
        bw_k: Number of top BigWig deltas to average

    Returns:
        DataFrame with added Impact_Score_BED, Impact_Score_BW, and
        Global_z_sum_log columns
    """
    df = df.copy()
    bed_cols, bw_cols = identify_delta_columns(df)

    # BED impact — use argpartition for O(n) top-K selection instead of O(n log n) sort
    if bed_cols:
        bed_matrix = np.abs(df[bed_cols].to_numpy())
        k_bed = min(bed_k, bed_matrix.shape[1])
        if k_bed > 0:
            # argpartition gives us indices that partition the array at position k
            # We want the k largest elements, so we partition at (n - k) to get them on the right
            indices = np.argpartition(bed_matrix, kth=-k_bed, axis=1)[:, -k_bed:]
            top_abs_bed = np.take_along_axis(bed_matrix, indices, axis=1)
            df["Impact_Score_BED"] = top_abs_bed.mean(axis=1)
        else:
            df["Impact_Score_BED"] = np.nan
    else:
        df["Impact_Score_BED"] = np.nan

    # BigWig impact — use argpartition for O(n) top-K selection
    if bw_cols:
        bw_matrix = np.abs(df[bw_cols].to_numpy())
        k_bw = min(bw_k, bw_matrix.shape[1])
        if k_bw > 0:
            indices = np.argpartition(bw_matrix, kth=-k_bw, axis=1)[:, -k_bw:]
            top_abs_bw = np.take_along_axis(bw_matrix, indices, axis=1)
            df["Impact_Score_BW"] = top_abs_bw.mean(axis=1)
        else:
            df["Impact_Score_BW"] = np.nan
    else:
        df["Impact_Score_BW"] = np.nan

    baseline = magi_baseline if magi_baseline is not None else load_magi_baseline()
    if baseline:
        df["Global_z_sum_log"] = df.apply(
            lambda row: compute_global_z_sum_log(row, baseline), axis=1
        )
    else:
        df["Global_z_sum_log"] = np.nan

    return df


def _normalize_variant_type_key(
    value: Optional[str],
    row: Optional[pd.Series] = None,
) -> str:
    """Normalize variant-type labels for MAGI baseline lookup."""
    text = (
        str(value).strip().lower() if value is not None and str(value).strip() else ""
    )
    if text in {"snp", "snv", "single nucleotide variant"}:
        return "snp"
    if text in {"deletion", "del"} or ("deletion" in text and "insert" not in text):
        return "deletion"
    if text in {"insertion", "ins"} or "insertion" in text:
        return "insertion"
    if text in {"indel", "delins"} or "delins" in text:
        return "indel"

    if row is not None:
        ref = str(row.get("ref", "") or "")
        alt = str(row.get("alt", "") or "")
        if len(ref) == 1 and len(alt) == 1:
            return "snp"
        if len(ref) > len(alt):
            return "deletion"
        if len(alt) > len(ref):
            return "insertion"
        indel_size = row.get("indel_size", np.nan)
        if pd.notna(indel_size) and float(indel_size) != 0:
            return "indel"

    return "indel"


def load_magi_baseline(
    path: Path = MAGI_BASELINE_FILE,
) -> Dict[str, Dict[str, pd.Series]]:
    """Load cached benign baseline stats used for MAGI z-scoring."""
    global _MAGI_BASELINE_CACHE

    if _MAGI_BASELINE_CACHE is not None:
        return _MAGI_BASELINE_CACHE
    if not path.exists():
        _MAGI_BASELINE_CACHE = {}
        return _MAGI_BASELINE_CACHE

    stats_df = pd.read_csv(path)
    baseline: Dict[str, Dict[str, pd.Series]] = {}
    for variant_type, subset in stats_df.groupby("variant_type"):
        key = str(variant_type).strip().lower()
        baseline[key] = {
            "mean": pd.Series(
                subset["mean"].astype(float).to_numpy(),
                index=subset["delta_col"].astype(str),
            ),
            "std": pd.Series(
                subset["std"].astype(float).to_numpy(),
                index=subset["delta_col"].astype(str),
            ),
        }

    _MAGI_BASELINE_CACHE = baseline
    return baseline


def compute_global_z_sum_log(
    row: pd.Series,
    baseline: Optional[Dict[str, Dict[str, pd.Series]]] = None,
    variant_type: Optional[str] = None,
) -> float:
    """Compute MAGI `Global_z_sum_log` from BED and BigWig delta columns."""
    baseline = baseline if baseline is not None else load_magi_baseline()
    if not baseline:
        return np.nan

    variant_key = _normalize_variant_type_key(
        variant_type if variant_type is not None else row.get("variant_type"),
        row=row,
    )
    lookup_order = [variant_key]
    if variant_key != "snp":
        lookup_order.append("indel")

    stats = None
    for key in lookup_order:
        stats = baseline.get(key)
        if stats:
            break
    if not stats:
        return np.nan

    delta_cols = [
        col
        for col in row.index
        if (col.startswith("D_BED_") or col.startswith("D_BW_"))
        and not pd.isna(row[col])
    ]
    if not delta_cols:
        return np.nan

    mean_series = stats["mean"]
    std_series = stats["std"]
    present_cols = [
        col
        for col in delta_cols
        if col in mean_series.index and col in std_series.index
    ]
    if not present_cols:
        return np.nan

    values = np.asarray([float(row[col]) for col in present_cols], dtype=np.float64)
    mu = mean_series.reindex(present_cols).to_numpy(dtype=np.float64)
    sigma = std_series.reindex(present_cols).to_numpy(dtype=np.float64)
    sigma = np.where(np.isfinite(sigma) & (sigma >= 1e-8), sigma, 1.0)

    z_scores = (values - mu) / sigma
    z_scores[~np.isfinite(z_scores)] = np.nan
    if not np.isfinite(z_scores).any():
        return np.nan

    return float(np.nansum(np.log1p(np.abs(z_scores))))


def extract_top_summary_signals(
    row: pd.Series,
    ranked: List[Dict[str, Any]],
    min_abs_threshold: float = 0.03,
    max_per_source: int = 5,
) -> Dict[str, List[Dict[str, Any]]]:
    """Extract compact top BED, BigWig, and MLM signals for the summary card."""
    summary = {"bed": [], "bigwig": [], "mlm": []}

    for item in ranked:
        if abs(float(item.get("delta", 0.0))) < min_abs_threshold:
            continue
        payload = {
            "label": item.get("display_name", item.get("track_id", "")),
            "track_id": item.get("track_id", ""),
            "delta": float(item.get("delta", np.nan)),
            "ref_val": item.get("ref_val", np.nan),
            "alt_val": item.get("alt_val", np.nan),
        }
        if item.get("track_type") == "BED" and len(summary["bed"]) < max_per_source:
            summary["bed"].append(payload)
        if (
            item.get("track_type") == "BigWig"
            and len(summary["bigwig"]) < max_per_source
        ):
            summary["bigwig"].append(payload)
        if (
            len(summary["bed"]) >= max_per_source
            and len(summary["bigwig"]) >= max_per_source
        ):
            break

    mlm_items: List[Dict[str, Any]] = []
    for column, label, signed in MLM_SIGNAL_SPECS:
        if column not in row.index or pd.isna(row[column]):
            continue
        value = float(row[column])
        magnitude = abs(value)
        if magnitude < min_abs_threshold:
            continue
        mlm_items.append(
            {
                "label": label,
                "column": column,
                "value": value,
                "magnitude": magnitude,
                "signed": signed,
            }
        )
    mlm_items.sort(key=lambda item: item["magnitude"], reverse=True)
    summary["mlm"] = mlm_items[:max_per_source]
    return summary


def get_top_features(
    row: pd.Series,
    bed_cols: List[str],
    bw_cols: Optional[List[str]] = None,
    k: int = 10,
) -> List[Tuple[str, float, str]]:
    """
    Extract top K disrupted features for a variant.

    Args:
        row: DataFrame row with delta columns
        bed_cols: BED delta column names
        bw_cols: BigWig delta column names (optional)
        k: Number of features to return

    Returns:
        List of (feature_name, delta, direction) tuples
        direction is 'GOF' (gain) or 'LOF' (loss)
    """
    deltas = {}

    # Accumulate BED deltas
    for col in bed_cols:
        if col in row.index and not pd.isna(row[col]):
            name = col.replace("D_BED_", "")
            deltas[name] = float(row[col])

    # Accumulate BigWig deltas (truncate long names)
    if bw_cols:
        for col in bw_cols:
            if col in row.index and not pd.isna(row[col]):
                name = col.replace("D_BW_", "")
                if len(name) > 40:
                    name = name[:37] + "..."
                deltas[name] = float(row[col])

    # Sort by absolute value
    top_items = sorted(deltas.items(), key=lambda x: abs(x[1]), reverse=True)[:k]

    # Add direction
    result = []
    for feat_name, delta_val in top_items:
        direction = "GOF" if delta_val > 0 else "LOF"
        result.append((feat_name, delta_val, direction))

    return result


def summarize_variant(
    row: pd.Series, bed_cols: List[str], bw_cols: Optional[List[str]] = None, k: int = 5
) -> MechanisticSummary:
    """
    Create a mechanistic summary for a variant.

    Args:
        row: DataFrame row
        bed_cols: BED delta column names
        bw_cols: BigWig delta column names
        k: Number of top features

    Returns:
        MechanisticSummary object
    """
    top_features = get_top_features(row, bed_cols, bw_cols, k=k)

    variant_id = f"{row.get('chrom', 'chr?')}:{row.get('pos', '?')} {row.get('ref', '?')}>{row.get('alt', '?')}"
    gene_name = row.get("gene_name", "")
    region_class = row.get("region_class", "OTHER")
    impact_bed = row.get("Impact_Score_BED", np.nan)
    llr = row.get("LLR", np.nan)

    return MechanisticSummary(
        variant_id=variant_id,
        gene_name=gene_name,
        region_class=region_class,
        top_features=top_features,
        impact_score_bed=impact_bed,
        llr=llr,
    )


# ============================================================================
# BIGWIG DESCRIPTION MAPPING
# ============================================================================
def get_top_bigwig_descriptions(
    row: pd.Series, metadata_df: pd.DataFrame, k: int = 5
) -> List[Tuple[str, float, str]]:
    """
    Get top BigWig track changes with human-readable descriptions.

    Args:
        row: DataFrame row with D_BW_* columns
        metadata_df: Track metadata with columns [file_id, tissue, assay, experiment_target]
        k: Number of tracks to return

    Returns:
        List of (description, delta, direction) tuples
        description format: "{tissue} | {assay} | {target}"
    """
    bw_cols = [c for c in row.index if c.startswith("D_BW_")]

    deltas = {}
    for col in bw_cols:
        if not pd.isna(row[col]):
            track_id = col.replace("D_BW_", "")
            delta = float(row[col])

            # Look up metadata
            meta = metadata_df[metadata_df["file_id"] == track_id]
            if not meta.empty:
                r = meta.iloc[0]
                tissue = r.get("tissue", "")
                assay = r.get("assay", "")
                target = r.get("experiment_target", "")
                desc = f"{tissue} | {assay}"
                if pd.notna(target) and str(target).strip():
                    desc += f" | {target}"
            else:
                desc = track_id

            deltas[desc] = delta

    # Sort by absolute value
    top_items = sorted(deltas.items(), key=lambda x: abs(x[1]), reverse=True)[:k]

    result = []
    for desc, delta in top_items:
        direction = "Gain" if delta > 0 else "Loss"
        result.append((desc, delta, direction))

    return result


def _clean_part(s: str) -> str:
    """Strip whitespace, trailing/leading commas, and collapse double commas."""
    s = s.strip().strip(",").strip()
    # collapse repeated commas with optional spaces
    import re
    s = re.sub(r",\s*,", ",", s)
    return s


_EMPTY_TARGETS = {"", "nan", "none", "n/a", "na", "null"}


def _resolve_bw_name(
    track_id: str,
    metadata_df: Optional[pd.DataFrame] = None,
    metadata_dict: Optional[Dict[str, Dict[str, str]]] = None,
    max_len: int = 55,
) -> str:
    """Resolve BigWig track_id → human-readable display name."""
    if metadata_dict:
        meta = metadata_dict.get(track_id)
        if meta:
            parts = [
                _clean_part(p)
                for p in [
                    meta.get("tissue", ""),
                    meta.get("assay", ""),
                    meta.get("target", ""),
                ]
                if _clean_part(p) and p.strip().lower() not in _EMPTY_TARGETS
            ]
            if parts:
                name = " | ".join(parts)
                return name[:max_len] if len(name) > max_len else name

    if metadata_df is not None:
        rows = metadata_df[metadata_df["file_id"] == track_id]
        if not rows.empty:
            r = rows.iloc[0]
            parts = [
                _clean_part(str(p))
                for p in [
                    r.get("tissue", ""),
                    r.get("assay", ""),
                    r.get("experiment_target", ""),
                ]
                if pd.notna(p) and _clean_part(str(p)) and str(p).strip().lower() not in _EMPTY_TARGETS
            ]
            if parts:
                name = " | ".join(parts)
                return name[:max_len] if len(name) > max_len else name

    return track_id[:40]


def rank_top_disrupted_tracks(
    row: pd.Series,
    bed_names: List[str],
    bw_names: Optional[List[str]] = None,
    metadata_df: Optional[pd.DataFrame] = None,
    metadata_dict: Optional[Dict[str, Dict[str, str]]] = None,
    top_k: Optional[int] = None,
) -> List[Dict]:
    """
    Single-source ranking of the most disrupted tracks across BED + BigWig.

    Returns an ordered list (by |delta| descending) of dicts:
        {track_id, display_name, delta, track_type ("BED"/"BigWig"),
         ref_val, alt_val}
    """
    items: List[Dict] = []

    for name in bed_names:
        d_col = f"D_BED_{name}"
        r_col = f"REF_BED_{name}"
        if d_col not in row.index or pd.isna(row[d_col]):
            continue
        delta = float(row[d_col])
        ref_v = (
            float(row[r_col])
            if r_col in row.index and not pd.isna(row[r_col])
            else np.nan
        )
        items.append(
            {
                "track_id": name,
                "display_name": name,
                "delta": delta,
                "track_type": "BED",
                "ref_val": ref_v,
                "alt_val": ref_v + delta if not np.isnan(ref_v) else np.nan,
            }
        )

    if bw_names:
        for track_id in bw_names:
            d_col = f"D_BW_{track_id}"
            r_col = f"REF_BW_{track_id}"
            if d_col not in row.index or pd.isna(row[d_col]):
                continue
            delta = float(row[d_col])
            ref_v = (
                float(row[r_col])
                if r_col in row.index and not pd.isna(row[r_col])
                else np.nan
            )
            items.append(
                {
                    "track_id": track_id,
                    "display_name": _resolve_bw_name(
                        track_id, metadata_df, metadata_dict
                    ),
                    "delta": delta,
                    "track_type": "BigWig",
                    "ref_val": ref_v,
                    "alt_val": ref_v + delta if not np.isnan(ref_v) else np.nan,
                }
            )

    items.sort(key=lambda x: abs(x["delta"]), reverse=True)
    for idx, item in enumerate(items, start=1):
        item["abs_delta"] = abs(item["delta"])
        item["rank"] = idx

    if top_k is None:
        return items
    return items[:top_k]


def build_top_track_table(
    ranked: List[Dict],
    max_rows: int = 15,
    min_rows_by_type: Optional[Dict[str, int]] = None,
) -> pd.DataFrame:
    """Build a single merged top-tracks table from the unified ranking."""
    if not ranked:
        return pd.DataFrame(
            [
                {
                    "Track": "No disruptions detected",
                    "Type": "",
                    "REF": "",
                    "ALT": "",
                    "Δ": "",
                    "Direction": "",
                }
            ]
        )

    min_rows_by_type = min_rows_by_type or {"BED": 4}
    max_rows = max(max_rows, sum(min_rows_by_type.values()))

    forced_keys = set()
    for track_type, min_rows in min_rows_by_type.items():
        count = 0
        for item in ranked:
            if item["track_type"] != track_type:
                continue
            forced_keys.add((item["track_type"], item["track_id"]))
            count += 1
            if count >= min_rows:
                break

    selected = []
    selected_keys = set()
    for item in ranked:
        key = (item["track_type"], item["track_id"])
        if key in forced_keys and key not in selected_keys:
            selected.append(item)
            selected_keys.add(key)

    for item in ranked:
        key = (item["track_type"], item["track_id"])
        if key in selected_keys:
            continue
        if len(selected) >= max_rows:
            break
        selected.append(item)
        selected_keys.add(key)

    selected.sort(key=lambda item: item.get("rank", 10**9))

    rows = []
    for item in selected[:max_rows]:
        ref_v = item["ref_val"]
        alt_v = item["alt_val"]
        rows.append(
            {
                "Track": item["display_name"],
                "Type": item["track_type"],
                "REF": f"{ref_v:.4f}" if not np.isnan(ref_v) else "N/A",
                "ALT": f"{alt_v:.4f}" if not np.isnan(alt_v) else "N/A",
                "Δ": f"{item['delta']:+.4f}",
                "Direction": "Gain" if item["delta"] > 0 else "Loss",
            }
        )
    return pd.DataFrame(rows)


# ============================================================================
# VISUALIZATION
# ============================================================================
def make_fingerprint_plot(
    ranked: List[Dict],
    top_k: int = 15,
    figsize: Tuple[float, float] = (10, 6),
) -> plt.Figure:
    """
    Horizontal bar chart of top disrupted features.

    Args:
        ranked: Pre-ranked list from rank_top_disrupted_tracks().
        top_k: Number of features to show (capped by len(ranked)).
        figsize: Figure size.
    """
    items = ranked[:top_k]

    if not items:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5, 0.5, "No feature data available", ha="center", va="center", fontsize=14
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        return fig

    names = [it["display_name"] for it in items]
    vals = [it["delta"] for it in items]
    colors = ["#d73027" if v > 0 else "#4575b4" for v in vals]

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(names))
    ax.barh(y_pos, vals, color=colors, alpha=0.8, edgecolor="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=1.5, linestyle="-", alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Δ Probability (Alt − Ref)", fontsize=11, fontweight="bold")
    ax.set_title(
        "Top Disrupted Genomic Features", fontsize=13, fontweight="bold", pad=15
    )
    ax.grid(True, axis="x", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)

    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#d73027", alpha=0.8, label="Gain of Function"),
        Patch(facecolor="#4575b4", alpha=0.8, label="Loss of Function"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", frameon=True, fontsize=9)

    plt.tight_layout()
    return fig


def format_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract clinically relevant columns for display.

    Args:
        df: Full results DataFrame

    Returns:
        DataFrame with selected columns for display
    """
    display_cols = ["chrom", "pos", "ref", "alt"]

    # Add annotation if available
    if "gene_name" in df.columns:
        display_cols.append("gene_name")
    if "region_class" in df.columns:
        display_cols.append("region_class")

    # Add impact scores
    if "Impact_Score_BED" in df.columns:
        display_cols.append("Impact_Score_BED")
    if "Impact_Score_BW" in df.columns:
        display_cols.append("Impact_Score_BW")
    if "Global_z_sum_log" in df.columns:
        display_cols.append("Global_z_sum_log")

    # Add MLM features
    for col in ["LLR", "MLM_KL_mean", "MLM_KL_max"]:
        if col in df.columns:
            display_cols.append(col)

    # Add indel size if present
    if "indel_size" in df.columns:
        display_cols.append("indel_size")

    # Filter to available columns
    available = [c for c in display_cols if c in df.columns]

    return df[available].copy()
