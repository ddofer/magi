"""Threshold-gated signal extraction (matches analysis_v2.ipynb legacy outputs)."""

from __future__ import annotations

import heapq
from typing import Optional

import numpy as np
import pandas as pd
import heapq

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# helpers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

################################################################################
# thresholds access getters
################################################################################

# --- helpers (fast + safe) ---
def _thr_row(thresholds, track_id, ref="benign", stat="std"):
    # returns Series with ["lo","hi"]
    return thresholds.loc[(track_id, ref, stat), ["lo", "hi"]]

def thr_hi(thresholds, track_id, ref="benign", stat="std"):
    return float(_thr_row(thresholds, track_id, ref, stat)["hi"])

def thr_lo(thresholds, track_id, ref="benign", stat="std"):
    return float(_thr_row(thresholds, track_id, ref, stat)["lo"])

def thr_mag(thresholds, track_id, ref="benign", stat="std"):
    r = _thr_row(thresholds, track_id, ref, stat)
    lo = float(r["lo"]); hi = float(r["hi"])
    return max(abs(lo), abs(hi))

################################################################################
# build threshold maps
################################################################################

def build_bed_threshold_maps(thresholds, prefix="D_BED_", ref="benign", stat="std"):
    """
    thresholds: DataFrame with MultiIndex (track, reference, stat) and columns lo/hi
    Returns dicts keyed by feature (without prefix): mag, hi, lo
    """
    # slice once: all BED rows for (ref, stat)
    sub = thresholds.xs((ref, stat), level=("reference", "stat"), drop_level=False)

    # keep only BED tracks
    sub = sub[sub.index.get_level_values("track").str.startswith(prefix)]

    # index track_id -> lo/hi as float numpy arrays
    track_ids = sub.index.get_level_values("track").to_numpy()
    lo = sub["lo"].to_numpy(dtype=float)
    hi = sub["hi"].to_numpy(dtype=float)
    mag = np.maximum(np.abs(lo), np.abs(hi))

    # convert track_id -> feature
    feats = np.char.replace(track_ids.astype(str), prefix, "")

    thr_mag = dict(zip(feats, mag))
    thr_hi  = dict(zip(feats, hi))
    thr_lo  = dict(zip(feats, lo))
    return thr_mag, thr_hi, thr_lo

def build_bw_threshold_maps(thresholds, prefix="D_BW_", ref="benign", stat="std"):
    sub = thresholds.xs((ref, stat), level=("reference", "stat"), drop_level=False)
    sub = sub[sub.index.get_level_values("track").str.startswith(prefix)]

    track_ids = sub.index.get_level_values("track").to_numpy()
    lo = sub["lo"].to_numpy(dtype=float)
    hi = sub["hi"].to_numpy(dtype=float)
    mag = np.maximum(np.abs(lo), np.abs(hi))

    tids = np.char.replace(track_ids.astype(str), prefix, "")

    thr_mag = dict(zip(tids, mag))
    thr_hi  = dict(zip(tids, hi))
    thr_lo  = dict(zip(tids, lo))
    return thr_mag, thr_hi, thr_lo


################################################################################
# top k selectors
################################################################################

def select_bed_topk(bed_data, bed_thr_mag, bed_thr_hi, bed_thr_lo, k=10, eps=0.01):
    cand_abs = []
    cand_gain = []
    cand_loss = []

    for x in bed_data:
        f = x["feature"]
        d = x["delta"]

        # lookup thresholds (defaults chosen so missing features are ignored)
        mag = bed_thr_mag.get(f, np.inf)
        hi  = bed_thr_hi.get(f, np.inf)
        lo  = bed_thr_lo.get(f, -np.inf)

        # abs gate
        if abs(d) > eps and abs(d) > mag:
            cand_abs.append(x)

        # gains gate
        if d > eps and d > hi:
            cand_gain.append(x)

        # losses gate
        if d < -eps and d < lo:
            cand_loss.append(x)

    bed_top_abs = heapq.nlargest(k, cand_abs, key=lambda x: abs(x["delta"]))
    bed_top_gains = heapq.nlargest(k, cand_gain, key=lambda x: x["delta"])
    bed_top_losses = heapq.nsmallest(k, cand_loss, key=lambda x: x["delta"])
    return bed_top_abs, bed_top_gains, bed_top_losses


def select_bw_topk(bw_data, bw_thr_mag, bw_thr_hi, bw_thr_lo, k=10, eps=0.01):
    cand_abs = []
    cand_gain = []
    cand_loss = []

    for x in bw_data:
        tid = x["track_id"]
        d = x["delta"]

        mag = bw_thr_mag.get(tid, np.inf)
        hi  = bw_thr_hi.get(tid, np.inf)
        lo  = bw_thr_lo.get(tid, -np.inf)

        if abs(d) > eps and abs(d) > mag:
            cand_abs.append(x)
        if d > eps and d > hi:
            cand_gain.append(x)
        if d < -eps and d < lo:
            cand_loss.append(x)

    bw_top_abs = heapq.nlargest(k, cand_abs, key=lambda x: abs(x["delta"]))
    bw_top_gains = heapq.nlargest(k, cand_gain, key=lambda x: x["delta"])
    bw_top_losses = heapq.nsmallest(k, cand_loss, key=lambda x: x["delta"])
    return bw_top_abs, bw_top_gains, bw_top_losses

def extract_signals_comprehensive(
    row,
    delta_cols,
    metadata_lookup=None,
    k=10,
    variant_type='snp',
    bed_thr_mag=None, bed_thr_hi=None, bed_thr_lo=None,
    bw_thr_mag=None,  bw_thr_hi=None,  bw_thr_lo=None
    ):   
    """
    Extract comprehensive signal information including:
    - BED signals (delta + reference)
    - MLM signals (all relevant columns)
    - BW signals (delta + reference + metadata annotation)
    
    Args:
        row: DataFrame row with all columns
        delta_cols: List of delta column names (D_BED_*, D_BW_*)
        metadata_df: Tracks metadata DataFrame with file_id -> biosample/tissue/assay mapping
        k: Number of top signals to extract per category
        variant_type: 'snp' or 'indel' (affects which MLM columns to extract)
    
    Returns:
        dict with structured signal information
    """
    
    # =========================================================================
    # 1. SEPARATE COLUMN TYPES
    # =========================================================================
    
    bed_delta_cols = [c for c in delta_cols if c.startswith('D_BED_')]
    bw_delta_cols = [c for c in delta_cols if c.startswith('D_BW_')]
    
    
    # =========================================================================
    # 2. EXTRACT BED SIGNALS (Delta + Reference)
    # =========================================================================
    
    bed_data = []
    for col in bed_delta_cols:
        delta_val = row.get(col, 0)
        if pd.isna(delta_val):
            delta_val = 0
        
        # Get corresponding reference value
        ref_col = col.replace('D_BED_', 'REF_BED_')
        ref_val = row.get(ref_col, np.nan)
        
        feature_name = col.replace('D_BED_', '')
        
        bed_data.append({
            'feature': feature_name,
            'delta': delta_val,
            'ref': ref_val if not pd.isna(ref_val) else None,
            'type': 'BED'
        })
    
    # Sort by absolute delta
    # Fast top-k selection using precomputed BED threshold maps
    bed_top_abs, bed_top_gains, bed_top_losses = select_bed_topk(
        bed_data,
        bed_thr_mag=bed_thr_mag or {},
        bed_thr_hi=bed_thr_hi or {},
        bed_thr_lo=bed_thr_lo or {},
        k=k,
        eps=0.01
    )
    
    # =========================================================================
    # 3. EXTRACT BW SIGNALS (Delta + Reference + Metadata)
    # =========================================================================
    
    # Build metadata lookup if available
    metadata_lookup = metadata_lookup or {}
    
    bw_data = []
    for col in bw_delta_cols:
        delta_val = row.get(col, 0)
        if pd.isna(delta_val):
            delta_val = 0
        
        # Get corresponding reference value
        ref_col = col.replace('D_BW_', 'REF_BW_')
        ref_val = row.get(ref_col, np.nan)
        
        # Extract track ID
        track_id = col.replace('D_BW_', '')
        
        # Get metadata annotation
        meta = metadata_lookup.get(track_id, {})
        
        # Create human-readable description
        if meta:
            description = f"{meta.get('tissue', 'Unknown')}|{meta.get('biosample_type', '')}|{meta.get('assay', '')}"
            if meta.get('experiment_target'):
                description += f"|{meta['experiment_target']}"
        else:
            description = track_id  # Fallback to ID if no metadata
        
        bw_data.append({
            'track_id': track_id,
            'description': description,
            'delta': delta_val,
            'ref': ref_val if not pd.isna(ref_val) else None,
            'type': 'BW',
            'metadata': meta
        })
    
    # Sort by absolute delta
    # Fast top-k selection using precomputed BW threshold maps
    bw_top_abs, bw_top_gains, bw_top_losses = select_bw_topk(
        bw_data,
        bw_thr_mag=bw_thr_mag or {},
        bw_thr_hi=bw_thr_hi or {},
        bw_thr_lo=bw_thr_lo or {},
        k=k,
        eps=0.01
    )
    
    # =========================================================================
    # 4. EXTRACT MLM SIGNALS
    # =========================================================================
    
    # Define MLM columns based on variant type
    if variant_type == 'snp':
        mlm_cols = ['REF_5mer', 'ALT_5mer', 'LLR', 'MLM_Prior', 'MLM_Delta']
    else:  # indel
        mlm_cols = [
            'LLR', 'MLM_Prior', 'MLM_Delta', 'REF_5mer', 'ALT_5mer',
            'MLM_KL_mean', 'MLM_KL_max', 'MLM_logprob_ref', 'MLM_logprob_alt',
            'MLM_logprob_delta', 'EMB_cosine_dist', 'EMB_l2_dist',
            'EMB_max_pos_dist', 'EMB_mean_pos_dist'
        ]
    
    mlm_data = {}
    for col in mlm_cols:
        val = row.get(col, np.nan)
        mlm_data[col] = val if not pd.isna(val) else None
    
    # =========================================================================
    # 5. FORMAT OUTPUT - AGGREGATE TOP BED, BW, MLM SIGNALS
    # =========================================================================
    
    # --- String representations for DataFrame columns ---
    
    # BED: "feature=delta(ref=X); ..."
    def format_bed_signal(sig):
        ref_str = f"(ref={sig['ref']:.3f})" if sig['ref'] is not None else ""
        return f"{sig['feature']}={sig['delta']:.3f}{ref_str}"
    
    bed_abs_str = "; ".join([format_bed_signal(s) for s in bed_top_abs]) if bed_top_abs else "None"
    bed_gain_str = "; ".join([format_bed_signal(s) for s in bed_top_gains]) if bed_top_gains else "None"
    bed_loss_str = "; ".join([format_bed_signal(s) for s in bed_top_losses]) if bed_top_losses else "None"
    
    # BW: "description=delta(ref=X); ..."
    def format_bw_signal(sig):
        ref_str = f"(ref={sig['ref']:.3f})" if sig['ref'] is not None else ""
        return f"{sig['description']}={sig['delta']:.3f}{ref_str}"
    
    bw_abs_str = "; ".join([format_bw_signal(s) for s in bw_top_abs]) if bw_top_abs else "None"
    bw_gain_str = "; ".join([format_bw_signal(s) for s in bw_top_gains]) if bw_top_gains else "None"
    bw_loss_str = "; ".join([format_bw_signal(s) for s in bw_top_losses]) if bw_top_losses else "None"
    
    # MLM: Key metrics as string
    mlm_summary_parts = []
    if mlm_data.get('LLR') is not None:
        mlm_summary_parts.append(f"LLR={mlm_data['LLR']:.3f}")
    if mlm_data.get('MLM_Delta') is not None:
        mlm_summary_parts.append(f"MLM_Delta={mlm_data['MLM_Delta']:.3f}")
    if mlm_data.get('MLM_Prior') is not None:
        mlm_summary_parts.append(f"MLM_Prior={mlm_data['MLM_Prior']:.3f}")
    if variant_type == 'indel':
        if mlm_data.get('EMB_cosine_dist') is not None:
            mlm_summary_parts.append(f"EMB_cos={mlm_data['EMB_cosine_dist']:.3f}")
        if mlm_data.get('MLM_KL_max') is not None:
            mlm_summary_parts.append(f"KL_max={mlm_data['MLM_KL_max']:.3f}")
        if mlm_data.get('MLM_logprob_delta') is not None:
             mlm_summary_parts.append(f"LogProb_Delta={mlm_data['MLM_logprob_delta']:.1f}")
    
    mlm_summary_str = "; ".join(mlm_summary_parts) if mlm_summary_parts else "None"
    
    # =========================================================================
    # 6. RETURN STRUCTURED RESULT
    # =========================================================================
    
    return {
        # --- String columns for DataFrame storage ---
        'BED_Top_Abs': bed_abs_str,
        'BED_Top_Gains': bed_gain_str,
        'BED_Top_Losses': bed_loss_str,
        'BW_Top_Abs': bw_abs_str,
        'BW_Top_Gains': bw_gain_str,
        'BW_Top_Losses': bw_loss_str,
        'MLM_Summary': mlm_summary_str,
        
        # --- Individual MLM columns (for analysis) ---
        **{f'MLM_{k}': v for k, v in mlm_data.items()},
        
        # --- Aggregate metrics ---
        'BED_Max_Abs_Delta': max([abs(s['delta']) for s in bed_data]) if bed_data else 0,
        'BW_Max_Abs_Delta': max([abs(s['delta']) for s in bw_data]) if bw_data else 0,
        'BED_N_Strong_Signals': sum(1 for s in bed_data if abs(s['delta']) > 0.3),
        'BW_N_Strong_Signals': sum(1 for s in bw_data if abs(s['delta']) > 0.3),
    }


def _build_metadata_lookup(metadata_df: Optional[pd.DataFrame]) -> dict:
    if metadata_df is None:
        return {}
    lookup = {}
    for _, meta_row in metadata_df.iterrows():
        file_id = meta_row.get("file_id", "")
        lookup[file_id] = {
            "biosample_type": meta_row.get("biosample_type", ""),
            "tissue": meta_row.get("tissue", ""),
            "assay": meta_row.get("assay", ""),
            "experiment_target": meta_row.get("experiment_target", ""),
            "dataset": meta_row.get("dataset", ""),
        }
    return lookup


def _load_threshold_index(thresholds_path: str) -> Optional[pd.DataFrame]:
    try:
        df_thr = pd.read_csv(thresholds_path)
        return df_thr.set_index(["track", "reference", "stat"]).sort_index()
    except OSError as e:
        print(f"Warning: Could not load thresholds ({e})")
        return None


def extract_signals_with_thresholds(
    variant_df: pd.DataFrame,
    metadata_df: Optional[pd.DataFrame],
    thresholds_path: str,
    k: int = 5,
    variant_type: Optional[str] = None,
) -> pd.DataFrame:
    """Row-wise extraction using benign/std thresholds (legacy analysis_v2 path)."""
    if variant_type is None:
        variant_type = "indel" if "EMB_cosine_dist" in variant_df.columns else "snp"

    metadata_lookup = _build_metadata_lookup(metadata_df)
    thresholds = _load_threshold_index(thresholds_path)

    bed_thr_mag = bed_thr_hi = bed_thr_lo = {}
    bw_thr_mag = bw_thr_hi = bw_thr_lo = {}
    if thresholds is not None:
        bed_thr_mag, bed_thr_hi, bed_thr_lo = build_bed_threshold_maps(thresholds)
        bw_thr_mag, bw_thr_hi, bw_thr_lo = build_bw_threshold_maps(thresholds)

    delta_cols = [
        c for c in variant_df.columns if c.startswith("D_BED_") or c.startswith("D_BW_")
    ]
    out = variant_df.copy()
    for col in delta_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    print(f"Extracting signals (threshold-gated, k={k}, variant_type={variant_type})...")
    signal_results = out.apply(
        lambda row: extract_signals_comprehensive(
            row,
            delta_cols,
            metadata_lookup=metadata_lookup,
            k=k,
            variant_type=variant_type,
            bed_thr_mag=bed_thr_mag,
            bed_thr_hi=bed_thr_hi,
            bed_thr_lo=bed_thr_lo,
            bw_thr_mag=bw_thr_mag,
            bw_thr_hi=bw_thr_hi,
            bw_thr_lo=bw_thr_lo,
        ),
        axis=1,
        result_type="expand",
    )
    return pd.concat([out, signal_results], axis=1)
