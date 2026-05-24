"""Vectorized BED / BW / MLM signal extraction for LLM prompts."""

from __future__ import annotations

# ──────────────────────────────────────────────────────────────────────────────
# FAST VECTORIZED SIGNAL EXTRACTION FOR LLM
# ──────────────────────────────────────────────────────────────────────────────
#
# Replaces row-wise df.apply() with vectorized NumPy operations.
# Expected speedup: 50-100x on numeric aggregation, 5-10x on string formatting.
# ──────────────────────────────────────────────────────────────────────────────

"""
 Extract comprehensive signal information including:
    - BED signals (delta + reference)
    - MLM signals (all relevant columns)
    - BW signals (delta + reference + metadata annotation)
"""
import pandas as pd
import numpy as np
from typing import Optional, Dict, List, Tuple


# =============================================================================
# CORE VECTORIZED HELPERS
# =============================================================================

def _vectorized_top_k_indices(matrix: np.ndarray, k: int) -> np.ndarray:
    """
    Get top-k column indices per row by absolute value.
    Uses argpartition (O(n)) instead of full argsort (O(n log n)).

    Args:
        matrix: (n_rows, n_cols) array of values (already filled, no NaN)
        k: number of top signals to extract

    Returns:
        (n_rows, k) array of column indices, sorted by descending |value|
    """
    n_rows, n_cols = matrix.shape
    k_actual = min(k, n_cols)
    abs_matrix = np.abs(matrix)

    if k_actual >= n_cols:
        # Just argsort the whole thing
        top_indices = np.argsort(-abs_matrix, axis=1)[:, :k_actual]
        return top_indices

    # argpartition: O(n) to get top-k (unsorted)
    top_indices = np.argpartition(abs_matrix, -k_actual, axis=1)[:, -k_actual:]

    # Sort just the k elements per row (O(k log k) per row)
    rows = np.arange(n_rows)[:, None]
    top_vals = abs_matrix[rows, top_indices]
    sort_within = np.argsort(-top_vals, axis=1)
    top_indices_sorted = np.take_along_axis(top_indices, sort_within, axis=1)

    return top_indices_sorted


def _vectorized_top_k_signed(matrix: np.ndarray, k: int, direction: str) -> np.ndarray:
    """
    Get top-k column indices for gains (direction='gain') or losses (direction='loss').

    For gains: top-k by descending value where value > 0.01
    For losses: top-k by ascending value where value < -0.01

    Returns:
        (n_rows, k) array of column indices. Padded with -1 for rows with fewer signals.
    """
    n_rows, n_cols = matrix.shape
    k_actual = min(k, n_cols)

    if direction == 'gain':
        # Mask non-gains, sort descending
        work = matrix.copy()
        work[work <= 0.01] = -np.inf  # push non-gains to bottom
        top_indices = np.argsort(-work, axis=1)[:, :k_actual]
        # Mark invalid entries
        rows = np.arange(n_rows)[:, None]
        vals = matrix[rows, top_indices]
        top_indices[vals <= 0.01] = -1
    else:  # loss
        work = matrix.copy()
        work[work >= -0.01] = np.inf  # push non-losses to bottom
        top_indices = np.argsort(work, axis=1)[:, :k_actual]
        rows = np.arange(n_rows)[:, None]
        vals = matrix[rows, top_indices]
        top_indices[vals >= -0.01] = -1

    return top_indices


# =============================================================================
# STRING FORMATTING (BULK)
# =============================================================================

def _format_signal_strings(
    matrix: np.ndarray,
    ref_matrix: np.ndarray,
    top_k_indices: np.ndarray,
    feature_names: np.ndarray,
    include_ref: bool = True
) -> List[str]:
    """
    Build formatted signal strings for all rows at once.

    Format: "feature=delta(ref=X); feature2=delta2(ref=Y); ..."

    Args:
        matrix: (n_rows, n_cols) delta values
        ref_matrix: (n_rows, n_cols) reference values (may contain NaN)
        top_k_indices: (n_rows, k) indices into columns. -1 means skip.
        feature_names: (n_cols,) array of feature name strings
        include_ref: whether to include reference values

    Returns:
        List of formatted strings, one per row
    """
    n_rows, k = top_k_indices.shape
    results = []

    for i in range(n_rows):
        parts = []
        for j in range(k):
            col_idx = top_k_indices[i, j]
            if col_idx == -1:
                continue
            d = matrix[i, col_idx]
            if d == 0 and not include_ref:
                continue

            part = f"{feature_names[col_idx]}={d:.3f}"
            if include_ref:
                r = ref_matrix[i, col_idx]
                if not np.isnan(r):
                    part += f"(ref={r:.3f})"
            parts.append(part)

        results.append("; ".join(parts) if parts else "None")

    return results


def _format_bw_signal_strings(
    matrix: np.ndarray,
    ref_matrix: np.ndarray,
    top_k_indices: np.ndarray,
    descriptions: np.ndarray,
) -> List[str]:
    """
    Build formatted BW signal strings with metadata descriptions.
    Format: "tissue|biosample|assay=delta(ref=X); ..."
    """
    n_rows, k = top_k_indices.shape
    results = []

    for i in range(n_rows):
        parts = []
        for j in range(k):
            col_idx = top_k_indices[i, j]
            if col_idx == -1:
                continue
            d = matrix[i, col_idx]

            part = f"{descriptions[col_idx]}={d:.3f}"
            r = ref_matrix[i, col_idx]
            if not np.isnan(r):
                part += f"(ref={r:.3f})"
            parts.append(part)

        results.append("; ".join(parts) if parts else "None")

    return results


# =============================================================================
# METADATA HELPERS
# =============================================================================

def _build_metadata_lookup(metadata_df: Optional[pd.DataFrame]) -> Dict[str, Dict]:
    """Build file_id -> metadata dict."""
    if metadata_df is None:
        return {}

    lookup = {}
    for _, row in metadata_df.iterrows():
        file_id = row.get('file_id', '')
        lookup[file_id] = {
            'biosample_type': row.get('biosample_type', ''),
            'tissue': row.get('tissue', ''),
            'assay': row.get('assay', ''),
            'experiment_target': row.get('experiment_target', ''),
            'dataset': row.get('dataset', '')
        }
    return lookup


def _build_bw_descriptions(bw_delta_cols: List[str], metadata_lookup: Dict) -> np.ndarray:
    """Build description array for BW tracks."""
    descriptions = []
    for col in bw_delta_cols:
        track_id = col.replace('D_BW_', '')
        meta = metadata_lookup.get(track_id, {})
        if meta:
            desc = f"{meta.get('tissue', 'Unknown')}|{meta.get('biosample_type', '')}|{meta.get('assay', '')}"
            if meta.get('experiment_target'):
                desc += f"|{meta['experiment_target']}"
        else:
            desc = track_id
        descriptions.append(desc)
    return np.array(descriptions)


# =============================================================================
# MLM EXTRACTION (VECTORIZED)
# =============================================================================

def _extract_mlm_columns(df: pd.DataFrame, variant_type: str) -> Tuple[pd.DataFrame, List[str]]:
    """Extract MLM columns and return as sub-DataFrame."""
    if variant_type == 'snp':
        mlm_cols = ['REF_5mer', 'ALT_5mer', 'LLR', 'MLM_Prior', 'MLM_Delta']
    else:
        mlm_cols = [
            'LLR', 'MLM_Prior', 'MLM_Delta', 'REF_5mer', 'ALT_5mer',
            'MLM_KL_mean', 'MLM_KL_max', 'MLM_logprob_ref', 'MLM_logprob_alt',
            'MLM_logprob_delta', 'EMB_cosine_dist', 'EMB_l2_dist',
            'EMB_max_pos_dist', 'EMB_mean_pos_dist'
        ]

    existing_cols = [c for c in mlm_cols if c in df.columns]
    return df[existing_cols] if existing_cols else pd.DataFrame(index=df.index), existing_cols


def _format_mlm_summary_strings(df: pd.DataFrame, variant_type: str) -> List[str]:
    """Vectorized MLM summary string construction."""
    n = len(df)
    results = [""] * n

    # Define which columns go into the summary and their format
    core_cols = [
        ('LLR', 'LLR', '.3f'),
        ('MLM_Delta', 'MLM_Delta', '.3f'),
        ('MLM_Prior', 'MLM_Prior', '.3f'),
        ('MLM_logprob_delta', 'LogProb_Delta', '.1f'),
        ('MLM_logprob_ref', 'LogProb_Ref', '.1f'),

    ]

    indel_cols = [
        ('EMB_cosine_dist', 'EMB_cos', '.3f'),
        ('MLM_KL_max', 'KL_max', '.3f'),
        ('MLM_KL_mean', 'KL_mean', '.3f'),
    ]

    cols_to_format = core_cols + (indel_cols if variant_type == 'indel' else [])

    # Pre-extract arrays
    arrays = {}
    for col_name, _, _ in cols_to_format:
        if col_name in df.columns:
            arrays[col_name] = df[col_name].values

    for i in range(n):
        parts = []
        for col_name, label, fmt in cols_to_format:
            if col_name in arrays:
                val = arrays[col_name][i]
                if not pd.isna(val):
                    parts.append(f"{label}={val:{fmt}}")
        results[i] = "; ".join(parts) if parts else "None"

    return results


# =============================================================================
# MAIN FAST EXTRACTION
# =============================================================================

def extract_signals_comprehensive_fast(
    df: pd.DataFrame,
    delta_cols: List[str],
    metadata_df: Optional[pd.DataFrame] = None,
    k: int = 5,
    variant_type: str = 'snp'
) -> pd.DataFrame:
    """
    Vectorized replacement for row-wise extract_signals_comprehensive.
    Operates on entire DataFrame at once using NumPy.

    Args:
        df: Full variant DataFrame
        delta_cols: List of delta column names (D_BED_*, D_BW_*)
        metadata_df: Track metadata DataFrame
        k: Number of top signals per category
        variant_type: 'snp' or 'indel'

    Returns:
        DataFrame with signal columns to join back
    """
    n_rows = len(df)
    result_df = pd.DataFrame(index=df.index)

    # -----------------------------------------------------------------
    # Separate column types
    # -----------------------------------------------------------------
    bed_delta_cols = [c for c in delta_cols if c.startswith('D_BED_')]
    bw_delta_cols = [c for c in delta_cols if c.startswith('D_BW_')]

    # -----------------------------------------------------------------
    # BED SIGNALS
    # -----------------------------------------------------------------
    if bed_delta_cols:
        bed_feature_names = np.array([c.replace('D_BED_', '') for c in bed_delta_cols])
        bed_ref_cols = [c.replace('D_BED_', 'REF_BED_') for c in bed_delta_cols]
        # Ensure ref cols exist
        bed_ref_cols_exist = [c for c in bed_ref_cols if c in df.columns]

        # Extract matrices
        bed_matrix = df[bed_delta_cols].fillna(0).values.astype(np.float64)

        # Build ref matrix (align columns)
        if bed_ref_cols_exist:
            bed_ref_matrix = np.full_like(bed_matrix, np.nan)
            for i, delta_col in enumerate(bed_delta_cols):
                ref_col = delta_col.replace('D_BED_', 'REF_BED_')
                if ref_col in df.columns:
                    bed_ref_matrix[:, i] = pd.to_numeric(df[ref_col], errors='coerce').values
        else:
            bed_ref_matrix = np.full_like(bed_matrix, np.nan)

        # Aggregate metrics (fully vectorized — no Python loops)
        abs_bed = np.abs(bed_matrix)
        result_df['BED_Max_Abs_Delta'] = abs_bed.max(axis=1)
        result_df['BED_N_Strong_Signals'] = (abs_bed > 0.3).sum(axis=1)

        # Top-K indices
        bed_top_abs_idx = _vectorized_top_k_indices(bed_matrix, k)
        bed_top_gain_idx = _vectorized_top_k_signed(bed_matrix, k, 'gain')
        bed_top_loss_idx = _vectorized_top_k_signed(bed_matrix, k, 'loss')

        # Format strings
        result_df['BED_Top_Abs'] = _format_signal_strings(
            bed_matrix, bed_ref_matrix, bed_top_abs_idx, bed_feature_names
        )
        result_df['BED_Top_Gains'] = _format_signal_strings(
            bed_matrix, bed_ref_matrix, bed_top_gain_idx, bed_feature_names
        )
        result_df['BED_Top_Losses'] = _format_signal_strings(
            bed_matrix, bed_ref_matrix, bed_top_loss_idx, bed_feature_names
        )
    else:
        result_df['BED_Max_Abs_Delta'] = 0
        result_df['BED_N_Strong_Signals'] = 0
        result_df['BED_Top_Abs'] = "None"
        result_df['BED_Top_Gains'] = "None"
        result_df['BED_Top_Losses'] = "None"

    # -----------------------------------------------------------------
    # BW SIGNALS (with metadata)
    # -----------------------------------------------------------------
    if bw_delta_cols:
        metadata_lookup = _build_metadata_lookup(metadata_df)
        bw_descriptions = _build_bw_descriptions(bw_delta_cols, metadata_lookup)
        bw_ref_cols = [c.replace('D_BW_', 'REF_BW_') for c in bw_delta_cols]

        bw_matrix = df[bw_delta_cols].fillna(0).values.astype(np.float64)

        bw_ref_matrix = np.full_like(bw_matrix, np.nan)
        for i, delta_col in enumerate(bw_delta_cols):
            ref_col = delta_col.replace('D_BW_', 'REF_BW_')
            if ref_col in df.columns:
                bw_ref_matrix[:, i] = pd.to_numeric(df[ref_col], errors='coerce').values

        # Aggregate metrics
        abs_bw = np.abs(bw_matrix)
        result_df['BW_Max_Abs_Delta'] = abs_bw.max(axis=1)
        result_df['BW_N_Strong_Signals'] = (abs_bw > 0.3).sum(axis=1)

        # Top-K indices
        bw_top_abs_idx = _vectorized_top_k_indices(bw_matrix, k)
        bw_top_gain_idx = _vectorized_top_k_signed(bw_matrix, k, 'gain')
        bw_top_loss_idx = _vectorized_top_k_signed(bw_matrix, k, 'loss')

        # Format strings
        result_df['BW_Top_Abs'] = _format_bw_signal_strings(
            bw_matrix, bw_ref_matrix, bw_top_abs_idx, bw_descriptions
        )
        result_df['BW_Top_Gains'] = _format_bw_signal_strings(
            bw_matrix, bw_ref_matrix, bw_top_gain_idx, bw_descriptions
        )
        result_df['BW_Top_Losses'] = _format_bw_signal_strings(
            bw_matrix, bw_ref_matrix, bw_top_loss_idx, bw_descriptions
        )
    else:
        result_df['BW_Max_Abs_Delta'] = 0
        result_df['BW_N_Strong_Signals'] = 0
        result_df['BW_Top_Abs'] = "None"
        result_df['BW_Top_Gains'] = "None"
        result_df['BW_Top_Losses'] = "None"

    # -----------------------------------------------------------------
    # MLM SIGNALS
    # -----------------------------------------------------------------
    mlm_sub, mlm_cols_used = _extract_mlm_columns(df, variant_type)

    # Individual MLM columns
    for col in mlm_cols_used:
        result_df[f'MLM_{col}'] = mlm_sub[col].values if col in mlm_sub.columns else np.nan

    # Summary string
    result_df['MLM_Summary'] = _format_mlm_summary_strings(df, variant_type)

    return result_df


# =============================================================================
# LLM PROMPT FORMATTER (SINGLE-ROW, ON-DEMAND)
# =============================================================================

def extract_signals_for_llm(
    row,
    delta_cols: List[str],
    metadata_df: Optional[pd.DataFrame] = None,
    k: int = 5,
    variant_type: str = 'snp'
) -> Dict:
    """
    Format signals for a SINGLE variant's LLM prompt.
    Call this at inference time, NOT during batch processing.

    Can work from either:
    - Pre-computed signal columns (if extract_signals_comprehensive_fast was run)
    - Raw delta columns (slower, self-contained fallback)
    """

    # --- Check if pre-computed columns exist ---
    if 'BED_Top_Abs' in row.index and 'BW_Top_Abs' in row.index:
        # Fast path: reformat pre-computed strings into LLM-friendly text
        bed_text = _reformat_compact_to_llm(row.get('BED_Top_Abs', 'None'), 'BED')
        bw_text = _reformat_compact_to_llm(row.get('BW_Top_Abs', 'None'), 'BW')
        mlm_text = _reformat_mlm_to_llm(row, variant_type)

        return {
            'BED_Signals_Text': bed_text,
            'BW_Signals_Text': bw_text,
            'MLM_Signals_Text': mlm_text,
            'Full_Signal_Block': f"""**BED Feature Changes (Top {k}):**
{bed_text}

**Epigenomic Track Changes (Top {k}):**
{bw_text}

**MLM Sequence Context:**
{mlm_text}"""
        }

    # --- Slow fallback: compute from raw columns ---
    # (Same logic as original, for standalone use)
    return _extract_signals_for_llm_from_raw(row, delta_cols, metadata_df, k, variant_type)


def _reformat_compact_to_llm(compact_str: str, signal_type: str) -> str:
    """Convert 'feature=delta(ref=X); ...' to verbose LLM format."""
    if compact_str == "None" or not compact_str:
        return "  None"

    lines = []
    for entry in compact_str.split("; "):
        # Parse: "feature=delta(ref=X)" or "feature=delta"
        if "=" not in entry:
            continue

        name_part, rest = entry.split("=", 1)

        # Check for ref value
        if "(ref=" in rest:
            delta_str, ref_part = rest.split("(ref=", 1)
            ref_str = ref_part.rstrip(")")
            lines.append(f"  - {name_part}: Δ={float(delta_str):+.3f} (REF={float(ref_str):.3f})")
        else:
            lines.append(f"  - {name_part}: Δ={float(rest):+.3f}")

    return "\n".join(lines) if lines else "  None"


def _reformat_mlm_to_llm(row, variant_type: str) -> str:
    """Reformat MLM columns into LLM-friendly text."""
    if variant_type == 'snp':
        cols = ['REF_5mer', 'ALT_5mer', 'LLR', 'MLM_Prior', 'MLM_Delta']
    else:
        cols = [
            'LLR', 'MLM_Prior', 'MLM_Delta', 'REF_5mer', 'ALT_5mer',
            'MLM_KL_mean', 'MLM_KL_max', 'MLM_logprob_ref', 'MLM_logprob_alt',
            'MLM_logprob_delta', 'EMB_cosine_dist', 'EMB_l2_dist',
            'EMB_max_pos_dist', 'EMB_mean_pos_dist'
        ]

    lines = []
    for col in cols:
        # Check both raw and prefixed versions
        val = row.get(f'MLM_{col}', row.get(col, np.nan))
        if not pd.isna(val):
            if isinstance(val, str):
                lines.append(f"  - {col}: {val}")
            else:
                lines.append(f"  - {col}: {val:.4f}")

    return "\n".join(lines) if lines else "  None"


def _extract_signals_for_llm_from_raw(row, delta_cols, metadata_df, k, variant_type):
    """Fallback: original row-wise extraction for single-variant use."""
    bed_delta_cols = [c for c in delta_cols if c.startswith('D_BED_')]
    bw_delta_cols = [c for c in delta_cols if c.startswith('D_BW_')]

    metadata_lookup = _build_metadata_lookup(metadata_df)

    # BED
    bed_entries = []
    for col in bed_delta_cols:
        delta = row.get(col, 0)
        if pd.isna(delta): delta = 0
        ref = row.get(col.replace('D_BED_', 'REF_BED_'), np.nan)
        feature = col.replace('D_BED_', '')
        bed_entries.append({'feature': feature, 'delta': delta,
                           'ref': ref if not pd.isna(ref) else None})

    bed_sorted = sorted(bed_entries, key=lambda x: abs(x['delta']), reverse=True)[:k]
    bed_lines = []
    for e in bed_sorted:
        line = f"  - {e['feature']}: Δ={e['delta']:+.3f}"
        if e['ref'] is not None: line += f" (REF={e['ref']:.3f})"
        bed_lines.append(line)
    bed_text = "\n".join(bed_lines) if bed_lines else "  None"

    # BW
    bw_entries = []
    for col in bw_delta_cols:
        delta = row.get(col, 0)
        if pd.isna(delta): delta = 0
        ref = row.get(col.replace('D_BW_', 'REF_BW_'), np.nan)
        track_id = col.replace('D_BW_', '')
        meta = metadata_lookup.get(track_id, {})
        if meta:
            desc = f"{meta.get('tissue', '?')} / {meta.get('biosample_type', '?')} / {meta.get('assay', '?')}"
            if meta.get('experiment_target'): desc += f" / {meta['experiment_target']}"
        else:
            desc = track_id
        bw_entries.append({'description': desc, 'delta': delta,
                          'ref': ref if not pd.isna(ref) else None})

    bw_sorted = sorted(bw_entries, key=lambda x: abs(x['delta']), reverse=True)[:k]
    bw_lines = []
    for e in bw_sorted:
        line = f"  - {e['description']}: Δ={e['delta']:+.3f}"
        if e['ref'] is not None: line += f" (REF={e['ref']:.3f})"
        bw_lines.append(line)
    bw_text = "\n".join(bw_lines) if bw_lines else "  None"

    # MLM
    mlm_text = _reformat_mlm_to_llm(row, variant_type)

    return {
        'BED_Signals_Text': bed_text,
        'BW_Signals_Text': bw_text,
        'MLM_Signals_Text': mlm_text,
        'Full_Signal_Block': f"""**BED Feature Changes (Top {k}):**
{bed_text}

**Epigenomic Track Changes (Top {k}):**
{bw_text}

**MLM Sequence Context:**
{mlm_text}"""
    }


# =============================================================================
