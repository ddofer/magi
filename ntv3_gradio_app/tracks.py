#!/usr/bin/env python3
"""
Track Visualization Module
===========================
Generates region-level fill-between plots from cached NTv3 track profiles.

Shows the continuous predicted probability of each genomic feature across the
analysis window, highlighting the variant position and the most impacted tracks.
Unified track selection is driven by the ranked list from analysis.py to
ensure consistency with the fingerprint bar chart.

Usage:
    from tracks import generate_region_tracks_plot
    fig = generate_region_tracks_plot(ranked_tracks=ranked, visible_radius_bp=1000)
"""

from typing import Optional, Dict, List

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

matplotlib.use("Agg")

# LOF/GOF palette — consistent with fingerprint bar chart
_CLR_REF = "#555555"  # neutral gray for reference
_CLR_GAIN = "#d73027"  # red — gain of function (ALT > REF)
_CLR_LOSS = "#2166ac"  # blue — loss of function (ALT < REF)
_CLR_VARIANT = "#333333"  # dark gray for variant position line
_BED_BG = "#f7f7f7"  # subtle background for BED group
_BW_BG = "#fffef5"  # subtle warm tint for BigWig group


def get_track_view_bounds() -> Dict[str, Optional[int]]:
    """Return the exact symmetric track-view bounds available for the current variant."""
    from inference import _LAST_TRACK_PROFILES

    profiles = _LAST_TRACK_PROFILES
    if not profiles:
        return {"max_radius": None, "window_start": None, "window_end": None}

    pos = int(profiles["pos"])
    vcenter = int(profiles["variant_center"])
    candidate_radii: List[int] = []

    for prefix in ("bed", "bw"):
        if profiles.get(f"{prefix}_ref") is None:
            continue
        track_len = profiles.get(f"{prefix}_track_len")
        track_start = profiles.get(f"{prefix}_track_start")
        if track_len is None or track_start is None:
            continue

        left_radius = max(0, int(vcenter - int(track_start)))
        right_radius = max(0, int(int(track_start) + int(track_len) - 1 - vcenter))
        candidate_radii.append(min(left_radius, right_radius))

    if not candidate_radii:
        return {"max_radius": None, "window_start": None, "window_end": None}

    max_radius = max(8, int(min(candidate_radii)))
    return {
        "max_radius": max_radius,
        "window_start": pos - max_radius,
        "window_end": pos + max_radius,
    }


def generate_region_tracks_plot(
    ranked_tracks: Optional[List[Dict]] = None,
    metadata_df: Optional[pd.DataFrame] = None,
    metadata_dict: Optional[Dict[str, Dict[str, str]]] = None,
    top_k_bed: int = 3,
    top_k_bw: int = 3,
    visible_radius_bp: int = 1000,
    max_ranked_tracks: int = 10,
    figsize_x: float = 14.0,
    row_height: float = 1.6,
) -> Optional[plt.Figure]:
    """
    Generate fill-between region view from cached track profiles.

    If *ranked_tracks* is provided (from analysis.rank_top_disrupted_tracks),
    those exact tracks are shown (ensuring consistency with fingerprint).
    Otherwise falls back to auto-selecting top-N by point delta.

    Args:
        ranked_tracks: Pre-ranked list of dicts with keys {track_id, track_type, display_name, delta}.
        metadata_df / metadata_dict: BigWig metadata (used only for fallback auto-select).
        top_k_bed / top_k_bw: Fallback auto-select counts (ignored when ranked_tracks given).
        visible_radius_bp: Half-width of visible window centred on variant (bp).
        max_ranked_tracks: Maximum number of ranked tracks to render.
        figsize_x: Figure width in inches.
        row_height: Height per subplot row.
    """
    from inference import _LAST_TRACK_PROFILES

    profiles = _LAST_TRACK_PROFILES
    if not profiles:
        return None

    chrom = profiles["chrom"]
    pos = profiles["pos"]
    ref = profiles["ref"]
    alt = profiles["alt"]
    vcenter = profiles["variant_center"]
    bed_names = profiles["bed_names"]
    bigwig_names = profiles["bigwig_names"]
    selected_bw = profiles["selected_bw_indices"]

    # ── Build per-track rendering list ──────────────────────────────────
    all_tracks: List[dict] = []

    if ranked_tracks:
        # Use unified ranking — pull continuous arrays from cache
        for item in ranked_tracks[: max(1, int(max_ranked_tracks))]:
            tid = item["track_id"]
            ttype = item["track_type"]

            if ttype == "BED":
                if profiles.get("bed_ref") is None:
                    continue
                bed_ref = profiles["bed_ref"]
                bed_alt = profiles["bed_alt"]
                track_len = profiles["bed_track_len"]
                track_start = profiles["bed_track_start"]
                try:
                    idx = bed_names.index(tid)
                except ValueError:
                    continue
                bed_pos = vcenter - track_start
                delta_at = (
                    float(bed_alt[bed_pos, idx] - bed_ref[bed_pos, idx])
                    if 0 <= bed_pos < track_len
                    else 0.0
                )
                all_tracks.append(
                    {
                        "name": item["display_name"],
                        "ref": bed_ref[:, idx],
                        "alt": bed_alt[:, idx],
                        "type": "BED",
                        "track_start": track_start,
                        "delta_at_variant": delta_at,
                    }
                )

            elif ttype == "BigWig":
                if profiles.get("bw_ref") is None:
                    continue
                bw_ref = profiles["bw_ref"]
                bw_alt = profiles["bw_alt"]
                track_len = profiles["bw_track_len"]
                track_start = profiles["bw_track_start"]
                try:
                    global_idx = bigwig_names.index(tid)
                except ValueError:
                    continue
                bw_pos = vcenter - track_start
                delta_at = (
                    float(bw_alt[bw_pos, global_idx] - bw_ref[bw_pos, global_idx])
                    if 0 <= bw_pos < track_len
                    else 0.0
                )
                all_tracks.append(
                    {
                        "name": item["display_name"],
                        "ref": bw_ref[:, global_idx],
                        "alt": bw_alt[:, global_idx],
                        "type": "BigWig",
                        "track_start": track_start,
                        "delta_at_variant": delta_at,
                    }
                )
    else:
        # Fallback: auto-select by point delta (legacy behaviour)
        all_tracks = _auto_select_tracks(
            profiles, metadata_df, metadata_dict, top_k_bed, top_k_bw
        )

    if not all_tracks:
        return None

    # ── Separate BED and BigWig groups for visual banding ───────────────
    bed_group = [t for t in all_tracks if t["type"] == "BED"]
    bw_group = [t for t in all_tracks if t["type"] == "BigWig"]
    ordered = bed_group + bw_group
    n_tracks = len(ordered)

    fig, axes = plt.subplots(
        n_tracks,
        1,
        figsize=(figsize_x, row_height * n_tracks + 1.0),
        sharex=True,
        squeeze=False,
    )
    axes = axes.flatten()

    genomic_origin = pos - vcenter  # genomic coord at token 0
    n_bed = len(bed_group)
    requested_radius_bp = max(int(visible_radius_bp), 8)
    view_bounds = get_track_view_bounds()
    max_radius_bp = view_bounds.get("max_radius")
    effective_radius_bp = (
        min(requested_radius_bp, max_radius_bp)
        if max_radius_bp is not None
        else requested_radius_bp
    )
    xleft = pos - effective_radius_bp
    xright = pos + effective_radius_bp

    for ax_idx, track_info in enumerate(ordered):
        ax = axes[ax_idx]
        ts = track_info["track_start"]
        arr_len = len(track_info["ref"])
        x_genomic = np.arange(ts, ts + arr_len) + genomic_origin

        window_mask = (x_genomic >= xleft) & (x_genomic <= xright)
        if not np.any(window_mask):
            nearest_idx = int(np.argmin(np.abs(x_genomic - pos)))
            lo = max(0, nearest_idx - 1)
            hi = min(arr_len, nearest_idx + 2)
            window_mask = np.zeros(arr_len, dtype=bool)
            window_mask[lo:hi] = True

        x_window = x_genomic[window_mask]
        ref_y = track_info["ref"][window_mask]
        alt_y = track_info["alt"][window_mask]

        # Subtle group background
        bg = _BED_BG if ax_idx < n_bed else _BW_BG
        ax.set_facecolor(bg)

        # ALT coloured by delta sign at variant (draw first, behind REF)
        delta = track_info["delta_at_variant"]
        alt_clr = _CLR_GAIN if delta > 0 else _CLR_LOSS
        ax.plot(x_window, alt_y, color=alt_clr, linewidth=1.0, alpha=0.85, label="ALT")

        # Directional delta fill: red where gain, blue where loss
        ax.fill_between(
            x_window,
            ref_y,
            alt_y,
            where=(alt_y > ref_y),
            color=_CLR_GAIN,
            alpha=0.22,
            interpolate=True,
        )
        ax.fill_between(
            x_window,
            ref_y,
            alt_y,
            where=(alt_y <= ref_y),
            color=_CLR_LOSS,
            alpha=0.22,
            interpolate=True,
        )

        # REF as gray dashed line — drawn on top so it stays visible
        ax.plot(x_window, ref_y, color=_CLR_REF, linewidth=1.1, alpha=0.7,
                linestyle="--", dashes=(4, 2), label="REF")

        # Variant position marker
        ax.axvline(pos, color=_CLR_VARIANT, linewidth=1.2, linestyle="--", alpha=0.7)

        # Track label
        direction = "↑ Gain" if delta > 0 else "↓ Loss"
        label_txt = f"{track_info['name']}  (Δ = {delta:+.4f}  {direction})"
        ax.set_title(label_txt, fontsize=8.5, fontweight="bold", loc="left", pad=3)
        ax.set_ylabel("P", fontsize=7, labelpad=1)
        ax.tick_params(axis="both", labelsize=6.5)
        ax.set_ylim(bottom=0)

        # Minimal Tufte-style axes — only left spine + bottom on last
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax_idx < n_tracks - 1:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="x", length=0)

    # Group separator line between BED and BigWig
    if n_bed > 0 and len(bw_group) > 0:
        # Use the axis position of the first BigWig row to draw a thin separator
        sep_ax = axes[n_bed]
        sep_ax.annotate(
            "",
            xy=(0, 1),
            xycoords="axes fraction",
            xytext=(1, 1),
            textcoords="axes fraction",
            arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.8),
        )

    # X-axis
    axes[-1].set_xlabel(f"Genomic position ({chrom})", fontsize=9)
    axes[-1].ticklabel_format(axis="x", style="plain", useOffset=False)
    axes[-1].set_xlim(xleft, xright)
    tick_positions = np.linspace(xleft, xright, num=5)
    tick_positions = np.unique(np.rint(tick_positions).astype(int))
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels([f"{tick:,}" for tick in tick_positions], fontsize=6.5)

    # Suptitle
    fig.suptitle(
        f"Region Track View — {chrom}:{pos:,} {ref}>{alt}",
        fontsize=11,
        fontweight="bold",
        y=1.0,
    )

    # Legend
    legend_elements = [
        plt.Line2D([0], [0], color=_CLR_REF, linewidth=1.2, alpha=0.7,
                   linestyle="--", dashes=(4, 2), label="REF"),
        Patch(facecolor=_CLR_GAIN, alpha=0.3, label="Gain (ALT > REF)"),
        Patch(facecolor=_CLR_LOSS, alpha=0.3, label="Loss (ALT < REF)"),
        plt.Line2D(
            [0], [0], color=_CLR_VARIANT, linewidth=1.2, linestyle="--", label="Variant"
        ),
    ]
    fig.legend(
        handles=legend_elements,
        loc="upper right",
        fontsize=7,
        frameon=True,
        ncol=4,
        bbox_to_anchor=(0.99, 0.99),
    )

    window_note = (
        f"Visible window: {chrom}:{xleft:,}-{xright:,} "
        f"(radius {effective_radius_bp:,} bp)"
    )
    if effective_radius_bp != requested_radius_bp:
        window_note += (
            f" | requested {requested_radius_bp:,} bp, limited by available track signal"
        )
    fig.text(0.5, 0.016, window_note, ha="center", fontsize=7.5, color="#555555")

    # Footnote: P = predicted probability
    fig.text(
        0.01, 0.002,
        "P = predicted probability of genomic feature",
        fontsize=7, color="#777777", style="italic",
    )

    plt.tight_layout(rect=(0, 0.03, 1, 0.96))
    return fig


# ── Fallback auto-select (preserves legacy behaviour) ──────────────────
def _auto_select_tracks(profiles, metadata_df, metadata_dict, top_k_bed, top_k_bw):
    """Select top tracks by point-delta when no pre-ranked list is given."""
    tracks = []
    vcenter = profiles["variant_center"]
    bed_names = profiles["bed_names"]
    bigwig_names = profiles["bigwig_names"]
    selected_bw = profiles["selected_bw_indices"]

    if profiles.get("bed_ref") is not None:
        bed_ref = profiles["bed_ref"]
        bed_alt = profiles["bed_alt"]
        track_len = profiles["bed_track_len"]
        track_start = profiles["bed_track_start"]
        bed_pos = vcenter - track_start
        if 0 <= bed_pos < track_len:
            deltas = np.abs(bed_alt[bed_pos] - bed_ref[bed_pos])
            top_idx = np.argsort(deltas)[-top_k_bed:][::-1]
            for idx in top_idx:
                name = bed_names[idx] if idx < len(bed_names) else f"BED_{idx}"
                tracks.append(
                    {
                        "name": name,
                        "ref": bed_ref[:, idx],
                        "alt": bed_alt[:, idx],
                        "type": "BED",
                        "track_start": track_start,
                        "delta_at_variant": float(
                            bed_alt[bed_pos, idx] - bed_ref[bed_pos, idx]
                        ),
                    }
                )

    if profiles.get("bw_ref") is not None:
        bw_ref = profiles["bw_ref"]
        bw_alt = profiles["bw_alt"]
        track_len = profiles["bw_track_len"]
        track_start = profiles["bw_track_start"]
        bw_pos = vcenter - track_start
        if 0 <= bw_pos < track_len and selected_bw:
            sel_ref = bw_ref[bw_pos, selected_bw]
            sel_alt = bw_alt[bw_pos, selected_bw]
            abs_d = np.abs(sel_alt - sel_ref)
            top_local = np.argsort(abs_d)[-top_k_bw:][::-1]
            for li in top_local:
                gi = selected_bw[li]
                tid = bigwig_names[gi] if gi < len(bigwig_names) else f"BW_{gi}"
                display = _resolve_track_name(tid, metadata_df, metadata_dict)
                tracks.append(
                    {
                        "name": display,
                        "ref": bw_ref[:, gi],
                        "alt": bw_alt[:, gi],
                        "type": "BigWig",
                        "track_start": track_start,
                        "delta_at_variant": float(
                            bw_alt[bw_pos, gi] - bw_ref[bw_pos, gi]
                        ),
                    }
                )
    return tracks


def _resolve_track_name(
    track_id: str,
    metadata_df: Optional[pd.DataFrame] = None,
    metadata_dict: Optional[Dict[str, Dict[str, str]]] = None,
) -> str:
    """Resolve a BigWig track ID to a human-readable name."""
    if metadata_dict:
        meta = metadata_dict.get(track_id)
        if meta:
            parts = [
                p
                for p in [
                    meta.get("tissue", ""),
                    meta.get("assay", ""),
                    meta.get("target", ""),
                ]
                if p.strip() and p != "nan"
            ]
            if parts:
                name = " | ".join(parts)
                return name[:55] if len(name) > 55 else name

    if metadata_df is not None:
        rows = metadata_df[metadata_df["file_id"] == track_id]
        if not rows.empty:
            r = rows.iloc[0]
            parts = [
                str(p)
                for p in [
                    r.get("tissue", ""),
                    r.get("assay", ""),
                    r.get("experiment_target", ""),
                ]
                if pd.notna(p) and str(p).strip()
            ]
            if parts:
                name = " | ".join(parts)
                return name[:55] if len(name) > 55 else name

    return track_id[:40]
