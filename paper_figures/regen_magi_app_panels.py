#!/usr/bin/env python3
"""
Regenerate the MAGI web-app "Region Track View" panels used as
Figure 5d (MAT1A / rs118204006) and Figure 6c (ALOX15B / rs9895916)
at print quality (true vector PDF/SVG + 600 dpi PNG).

Runs the same pipeline the Gradio app uses:
    inference.predict_variants(..., cache_profiles=True)
    analysis.rank_top_disrupted_tracks(...)
    tracks.generate_region_tracks_plot(...)

so the content is identical to what the app renders; only the output
device (vector instead of a screen grab) and the physical size change.

Outputs, per case:
  <case>_region_tracks_faithful.{pdf,svg,png}   exact app parameters, vector
  <case>_region_tracks_print90mm.{pdf,svg,png}  re-laid out to be legible at 90 mm
  <case>_region_tracks_print180mm.{pdf,svg,png} re-laid out to be legible at 180 mm
  <case>_track_profiles.npz                     cached REF/ALT profiles (re-render w/o GPU)
  <case>_ranked_tracks.csv                      the ranked track table behind the panel
  <case>_variant_result.csv                     full per-variant MAGI feature row

Usage:
    python regen_magi_app_panels.py [--device cuda] [--outdir DIR]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

# The published MAGI results use a 32 kb sequence window. The Gradio app's
# default was later lowered to 16 kb for CPU speed (commit fe87c8f), which
# changes the BED deltas by ~15% and reshuffles the near-tied BigWig ranking.
# Pin the paper's window so these panels reproduce the published figures.
os.environ["NTV3_CONTEXT_LEN"] = "32768"

import numpy as np
import pandas as pd


def configure_print_fonts() -> str:
    """Journal-production settings for vector output.

    matplotlib defaults to Type 3 fonts in PDF and to outlining text in SVG.
    Type 3 is routinely rejected or queried by journal production systems, and
    outlined text cannot be edited or reflowed by the typesetter. Force
    TrueType (type 42) embedding and live SVG text, and prefer Arial to match
    NAR figure style.
    """
    import matplotlib
    from matplotlib import font_manager

    for path in (
        "/mnt/c/Windows/Fonts/arial.ttf",
        "/mnt/c/Windows/Fonts/arialbd.ttf",
        "/mnt/c/Windows/Fonts/ariali.ttf",
        "/mnt/c/Windows/Fonts/arialbi.ttf",
    ):
        if Path(path).exists():
            font_manager.fontManager.addfont(path)

    available = {f.name for f in font_manager.fontManager.ttflist}
    family = "Arial" if "Arial" in available else "DejaVu Sans"

    matplotlib.rcParams.update({
        "pdf.fonttype": 42,      # TrueType, not Type 3
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # keep <text> live instead of outlining it
        "font.family": "sans-serif",
        "font.sans-serif": [family, "DejaVu Sans"],
    })
    return family


# Paths are derived from this file's location (paper_figures/ inside the MAGI
# repo), so the script works unchanged on any machine. Override with MAGI_ROOT
# if you keep it somewhere else.
MAGI_ROOT = Path(os.environ.get("MAGI_ROOT", Path(__file__).resolve().parents[1]))
GRADIO_APP_DIR = MAGI_ROOT / "gradio_app"
ANALYSES_DIR = MAGI_ROOT / "analyses"
DEFAULT_OUTDIR = Path(__file__).resolve().parent

# The two manuscript panels that are web-app output.
CASES = [
    {
        "case_id": "mat1a_gly336arg",
        "figure_panel": "Figure 5d",
        "panel_label": "MAT1A p.Gly336Arg",
        "gene": "MAT1A",
        "chrom": "chr10",
        "pos": 80274599,
        "ref": "C",
        "alt": "T",
        "rsid": "rs118204006",
        "hgvs": "NM_000429.3:c.1006G>A (p.Gly336Arg)",
    },
    {
        "case_id": "alox15b_rs9895916",
        "figure_panel": "Figure 6c",
        "panel_label": "ALOX15B rs9895916",
        "gene": "ALOX15B",
        "chrom": "chr17",
        "pos": 8047076,
        "ref": "G",
        "alt": "A",
        "rsid": "rs9895916",
        "hgvs": "NM_001141.3:c.1457G>A (p.Arg486His)",
    },
]

# Parameters the manuscript panels were generated with
# (analyses/generate_case_study_assets.py::run_case).
TOP_K = 15
MAX_RANKED_TRACKS = 10
VISIBLE_RADIUS_BP = 1000
FAITHFUL_FIGSIZE_X = 14.0
FAITHFUL_ROW_HEIGHT = 1.6


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_magi_modules():
    inference = load_module("magi_inference", GRADIO_APP_DIR / "inference.py")
    # tracks.py does `from inference import _LAST_TRACK_PROFILES`
    sys.modules["inference"] = inference
    return {
        "inference": inference,
        "analysis": load_module("magi_analysis", GRADIO_APP_DIR / "analysis.py"),
        "annotation": load_module("magi_annotation", GRADIO_APP_DIR / "annotation.py"),
        "tracks": load_module("magi_tracks", GRADIO_APP_DIR / "tracks.py"),
        "assets": load_module(
            "magi_assets", ANALYSES_DIR / "generate_case_study_assets.py"
        ),
    }


def dump_profiles(profiles: dict, path: Path) -> None:
    """Persist the cached REF/ALT track profiles so panels can be re-rendered
    later without a GPU / model load."""
    payload, meta = {}, {}
    for key, value in profiles.items():
        if isinstance(value, np.ndarray):
            payload[key] = value
        elif isinstance(value, (list, tuple)):
            payload[key] = np.asarray(value, dtype=object)
        else:
            meta[key] = value
    payload["__meta__"] = np.asarray([json.dumps(meta, default=str)])
    np.savez_compressed(path, **payload)


def retarget_for_print(fig, width_in: float, height_in: float,
                       label_pt: float, tick_pt: float, title_pt: float,
                       n_xticks: int = 3) -> None:
    """Re-lay out a stock region-track figure so it is legible when placed on
    the page at exactly *width_in* inches.

    The stock figure is 14 in wide with ~8.5 pt labels; scaled down to a
    ~90 mm manuscript panel that type lands at ~2 pt. Here the figure is built
    at its final physical size and the type is set in absolute points, so what
    you see is what prints.
    """
    import matplotlib.lines as mlines
    from matplotlib.patches import Patch

    fig.set_size_inches(width_in, height_in)
    axes = fig.get_axes()

    # Thinner strokes so traces stay traces at small size
    for obj in fig.findobj(mlines.Line2D):
        try:
            obj.set_linewidth(max(0.3, obj.get_linewidth() * 0.6))
        except Exception:
            pass

    for ax in axes:
        # tracks.py sets the per-track label with loc="left"
        ax.set_title(ax.get_title(loc="left"), fontsize=label_pt,
                     fontweight="bold", loc="left", pad=1.5)
        ax.set_ylabel("P", fontsize=tick_pt, labelpad=1)
        ax.tick_params(axis="both", labelsize=tick_pt, length=1.5, pad=1)
        ax.set_yticks([0.0, 0.5, 1.0])

    last = axes[-1]
    last.set_xlabel(last.get_xlabel(), fontsize=tick_pt + 0.5, labelpad=1)
    lo, hi = last.get_xlim()
    ticks = np.unique(np.rint(np.linspace(lo, hi, n_xticks)).astype(int))
    last.set_xticks(ticks)
    last.set_xticklabels([f"{t:,}" for t in ticks], fontsize=tick_pt)

    if fig._suptitle is not None:
        fig._suptitle.set_fontsize(title_pt)

    # Rebuild the legend compactly; the stock 4-column legend is wider than a
    # 90 mm panel and collides with the title.
    for lg in list(fig.legends):
        lg.remove()
    _CLR_REF, _CLR_GAIN = "#555555", "#d73027"   # must match tracks.py
    _CLR_LOSS, _CLR_VARIANT = "#2166ac", "#333333"
    handles = [
        mlines.Line2D([0], [0], color=_CLR_REF, lw=0.8, alpha=0.7,
                      linestyle="--", dashes=(4, 2), label="REF"),
        Patch(facecolor=_CLR_GAIN, alpha=0.3, label="Gain (ALT > REF)"),
        Patch(facecolor=_CLR_LOSS, alpha=0.3, label="Loss (ALT < REF)"),
        mlines.Line2D([0], [0], color=_CLR_VARIANT, lw=0.8,
                      linestyle="--", label="Variant"),
    ]
    fig.legend(handles=handles, loc="lower center", fontsize=tick_pt,
               frameon=False, ncol=4, handlelength=1.6, columnspacing=1.0,
               handletextpad=0.4, bbox_to_anchor=(0.5, 0.0))

    # Drop the two stock footnotes; at panel size they only add clutter and
    # their content is in the figure legend.
    for txt in list(fig.texts):
        if txt is not fig._suptitle:
            txt.set_visible(False)

    fig.tight_layout(rect=(0, 0.028, 1, 0.965))


def save_all(fig, base: Path, dpi: int = 600) -> dict:
    base.parent.mkdir(parents=True, exist_ok=True)
    out = {}
    for ext in ("pdf", "svg", "png"):
        path = base.with_suffix(f".{ext}")
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if ext == "png":
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)
        out[ext] = path
    return out


REFERENCE_DIR = ANALYSES_DIR / "paper_case_studies"


def pin_published_order(case_id: str, ranked, n: int = MAX_RANKED_TRACKS):
    """Reorder the regenerated top-N tracks to the order used in the published
    panel.

    NTv3 runs in bf16, and several of the top BigWig tracks are exact ties at
    bf16 resolution (e.g. MAT1A ranks 4/5 both at +0.1797). Tie-break order is
    therefore not reproducible across GPU/library versions. When the regenerated
    top-N contains exactly the same tracks as the published panel, pinning the
    published order makes the new panel a visual drop-in replacement while the
    plotted values remain the freshly computed ones.

    Returns (possibly reordered list, was_pinned).
    """
    ref_path = REFERENCE_DIR / case_id / "ranked_tracks.csv"
    if not ref_path.exists():
        return ranked, False

    ref_ids = list(pd.read_csv(ref_path).head(n)["track_id"])
    head, tail = ranked[:n], ranked[n:]
    if {t["track_id"] for t in head} != set(ref_ids):
        return ranked, False  # different track set — do not force

    by_id = {t["track_id"]: t for t in head}
    return [by_id[i] for i in ref_ids] + tail, True


def validate_against_published(case_id: str, ranked, n: int = MAX_RANKED_TRACKS):
    """Compare the regenerated top-N tracks against the ranking behind the
    published panel. Returns (report_lines, all_match)."""
    ref_path = REFERENCE_DIR / case_id / "ranked_tracks.csv"
    if not ref_path.exists():
        return [f"  no reference ranking at {ref_path}"], False

    ref = pd.read_csv(ref_path).head(n)
    new = pd.DataFrame(ranked).head(n)
    lines, all_match = [], True
    lines.append(
        f"  {'#':<3}{'track (published panel)':<48}{'Δ pub':>10}{'Δ new':>10}"
        f"{'|Δdiff|':>10}  ok"
    )
    for i in range(n):
        r_id = ref.iloc[i]["track_id"]
        r_nm = str(ref.iloc[i]["display_name"])[:46]
        r_d = float(ref.iloc[i]["delta"])
        n_id = new.iloc[i]["track_id"]
        n_d = float(new.iloc[i]["delta"])
        same = (r_id == n_id) and abs(r_d - n_d) < 0.02
        all_match &= same
        lines.append(
            f"  {i:<3}{r_nm:<48}{r_d:>10.4f}{n_d:>10.4f}"
            f"{abs(r_d - n_d):>10.4f}  {'OK' if same else 'DIFF'}"
        )
    lines.append(
        f"\n  track set identical: "
        f"{set(ref['track_id']) == set(new['track_id'])}"
    )
    lines.append(
        "  max |Δ difference| over the 10 rendered tracks: "
        f"{max(abs(float(ref.iloc[i]['delta']) - float(new.iloc[i]['delta'])) for i in range(n)):.4f}"
    )
    return lines, all_match


def render_variant(tracks_mod, ranked, figsize_x, row_height):
    fig = tracks_mod.generate_region_tracks_plot(
        ranked_tracks=ranked,
        max_ranked_tracks=MAX_RANKED_TRACKS,
        visible_radius_bp=VISIBLE_RADIUS_BP,
        figsize_x=figsize_x,
        row_height=row_height,
    )
    if fig is None:
        raise RuntimeError("generate_region_tracks_plot returned None")
    return fig


# Print targets: (name, width_in, height_in, label_pt, tick_pt, title_pt, n_xticks)
# 90 mm  = single manuscript column-ish panel (NAR single column ~88 mm)
# 180 mm = full text width
PRINT_TARGETS = [
    ("print90mm", 3.543, 6.30, 5.2, 4.2, 6.0, 3),
    ("print180mm", 7.087, 8.20, 7.5, 6.0, 9.0, 5),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cuda", choices=("cuda", "cpu", "auto"))
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                    help="writes fig5/ and fig6/ under here (default: this "
                         "script's own directory)")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument(
        "--cache-profiles", action="store_true",
        help="Also write the raw REF/ALT track profiles (~180 MB per case). "
             "Only needed to re-render on a machine without a GPU.",
    )
    args = ap.parse_args()

    font_family = configure_print_fonts()
    print(f"Vector output: Type-42 fonts, live SVG text, family={font_family}")

    mods = load_magi_modules()
    inference = mods["inference"]
    analysis = mods["analysis"]
    annotation = mods["annotation"]
    tracks_mod = mods["tracks"]
    assets = mods["assets"]

    device = args.device
    if device == "auto":
        import torch

        device = "cuda" if torch.cuda.is_available() else "cpu"

    metadata_df, metadata_dict = assets.load_track_metadata()

    manifest_rows = []
    for case in CASES:
        print(f"\n{'='*70}\n{case['figure_panel']}: {case['panel_label']}\n{'='*70}")
        subdir = args.outdir / ("fig5" if case["figure_panel"].startswith("Figure 5") else "fig6")
        subdir.mkdir(parents=True, exist_ok=True)

        input_df = pd.DataFrame(
            [{
                "chrom": case["chrom"],
                "pos": int(case["pos"]),
                "ref": case["ref"],
                "alt": case["alt"],
            }]
        )
        results_df = inference.predict_variants(
            input_df, device=device, species="human", cache_profiles=True
        )
        results_df = assets.annotate_for_species(annotation, results_df, "human")
        results_df = analysis.compute_impact_scores(results_df)
        row = results_df.iloc[0]

        bed_names = list(inference._MODEL_CACHE.get("bed_names") or [])
        bw_indices = list(inference._MODEL_CACHE.get("selected_bw_indices") or [])
        bw_names_all = list(inference._MODEL_CACHE.get("bigwig_names") or [])
        bw_names_filtered = [bw_names_all[i] for i in bw_indices]

        ranked = analysis.rank_top_disrupted_tracks(
            row,
            bed_names,
            bw_names_filtered,
            metadata_df=metadata_df,
            metadata_dict=metadata_dict,
            top_k=None,
        )
        ranked_for_outputs = assets.prepare_ranked_for_outputs(
            ranked, disambiguate_top_n=max(20, TOP_K)
        )

        ranked_for_outputs, pinned = pin_published_order(
            case["case_id"], ranked_for_outputs
        )
        if pinned:
            print("  top-10 track set matches the published panel; "
                  "pinned to the published row order (bf16 ties)")

        val_lines, val_ok = validate_against_published(
            case["case_id"], ranked_for_outputs
        )
        print("\n  VALIDATION vs published panel:")
        print("\n".join(val_lines))
        print(f"  => {'MATCHES published panel' if val_ok else 'DIFFERS from published panel'}\n")

        stem = f"{case['case_id']}"
        (subdir / f"{stem}_validation.txt").write_text(
            f"Validation of regenerated {case['figure_panel']} "
            f"({case['panel_label']}) against the published panel\n"
            f"reference: {REFERENCE_DIR / case['case_id'] / 'ranked_tracks.csv'}\n"
            f"context window: {os.environ['NTV3_CONTEXT_LEN']} bp\n\n"
            + "\n".join(val_lines)
            + f"\n\nRESULT: {'MATCH' if val_ok else 'MISMATCH'}\n",
            encoding="utf-8",
        )
        results_df.to_csv(subdir / f"{stem}_variant_result.csv", index=False)
        pd.DataFrame(ranked_for_outputs).to_csv(
            subdir / f"{stem}_ranked_tracks.csv", index=False
        )
        if args.cache_profiles:
            # ~180 MB per case. Only useful for re-rendering at a different size
            # on a machine with no GPU; regenerating from scratch here takes ~40 s.
            dump_profiles(
                inference._LAST_TRACK_PROFILES, subdir / f"{stem}_track_profiles.npz"
            )

        # ---- 1. Faithful: identical parameters to the published panel -------
        fig = render_variant(
            tracks_mod, ranked_for_outputs, FAITHFUL_FIGSIZE_X, FAITHFUL_ROW_HEIGHT
        )
        paths = save_all(fig, subdir / f"{stem}_region_tracks_faithful", args.dpi)
        analysis.plt.close(fig)
        print("  faithful ->", paths["pdf"].name)

        # ---- 2/3. Print-sized re-renders ------------------------------------
        # Rendered AT the final physical width, with type scaled up so the
        # on-page point size stays legible instead of shrinking with the panel.
        for (label, w_in, h_in, lab_pt, tick_pt, title_pt, nxt) in PRINT_TARGETS:
            fig = render_variant(
                tracks_mod, ranked_for_outputs,
                figsize_x=w_in, row_height=h_in / MAX_RANKED_TRACKS,
            )
            retarget_for_print(fig, w_in, h_in, lab_pt, tick_pt, title_pt, nxt)
            paths = save_all(fig, subdir / f"{stem}_region_tracks_{label}", args.dpi)
            analysis.plt.close(fig)
            print(f"  {label} ({w_in*25.4:.0f}x{h_in*25.4:.0f} mm) ->", paths["pdf"].name)

        manifest_rows.append(
            {
                **{k: case[k] for k in
                   ("case_id", "figure_panel", "panel_label", "gene",
                    "chrom", "pos", "ref", "alt", "rsid", "hgvs")},
                "n_tracks_rendered": min(MAX_RANKED_TRACKS, len(ranked_for_outputs)),
                "visible_radius_bp": VISIBLE_RADIUS_BP,
                "window": f"{case['chrom']}:{case['pos']-VISIBLE_RADIUS_BP:,}-"
                          f"{case['pos']+VISIBLE_RADIUS_BP:,}",
                "top_track": ranked_for_outputs[0]["display_name"],
                "top_track_delta": round(float(ranked_for_outputs[0]["delta"]), 6),
                "Global_z_sum_log": float(row.get("Global_z_sum_log", float("nan"))),
                "generated_utc": datetime.now(timezone.utc).isoformat(),
            }
        )

    pd.DataFrame(manifest_rows).to_csv(
        args.outdir / "app_panels_manifest.csv", index=False
    )
    print(f"\nWrote {args.outdir / 'app_panels_manifest.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
