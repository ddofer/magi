#!/usr/bin/env python3
"""Generate MAGI case-study assets for manuscript figures."""

from __future__ import annotations

import argparse
from collections import Counter
import importlib.util
import json
import os
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
GRADIO_APP_DIR = REPO_ROOT / "gradio_app"
DEFAULT_MANIFEST = Path(__file__).with_name("case_study_manifest.csv")
DEFAULT_OUTPUT_DIR = Path(__file__).with_name("paper_case_studies")
TRACK_METADATA_FILE = GRADIO_APP_DIR / "data" / "functional_tracks_metadata_human.csv"

DISPLAY_NAME_OVERRIDES = {
    "always_on_exon": "Always-on exon",
    "protein_coding_gene": "Protein-coding gene",
}

MANIFEST_COLUMNS = [
    "case_id",
    "panel_label",
    "species",
    "gene",
    "chrom",
    "pos",
    "ref",
    "alt",
    "transcript",
    "transcript_hgvs",
    "protein_hgvs",
    "source_type",
    "source_id",
    "status",
    "provenance_source",
    "provenance_note",
    "manuscript_note",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate MAGI fingerprint panels and track tables from a case manifest."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="CSV manifest describing manuscript case-study variants.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where per-case outputs and review files will be written.",
    )
    parser.add_argument(
        "--case-id",
        action="append",
        dest="case_ids",
        default=None,
        help="Optional case_id filter. Repeat to run more than one case.",
    )
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "cuda"),
        default="auto",
        help="Inference device. Defaults to auto.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=15,
        help="Number of ranked tracks to show in each fingerprint panel.",
    )
    parser.add_argument(
        "--figure-dpi",
        type=int,
        default=600,
        help="DPI for raster figure export.",
    )
    return parser.parse_args()


def load_module(module_name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Failed to create import spec for {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def load_magi_modules() -> Dict[str, Any]:
    os.environ.setdefault("MPLBACKEND", "Agg")
    inference_module = load_module(
        "magi_gradio_inference", GRADIO_APP_DIR / "inference.py"
    )
    sys.modules.setdefault("inference", inference_module)

    return {
        "inference": inference_module,
        "analysis": load_module("magi_gradio_analysis", GRADIO_APP_DIR / "analysis.py"),
        "annotation": load_module(
            "magi_gradio_annotation", GRADIO_APP_DIR / "annotation.py"
        ),
        "interpretation": load_module(
            "magi_gradio_interpretation", GRADIO_APP_DIR / "interpretation.py"
        ),
        "tracks": load_module("magi_gradio_tracks", GRADIO_APP_DIR / "tracks.py"),
    }


def load_manifest(manifest_path: Path) -> pd.DataFrame:
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    manifest = pd.read_csv(manifest_path, dtype=str).fillna("")
    for column in MANIFEST_COLUMNS:
        if column not in manifest.columns:
            manifest[column] = ""

    manifest = manifest[MANIFEST_COLUMNS].copy()
    manifest["pos"] = pd.to_numeric(manifest["pos"], errors="coerce")

    string_columns = [column for column in MANIFEST_COLUMNS if column != "pos"]
    for column in string_columns:
        manifest[column] = manifest[column].astype(str).str.strip()

    manifest["ref"] = manifest["ref"].str.upper()
    manifest["alt"] = manifest["alt"].str.upper()
    return manifest


def normalize_chrom(chrom: str) -> str:
    cleaned = str(chrom).strip()
    if not cleaned:
        return cleaned
    return cleaned if cleaned.startswith("chr") else f"chr{cleaned}"


def row_has_coordinates(case: pd.Series) -> bool:
    return bool(
        str(case.get("chrom", "")).strip()
        and pd.notna(case.get("pos"))
        and str(case.get("ref", "")).strip()
        and str(case.get("alt", "")).strip()
    )


def is_resolved_case(case: pd.Series) -> bool:
    return str(
        case.get("status", "")
    ).strip().lower() == "resolved" and row_has_coordinates(case)


def to_builtin(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): to_builtin(val) for key, val in value.items()}
    if isinstance(value, set):
        return [to_builtin(item) for item in sorted(value)]
    if isinstance(value, (list, tuple)):
        return [to_builtin(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    return value


def stringify_cell(value: Any) -> Any:
    if isinstance(value, set):
        return ";".join(str(item) for item in sorted(value))
    if isinstance(value, (list, tuple)):
        return ";".join(str(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(to_builtin(value), sort_keys=True)
    return value


def make_csv_safe(df: pd.DataFrame) -> pd.DataFrame:
    safe = df.copy()
    for column in safe.columns:
        safe[column] = safe[column].map(stringify_cell)
    return safe


def relpath(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def choose_device(inference_module, requested_device: str) -> str:
    if requested_device == "cpu":
        return "cpu"
    if requested_device == "cuda":
        return "cuda"
    return "cuda" if inference_module.torch.cuda.is_available() else "cpu"


def load_track_metadata() -> tuple[pd.DataFrame, Dict[str, Dict[str, str]]]:
    metadata_df = pd.read_csv(TRACK_METADATA_FILE)
    metadata_dict = {
        row["file_id"]: {
            "tissue": str(row.get("tissue", "") or ""),
            "assay": str(row.get("assay", "") or ""),
            "target": str(row.get("experiment_target", "") or ""),
        }
        for _, row in metadata_df.iterrows()
    }
    return metadata_df, metadata_dict


def top_signal(items: Iterable[Dict[str, Any]]) -> tuple[str, Optional[float]]:
    items = list(items)
    if not items:
        return "", None
    first = items[0]
    label = str(
        first.get("label") or first.get("display_name") or first.get("track_id") or ""
    )
    delta = first.get("delta")
    if delta is None:
        delta = first.get("value")
    try:
        delta_value = float(delta)
    except (TypeError, ValueError):
        delta_value = None
    return label, delta_value


def prettify_display_name(label: str) -> str:
    text = str(label or "").strip()
    if not text:
        return text
    if text in DISPLAY_NAME_OVERRIDES:
        return DISPLAY_NAME_OVERRIDES[text]
    if text.upper() == text and "_" not in text:
        return text
    return text.replace("_", " ")


def prepare_ranked_for_outputs(
    ranked: List[Dict[str, Any]],
    disambiguate_top_n: Optional[int] = None,
) -> List[Dict[str, Any]]:
    base_names = [
        prettify_display_name(item.get("display_name") or item.get("track_id") or "")
        for item in ranked
    ]
    preview_count = (
        len(base_names) if disambiguate_top_n is None else max(0, disambiguate_top_n)
    )
    counts = Counter(base_names[:preview_count])
    seen: Dict[str, int] = {}
    prepared: List[Dict[str, Any]] = []

    for index, (item, base_name) in enumerate(zip(ranked, base_names)):
        display_name = base_name
        if index < preview_count and counts[base_name] > 1:
            seen[base_name] = seen.get(base_name, 0) + 1
            display_name = f"{base_name} [{seen[base_name]}]"

        updated_item = dict(item)
        updated_item["display_name"] = display_name
        prepared.append(updated_item)

    return prepared


def polish_fingerprint_figure(fig, panel_label: str) -> None:
    if not getattr(fig, "axes", None):
        return

    ax = fig.axes[0]
    ax.set_title(
        f"{panel_label}\nTop disrupted MAGI signals",
        fontsize=14,
        fontweight="bold",
        pad=14,
    )
    ax.set_xlabel("Delta probability (Alt - Ref)", fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", labelsize=9)

    from matplotlib.patches import Patch

    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

    legend_elements = [
        Patch(facecolor="#d73027", alpha=0.8, label="Gain of function"),
        Patch(facecolor="#4575b4", alpha=0.8, label="Loss of function"),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=False,
        fontsize=9,
        title="Direction",
        title_fontsize=9,
    )
    fig.tight_layout()


def save_figure(fig, output_base: Path, dpi: int) -> Dict[str, Path]:
    output_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_base.with_suffix(".png")
    pdf_path = output_base.with_suffix(".pdf")
    svg_path = output_base.with_suffix(".svg")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    return {"png": png_path, "pdf": pdf_path, "svg": svg_path}


def annotate_for_species(
    annotation_module, results_df: pd.DataFrame, species: str
) -> pd.DataFrame:
    if str(species).strip() == "human":
        return annotation_module.annotate_dataframe(results_df)

    annotated = results_df.copy()
    annotation_columns = list(getattr(annotation_module, "ANNOTATION_COLUMNS", []))
    for column in annotation_columns:
        annotated[column] = 0
    annotated["transcript_set"] = [set() for _ in range(len(annotated))]
    annotated["promoter_transcript_set"] = [set() for _ in range(len(annotated))]
    annotated["region"] = "non_human_unavailable"
    annotated["region_class"] = "NON_HUMAN"
    annotated["gene_name"] = "N/A (non-human)"
    return annotated


def build_case_provenance(case: pd.Series) -> Dict[str, Any]:
    return {
        "case_id": case["case_id"],
        "panel_label": case["panel_label"],
        "species": case["species"],
        "gene": case["gene"],
        "chrom": case["chrom"],
        "pos": to_builtin(case.get("pos")),
        "ref": case["ref"],
        "alt": case["alt"],
        "transcript": case["transcript"],
        "transcript_hgvs": case["transcript_hgvs"],
        "protein_hgvs": case["protein_hgvs"],
        "source_type": case["source_type"],
        "source_id": case["source_id"],
        "status": case["status"],
        "provenance_source": case["provenance_source"],
        "provenance_note": case["provenance_note"],
        "manuscript_note": case["manuscript_note"],
        "assembly": "GRCh38"
        if str(case.get("species", "")).strip() == "human"
        else None,
        "recorded_at": datetime.now(timezone.utc).isoformat(),
    }


def build_review_row(
    case: pd.Series,
    summary_row: Dict[str, Any],
    top_k: int,
    figure_dpi: int,
    fingerprint_paths: Dict[str, Path],
    region_tracks_paths: Optional[Dict[str, Path]] = None,
) -> Dict[str, Any]:
    top_abs_delta = summary_row.get("top_abs_delta")
    llr = summary_row.get("LLR")
    relevance_review = "pass"
    review_note = "Signals and exports are suitable for manuscript assembly."
    if top_abs_delta is None or (
        top_abs_delta < 0.03 and (llr is None or abs(llr) < 0.5)
    ):
        relevance_review = "manual_review"
        review_note = "Signals are weak enough that narrative alignment should be checked manually."

    return {
        "case_id": case["case_id"],
        "panel_label": case["panel_label"],
        "status": summary_row.get("status", "generated"),
        "resolution_review": "pass" if figure_dpi >= 300 else "manual_review",
        "subfigure_review": "pass"
        if summary_row.get("fingerprint_label_count", 0) <= top_k
        else "manual_review",
        "relevance_review": relevance_review,
        "figure_dpi": figure_dpi,
        "figure_size_inches": "10x6",
        "label_count": summary_row.get("fingerprint_label_count", 0),
        "fingerprint_png": relpath(fingerprint_paths["png"]),
        "fingerprint_pdf": relpath(fingerprint_paths["pdf"]),
        "fingerprint_svg": relpath(fingerprint_paths["svg"]),
        "region_tracks_png": relpath(region_tracks_paths["png"])
        if region_tracks_paths
        else "",
        "region_tracks_pdf": relpath(region_tracks_paths["pdf"])
        if region_tracks_paths
        else "",
        "region_tracks_svg": relpath(region_tracks_paths["svg"])
        if region_tracks_paths
        else "",
        "review_note": review_note,
    }


def build_review_markdown(review_rows: List[Dict[str, Any]]) -> str:
    lines = [
        "# Paper Case Study Review",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
    ]

    for row in review_rows:
        lines.extend(
            [
                f"## {row['panel_label']}",
                f"Status: {row['status']}",
                f"Relevance: {row['relevance_review']}",
                f"Resolution: {row['resolution_review']} ({row['figure_dpi']} DPI PNG plus PDF and SVG)",
                f"Subfigure readiness: {row['subfigure_review']} ({row['label_count']} ranked labels)",
                f"PNG: {row.get('fingerprint_png', '')}",
                f"PDF: {row.get('fingerprint_pdf', '')}",
                f"SVG: {row.get('fingerprint_svg', '')}",
                f"Region PNG: {row.get('region_tracks_png', '')}",
                f"Notes: {row['review_note']}",
                "",
            ]
        )

    return "\n".join(lines).rstrip() + "\n"


def skipped_case_row(case: pd.Series, reason: str) -> Dict[str, Any]:
    return {
        "case_id": case["case_id"],
        "panel_label": case["panel_label"],
        "status": "skipped",
        "species": case["species"],
        "gene": case["gene"],
        "chrom": case["chrom"],
        "pos": to_builtin(case.get("pos")),
        "ref": case["ref"],
        "alt": case["alt"],
        "reason": reason,
    }


def render_case_summary_markdown(summary_row: Dict[str, Any]) -> str:
    lines = [
        f"# {summary_row['panel_label']}",
        "",
        f"Status: {summary_row['status']}",
        f"Variant: {summary_row['chrom']}:{summary_row['pos']} {summary_row['ref']}>{summary_row['alt']}",
        f"Gene: {summary_row['gene_name']}",
    ]

    if summary_row.get("transcript"):
        lines.append(f"Transcript: {summary_row['transcript']}")
    if summary_row.get("transcript_hgvs"):
        lines.append(f"Transcript HGVS: {summary_row['transcript_hgvs']}")
    if summary_row.get("protein_hgvs"):
        lines.append(f"Protein HGVS: {summary_row['protein_hgvs']}")

    lines.extend(
        [
            f"Region class: {summary_row['region_class']}",
            f"Variant type: {summary_row['variant_type']}",
            f"MAGI score: {summary_row.get('Global_z_sum_log')}",
            f"BED impact: {summary_row.get('Impact_Score_BED')}",
            f"BigWig impact: {summary_row.get('Impact_Score_BW')}",
            f"LLR: {summary_row.get('LLR')}",
            f"Top BED signal: {summary_row.get('top_bed_signal', '')} ({summary_row.get('top_bed_delta')})",
            f"Top BigWig signal: {summary_row.get('top_bigwig_signal', '')} ({summary_row.get('top_bigwig_delta')})",
            f"Top MLM signal: {summary_row.get('top_mlm_signal', '')} ({summary_row.get('top_mlm_value')})",
        ]
    )
    return "\n".join(lines) + "\n"


def run_case(
    case: pd.Series,
    modules: Dict[str, Any],
    metadata_df: pd.DataFrame,
    metadata_dict: Dict[str, Dict[str, str]],
    output_dir: Path,
    device: str,
    top_k: int,
    figure_dpi: int,
) -> tuple[Dict[str, Any], Dict[str, Any]]:
    inference_module = modules["inference"]
    analysis_module = modules["analysis"]
    annotation_module = modules["annotation"]
    interpretation_module = modules["interpretation"]
    tracks_module = modules["tracks"]

    case_dir = output_dir / case["case_id"]
    case_dir.mkdir(parents=True, exist_ok=True)

    species = str(case["species"] or "human").strip()
    chrom = normalize_chrom(case["chrom"])
    pos = int(case["pos"])
    ref = str(case["ref"]).upper().strip()
    alt = str(case["alt"]).upper().strip()

    input_df = pd.DataFrame([{"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}])
    results_df = inference_module.predict_variants(
        input_df,
        device=device,
        species=species,
        cache_profiles=True,
    )
    results_df = annotate_for_species(annotation_module, results_df, species)
    results_df = analysis_module.compute_impact_scores(results_df)

    row = results_df.iloc[0]
    bed_names = list(inference_module._MODEL_CACHE.get("bed_names") or [])
    bw_indices = list(inference_module._MODEL_CACHE.get("selected_bw_indices") or [])
    bw_names_all = list(inference_module._MODEL_CACHE.get("bigwig_names") or [])
    bw_names_filtered = (
        [bw_names_all[index] for index in bw_indices] if species == "human" else None
    )

    ranked = analysis_module.rank_top_disrupted_tracks(
        row,
        bed_names,
        bw_names_filtered,
        metadata_df=metadata_df,
        metadata_dict=metadata_dict,
        top_k=None,
    )
    ranked_for_outputs = prepare_ranked_for_outputs(
        ranked,
        disambiguate_top_n=max(20, top_k),
    )
    top_tracks_df = analysis_module.build_top_track_table(
        ranked_for_outputs,
        max_rows=max(20, top_k),
        min_rows_by_type={"BED": 5},
    )
    fingerprint_fig = analysis_module.make_fingerprint_plot(
        ranked_for_outputs,
        top_k=top_k,
        figsize=(11, 6.5),
    )
    polish_fingerprint_figure(fingerprint_fig, case["panel_label"])
    region_tracks_fig = tracks_module.generate_region_tracks_plot(
        ranked_tracks=ranked_for_outputs,
        max_ranked_tracks=min(10, max(6, top_k)),
        visible_radius_bp=1000,
    )

    variant_type = "Indel" if int(row.get("indel_size", 0) or 0) != 0 else "SNP"
    interpretation_md = interpretation_module.build_signal_interpretation(
        row,
        ranked_for_outputs,
        variant_type,
    )
    summary_signals = analysis_module.extract_top_summary_signals(
        row,
        ranked_for_outputs,
        min_abs_threshold=0.03,
    )

    fingerprint_paths = save_figure(
        fingerprint_fig, case_dir / "fingerprint", figure_dpi
    )
    analysis_module.plt.close(fingerprint_fig)
    region_tracks_paths = None
    if region_tracks_fig is not None:
        region_tracks_paths = save_figure(
            region_tracks_fig,
            case_dir / "region_tracks",
            figure_dpi,
        )
        analysis_module.plt.close(region_tracks_fig)

    ranked_df = pd.DataFrame(ranked_for_outputs)
    provenance = build_case_provenance(case)

    top_bed_signal, top_bed_delta = top_signal(summary_signals.get("bed", []))
    top_bigwig_signal, top_bigwig_delta = top_signal(summary_signals.get("bigwig", []))
    top_mlm_signal, top_mlm_value = top_signal(summary_signals.get("mlm", []))
    top_abs_delta = float(abs(ranked[0]["delta"])) if ranked else None

    summary_row = {
        "case_id": case["case_id"],
        "panel_label": case["panel_label"],
        "status": "generated",
        "species": species,
        "gene": case["gene"],
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "transcript": case["transcript"],
        "transcript_hgvs": case["transcript_hgvs"],
        "protein_hgvs": case["protein_hgvs"],
        "gene_name": row.get("gene_name", case["gene"]),
        "region_class": row.get("region_class", ""),
        "variant_type": variant_type,
        "Impact_Score_BED": to_builtin(row.get("Impact_Score_BED")),
        "Impact_Score_BW": to_builtin(row.get("Impact_Score_BW")),
        "Global_z_sum_log": to_builtin(row.get("Global_z_sum_log")),
        "LLR": to_builtin(row.get("LLR")),
        "MLM_KL_mean": to_builtin(row.get("MLM_KL_mean")),
        "MLM_KL_max": to_builtin(row.get("MLM_KL_max")),
        "top_abs_delta": top_abs_delta,
        "top_bed_signal": top_bed_signal,
        "top_bed_delta": top_bed_delta,
        "top_bigwig_signal": top_bigwig_signal,
        "top_bigwig_delta": top_bigwig_delta,
        "top_mlm_signal": top_mlm_signal,
        "top_mlm_value": top_mlm_value,
        "fingerprint_label_count": min(top_k, len(ranked_for_outputs)),
        "fingerprint_png": relpath(fingerprint_paths["png"]),
        "fingerprint_pdf": relpath(fingerprint_paths["pdf"]),
        "fingerprint_svg": relpath(fingerprint_paths["svg"]),
        "region_tracks_png": relpath(region_tracks_paths["png"])
        if region_tracks_paths
        else "",
        "region_tracks_pdf": relpath(region_tracks_paths["pdf"])
        if region_tracks_paths
        else "",
        "region_tracks_svg": relpath(region_tracks_paths["svg"])
        if region_tracks_paths
        else "",
        "variant_result_csv": relpath(case_dir / "variant_result.csv"),
        "top_tracks_csv": relpath(case_dir / "top_tracks.csv"),
        "ranked_tracks_csv": relpath(case_dir / "ranked_tracks.csv"),
        "interpretation_md": relpath(case_dir / "interpretation.md"),
        "summary_md": relpath(case_dir / "case_summary.md"),
        "summary_json": relpath(case_dir / "summary.json"),
        "provenance_json": relpath(case_dir / "provenance.json"),
        "manuscript_note": case["manuscript_note"],
    }

    make_csv_safe(results_df).to_csv(case_dir / "variant_result.csv", index=False)
    top_tracks_df.to_csv(case_dir / "top_tracks.csv", index=False)
    ranked_df.to_csv(case_dir / "ranked_tracks.csv", index=False)
    (case_dir / "interpretation.md").write_text(
        interpretation_md + "\n", encoding="utf-8"
    )
    (case_dir / "case_summary.md").write_text(
        render_case_summary_markdown(summary_row),
        encoding="utf-8",
    )
    (case_dir / "summary.json").write_text(
        json.dumps(to_builtin(summary_row), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (case_dir / "provenance.json").write_text(
        json.dumps(to_builtin(provenance), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    review_row = build_review_row(
        case,
        summary_row,
        top_k,
        figure_dpi,
        fingerprint_paths,
        region_tracks_paths=region_tracks_paths,
    )
    return summary_row, review_row


def main() -> int:
    args = parse_args()
    manifest = load_manifest(args.manifest)
    selected_case_ids = set(args.case_ids or [])

    modules = load_magi_modules()
    inference_module = modules["inference"]
    annotation_module = modules["annotation"]
    device = choose_device(inference_module, args.device)
    metadata_df, metadata_dict = load_track_metadata()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading MAGI resources on {device}...")
    inference_module.load_model_and_resources(device=device)
    annotation_module.load_mane_data()

    summaries: List[Dict[str, Any]] = []
    reviews: List[Dict[str, Any]] = []

    for _, case in manifest.iterrows():
        if selected_case_ids and case["case_id"] not in selected_case_ids:
            continue

        if not is_resolved_case(case):
            reason = (
                "pending exact variant input"
                if str(case.get("status", "")).startswith("pending")
                else "missing coordinates"
            )
            summaries.append(skipped_case_row(case, reason))
            reviews.append(
                {
                    "case_id": case["case_id"],
                    "panel_label": case["panel_label"],
                    "status": "skipped",
                    "resolution_review": "pending",
                    "subfigure_review": "pending",
                    "relevance_review": "pending",
                    "figure_dpi": args.figure_dpi,
                    "figure_size_inches": "10x6",
                    "label_count": 0,
                    "fingerprint_png": "",
                    "fingerprint_pdf": "",
                    "fingerprint_svg": "",
                    "review_note": reason,
                }
            )
            print(f"Skipping {case['case_id']}: {reason}")
            continue

        print(
            f"Generating assets for {case['case_id']} ({case['chrom']}:{int(case['pos'])} {case['ref']}>{case['alt']})..."
        )
        try:
            summary_row, review_row = run_case(
                case,
                modules,
                metadata_df,
                metadata_dict,
                args.output_dir,
                device,
                args.top_k,
                args.figure_dpi,
            )
            summaries.append(summary_row)
            reviews.append(review_row)
        except Exception as exc:
            case_dir = args.output_dir / case["case_id"]
            case_dir.mkdir(parents=True, exist_ok=True)
            error_path = case_dir / "error.txt"
            error_path.write_text(traceback.format_exc(), encoding="utf-8")
            summaries.append(
                {
                    "case_id": case["case_id"],
                    "panel_label": case["panel_label"],
                    "status": "error",
                    "species": case["species"],
                    "gene": case["gene"],
                    "chrom": case["chrom"],
                    "pos": to_builtin(case.get("pos")),
                    "ref": case["ref"],
                    "alt": case["alt"],
                    "reason": str(exc),
                    "error_file": relpath(error_path),
                }
            )
            reviews.append(
                {
                    "case_id": case["case_id"],
                    "panel_label": case["panel_label"],
                    "status": "error",
                    "resolution_review": "failed",
                    "subfigure_review": "failed",
                    "relevance_review": "failed",
                    "figure_dpi": args.figure_dpi,
                    "figure_size_inches": "10x6",
                    "label_count": 0,
                    "fingerprint_png": "",
                    "fingerprint_pdf": "",
                    "fingerprint_svg": "",
                    "review_note": str(exc),
                }
            )
            print(f"Failed for {case['case_id']}: {exc}")

    summary_df = pd.DataFrame(summaries)
    review_df = pd.DataFrame(reviews)

    summary_df.to_csv(args.output_dir / "case_study_summary.csv", index=False)
    review_df.to_csv(args.output_dir / "review.csv", index=False)
    (args.output_dir / "review.md").write_text(
        build_review_markdown(reviews),
        encoding="utf-8",
    )

    print(f"Wrote summary to {args.output_dir / 'case_study_summary.csv'}")
    print(f"Wrote review to {args.output_dir / 'review.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
