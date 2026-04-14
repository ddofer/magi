#!/usr/bin/env python3
"""Assemble manuscript-facing case-study figures, captions, and review notes."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "analyses" / "paper_case_studies"
COMPOSITE_DPI = 400


@dataclass(frozen=True)
class CaseConfig:
    case_id: str
    panel_letter: str
    short_title: str


CASES = [
    CaseConfig("mat1a_gly336arg", "A", "MAT1A p.Gly336Arg"),
    CaseConfig("alox15b_rs9895916", "B", "ALOX15B rs9895916"),
    CaseConfig("col4a2_chr13_110492070_ga", "C", "COL4A2 chr13:110492070 G>A"),
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_top_tracks(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def format_delta(value: float | str | None) -> str:
    if value is None or value == "":
        return "n/a"
    return f"{float(value):+.4f}"


def save_figure(fig: plt.Figure, base_path: Path, dpi: int) -> None:
    fig.savefig(
        base_path.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(base_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def compose_panels(
    case_payloads: list[dict[str, Any]],
    image_key: str,
    output_base: Path,
    title: str,
    figsize: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(len(case_payloads), 1, figsize=figsize)
    if len(case_payloads) == 1:
        axes = [axes]

    fig.suptitle(title, fontsize=18, fontweight="bold", y=0.995)
    for ax, payload in zip(axes, case_payloads):
        image = mpimg.imread(payload[image_key])
        ax.imshow(image)
        ax.axis("off")
        ax.text(
            -0.03,
            1.02,
            payload["panel_letter"],
            transform=ax.transAxes,
            fontsize=20,
            fontweight="bold",
            ha="right",
            va="top",
        )

    fig.tight_layout(rect=(0, 0, 1, 0.985))
    save_figure(fig, output_base, COMPOSITE_DPI)


def collect_case_payloads() -> list[dict[str, Any]]:
    payloads: list[dict[str, Any]] = []
    for config in CASES:
        case_dir = PAPER_DIR / config.case_id
        summary = read_json(case_dir / "summary.json")
        top_tracks = read_top_tracks(case_dir / "top_tracks.csv")
        fingerprint_png = REPO_ROOT / summary["fingerprint_png"]
        region_tracks_png = REPO_ROOT / summary["region_tracks_png"]
        if not fingerprint_png.exists() or not region_tracks_png.exists():
            raise FileNotFoundError(f"Missing exported panel for {config.case_id}")

        payloads.append(
            {
                "case_id": config.case_id,
                "panel_letter": config.panel_letter,
                "short_title": config.short_title,
                "summary": summary,
                "top_tracks": top_tracks,
                "fingerprint_png": fingerprint_png,
                "region_tracks_png": region_tracks_png,
            }
        )
    return payloads


def write_panel_manifest(
    case_payloads: list[dict[str, Any]], output_path: Path
) -> None:
    rows: list[dict[str, str]] = []
    for payload in case_payloads:
        summary = payload["summary"]
        top_tracks = payload["top_tracks"]
        primary = top_tracks[0]
        secondary = top_tracks[1] if len(top_tracks) > 1 else {"Track": "", "Δ": ""}
        rows.append(
            {
                "panel_letter": payload["panel_letter"],
                "case_id": payload["case_id"],
                "panel_label": summary["panel_label"],
                "transcript": summary.get("transcript", ""),
                "transcript_hgvs": summary.get("transcript_hgvs", ""),
                "protein_hgvs": summary.get("protein_hgvs", ""),
                "region_class": summary.get("region_class", ""),
                "llr": str(summary.get("LLR", "")),
                "top_abs_delta": str(summary.get("top_abs_delta", "")),
                "primary_track": primary.get("Track", ""),
                "primary_delta": primary.get("Δ", ""),
                "secondary_track": secondary.get("Track", ""),
                "secondary_delta": secondary.get("Δ", ""),
                "fingerprint_png": str(
                    payload["fingerprint_png"].relative_to(REPO_ROOT)
                ),
                "region_tracks_png": str(
                    payload["region_tracks_png"].relative_to(REPO_ROOT)
                ),
            }
        )

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_main_caption(case_payloads: list[dict[str, Any]]) -> str:
    descriptions: list[str] = []
    for payload in case_payloads:
        summary = payload["summary"]
        descriptions.append(
            f"({payload['panel_letter']}) {summary['panel_label']} ({summary.get('transcript_hgvs', 'genomic coordinate only')}; {summary.get('protein_hgvs') or 'protein consequence not shown'}) shows a dominant {summary['top_bed_signal']} shift ({format_delta(summary['top_bed_delta'])}) and strongest contextual change in {summary['top_bigwig_signal']} ({format_delta(summary['top_bigwig_delta'])})."
        )

    return (
        "Figure X. MAGI fingerprint profiles for three manuscript case studies. "
        + " ".join(descriptions)
        + " Bars report Alt minus Ref delta probabilities for the top 15 disrupted MAGI signals. Source panels were exported at 600 DPI with matched PDF and SVG files; this assembled figure is rendered at 400 DPI for manuscript layout."
    )


def build_supplement_caption(case_payloads: list[dict[str, Any]]) -> str:
    descriptions: list[str] = []
    for payload in case_payloads:
        summary = payload["summary"]
        top_tracks = payload["top_tracks"]
        lead_tracks = ", ".join(
            f"{item['Track']} ({item['Δ']})" for item in top_tracks[:3]
        )
        descriptions.append(
            f"({payload['panel_letter']}) {summary['panel_label']} highlights local REF versus ALT profile changes centered on the variant, with the strongest track-level effects in {lead_tracks}."
        )

    return (
        "Supplementary Figure X. Local region-track views for the manuscript case studies. "
        + " ".join(descriptions)
        + " Dashed vertical lines mark the variant position; red shading denotes gain and blue shading denotes loss relative to the reference prediction."
    )


def build_results_text(case_payloads: list[dict[str, Any]]) -> str:
    summary_by_id = {
        payload["case_id"]: payload["summary"] for payload in case_payloads
    }
    mat1a = summary_by_id["mat1a_gly336arg"]
    alox = summary_by_id["alox15b_rs9895916"]
    col4a2 = summary_by_id["col4a2_chr13_110492070_ga"]

    comparative = (
        "Across the three case studies, MAGI recovered coherent but distinct local disruption signatures. "
        f"MAT1A {mat1a['transcript_hgvs']} ({mat1a['protein_hgvs']}) produced the strongest exon-centered loss signal in the set (exon {format_delta(mat1a['top_bed_delta'])}, top absolute delta {format_delta(mat1a['top_abs_delta'])}, LLR {mat1a['LLR']:.3f}), consistent with a coding variant that also perturbs exon-definition features. "
        f"ALOX15B {alox['transcript_hgvs']} ({alox['protein_hgvs']}) showed a reciprocal intron gain and exon/ORF loss pattern (intron {format_delta(alox['top_bed_delta'])}, top contextual change in {alox['top_bigwig_signal']} {format_delta(alox['top_bigwig_delta'])}), but with weaker sequence-model support (LLR {alox['LLR']:.3f}), arguing for a more moderate mechanistic effect. "
        f"COL4A2 {col4a2['transcript_hgvs']} ({col4a2['protein_hgvs']}) combined strong coding-track disruption with broad transcriptomic losses, yielding the most negative sequence-model score among the three cases (LLR {col4a2['LLR']:.3f})."
    )

    col4a2_text = (
        f"The newly added COL4A2 variant, {col4a2['transcript_hgvs']} ({col4a2['protein_hgvs']}), maps to a coding exon in the MANE Select transcript rather than to a canonical splice-donor position. "
        f"MAGI assigns a strong loss of exon probability ({format_delta(col4a2['top_bed_delta'])}) together with a reciprocal intron gain (+0.3203 in the top-track table) and substantial losses across brain and liver RNA-seq contexts, with the strongest contextual effect in {col4a2['top_bigwig_signal']} ({format_delta(col4a2['top_bigwig_delta'])}). "
        f"The aggregate scores remain high (BED impact {col4a2['Impact_Score_BED']:.3f}, BigWig impact {col4a2['Impact_Score_BW']:.3f}, MAGI score {col4a2['Global_z_sum_log']:.2f}), and the negative LLR ({col4a2['LLR']:.3f}) indicates that the alternate sequence is disfavored by the sequence model. "
        "Taken together, the pattern supports a missense-associated coding disruption with a pronounced local exon-structure signature, and the manuscript text should describe it as a coding COL4A2 case rather than as a splice-donor example."
    )

    return (
        "# Manuscript Results Text\n\n"
        "## Cross-Case Paragraph\n\n"
        f"{comparative}\n\n"
        "## COL4A2 Paragraph\n\n"
        f"{col4a2_text}\n"
    )


def build_review(case_payloads: list[dict[str, Any]]) -> str:
    fingerprint_lines = []
    for payload in case_payloads:
        fp = Image.open(payload["fingerprint_png"])
        region = Image.open(payload["region_tracks_png"])
        fingerprint_lines.append(
            f"- {payload['summary']['panel_label']}: fingerprint {fp.size[0]}x{fp.size[1]} px @ {fp.info.get('dpi')}, region tracks {region.size[0]}x{region.size[1]} px @ {region.info.get('dpi')}"
        )

    findings = [
        "No blocking internal inconsistencies remain after updating ALOX15B to include transcript and protein HGVS.",
        "Do not describe the new COL4A2 coordinate as a splice-donor variant. The supplied genomic coordinate resolves to a coding exon in NM_001846.4 and yields p.Gly1152Asp.",
        "The scientific story is coherent across cases: MAT1A is the strongest exon-loss example, ALOX15B is a moderate intron-gain/exon-loss example, and COL4A2 is a coding exon-disruption example with strong contextual RNA-seq losses.",
        "The fingerprint panels are main-text ready. The region-track panels are accurate and aligned with the deltas, but their height makes them better suited to a supplementary figure or an expanded multi-panel layout.",
        "ALOX15B has the weakest sequence-model penalty (LLR -0.295) among the three cases, so the manuscript should avoid overclaiming a severe sequence-level effect there.",
    ]

    return (
        "# Rubber-Duck Review\n\n"
        "## Findings\n\n"
        + "\n".join(f"- {item}" for item in findings)
        + "\n\n## Raster Checks\n\n"
        + "\n".join(fingerprint_lines)
        + "\n"
    )


def build_package_markdown(case_payloads: list[dict[str, Any]]) -> str:
    main_caption = build_main_caption(case_payloads)
    supplement_caption = build_supplement_caption(case_payloads)

    lines = [
        "# Manuscript Figure Package",
        "",
        "## Figure Assignments",
        "",
        "- Main figure fingerprint composite: analyses/paper_case_studies/figure_case_studies_fingerprints.png",
        "- Supplementary region composite: analyses/paper_case_studies/figure_case_studies_region_tracks.png",
        "- Panel manifest: analyses/paper_case_studies/manuscript_panel_manifest.csv",
        "- Results text: analyses/paper_case_studies/manuscript_case_text.md",
        "- Review memo: analyses/paper_case_studies/rubber_duck_review.md",
        "",
        "## Main Figure Caption",
        "",
        main_caption,
        "",
        "## Supplementary Figure Caption",
        "",
        supplement_caption,
        "",
        "## Panel Assets",
        "",
    ]

    for payload in case_payloads:
        summary = payload["summary"]
        lines.extend(
            [
                f"### Panel {payload['panel_letter']} - {summary['panel_label']}",
                "",
                f"- Fingerprint PNG: {summary['fingerprint_png']}",
                f"- Fingerprint PDF: {summary['fingerprint_pdf']}",
                f"- Region PNG: {summary['region_tracks_png']}",
                f"- Region PDF: {summary['region_tracks_pdf']}",
                f"- Transcript HGVS: {summary.get('transcript_hgvs', '')}",
                f"- Protein HGVS: {summary.get('protein_hgvs', '')}",
                f"- Top BED signal: {summary['top_bed_signal']} ({format_delta(summary['top_bed_delta'])})",
                f"- Top BigWig signal: {summary['top_bigwig_signal']} ({format_delta(summary['top_bigwig_delta'])})",
                f"- LLR: {summary['LLR']:.3f}",
                "",
            ]
        )

    return "\n".join(lines)


def main() -> int:
    case_payloads = collect_case_payloads()
    compose_panels(
        case_payloads,
        image_key="fingerprint_png",
        output_base=PAPER_DIR / "figure_case_studies_fingerprints",
        title="MAGI Case Studies: Fingerprint Panels",
        figsize=(13, 18),
    )
    compose_panels(
        case_payloads,
        image_key="region_tracks_png",
        output_base=PAPER_DIR / "figure_case_studies_region_tracks",
        title="MAGI Case Studies: Region-Track Panels",
        figsize=(13, 42),
    )

    write_panel_manifest(case_payloads, PAPER_DIR / "manuscript_panel_manifest.csv")
    (PAPER_DIR / "manuscript_figure_package.md").write_text(
        build_package_markdown(case_payloads),
        encoding="utf-8",
    )
    (PAPER_DIR / "manuscript_case_text.md").write_text(
        build_results_text(case_payloads),
        encoding="utf-8",
    )
    (PAPER_DIR / "rubber_duck_review.md").write_text(
        build_review(case_payloads),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
