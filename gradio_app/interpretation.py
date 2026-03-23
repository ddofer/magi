#!/usr/bin/env python3
"""
Rule-based signal interpretation for the MAGI Gradio app.

This module is intentionally self-contained. It converts already-available
single-variant outputs into a short interpretation panel without depending on
notebook code or outer-folder runtime imports.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


STRONG_DELTA = 0.20
MODERATE_DELTA = 0.10
WEAK_DELTA = 0.05

LLR_STRONG_NEG = -3.0
LLR_MODERATE_NEG = -1.5
LLR_POSITIVE = 1.0

KL_STRONG = 0.30
KL_MODERATE = 0.15

BED_FEATURE_LABELS: Dict[str, Tuple[str, str]] = {
    "mRNA_splice": ("splice-site recognition", "splice disruption"),
    "coding_sequence": ("coding sequence identity", "coding disruption"),
    "mRNA_exon": ("exonic structure", "coding or exon-level disruption"),
    "start_codon": ("translation initiation", "start-codon disruption"),
    "stop_codon": ("translation termination", "stop-codon disruption"),
    "mRNA_promoter": ("promoter activity", "promoter-regulatory disruption"),
    "five_prime_UTR": ("5' UTR regulation", "UTR-level regulation change"),
    "three_prime_UTR": ("3' UTR regulation", "UTR-level regulation change"),
    "mRNA_intron": ("intronic transcript context", "intronic transcript disruption"),
    "gene": ("genic context", "genic structural disruption"),
    "other": ("annotated genomic context", "localized genomic disruption"),
}


def _is_missing(value) -> bool:
    return value is None or (isinstance(value, float) and np.isnan(value))


def _fmt_value(value, fmt: str = ".3f", na: str = "N/A") -> str:
    if value is None:
        return na
    try:
        if pd.isna(value):
            return na
    except TypeError:
        pass
    try:
        return format(float(value), fmt)
    except (TypeError, ValueError):
        return na


def _strength_from_delta(delta: float) -> str:
    magnitude = abs(delta)
    if magnitude >= STRONG_DELTA:
        return "strong"
    if magnitude >= MODERATE_DELTA:
        return "moderate"
    if magnitude >= WEAK_DELTA:
        return "subtle"
    return "weak"


def _direction_word(delta: float) -> str:
    return "gain" if delta > 0 else "loss"


def _top_ranked_by_type(ranked: List[Dict], track_type: str, k: int = 2) -> List[Dict]:
    return [item for item in ranked if item.get("track_type") == track_type][:k]


def _humanize_bw_context(display_name: str) -> str:
    parts = [part.strip() for part in str(display_name).split("|") if part.strip()]
    if not parts:
        return str(display_name)
    if len(parts) == 1:
        return parts[0]
    if len(parts) == 2:
        return f"{parts[0]} ({parts[1]})"
    return f"{parts[0]} ({parts[1]}, {parts[2]})"


def _bed_mechanism(track_id: str) -> Tuple[str, str]:
    return BED_FEATURE_LABELS.get(
        track_id,
        (track_id.replace("_", " "), "localized structural disruption"),
    )


def _primary_mechanism(
    ranked: List[Dict], row: pd.Series, variant_type: str
) -> Tuple[str, str]:
    bed_ranked = _top_ranked_by_type(ranked, "BED", k=5)
    for item in bed_ranked:
        if abs(float(item.get("delta", 0.0))) < WEAK_DELTA:
            continue
        label, mechanism = _bed_mechanism(str(item.get("track_id", "other")))
        return mechanism, f"The strongest BED signal points to {label}."

    bw_ranked = _top_ranked_by_type(ranked, "BigWig", k=3)
    if bw_ranked and abs(float(bw_ranked[0].get("delta", 0.0))) >= MODERATE_DELTA:
        context = _humanize_bw_context(str(bw_ranked[0].get("display_name", "track")))
        return (
            "context-specific regulatory change",
            f"The strongest ranked context signal suggests a shift in {context}.",
        )

    llr = row.get("LLR", np.nan)
    kl_mean = row.get("MLM_KL_mean", np.nan)
    if not _is_missing(llr) and float(llr) <= LLR_MODERATE_NEG:
        return (
            "sequence-constraint signal",
            "The sequence model ranks the alternate sequence as less likely even without a dominant BED or BigWig signal.",
        )
    if not _is_missing(kl_mean) and float(kl_mean) >= KL_MODERATE:
        return (
            "sequence-disruption signal",
            "Token-level sequence distributions shift even though no single BED or BigWig track dominates.",
        )

    return (
        "mixed or weak evidence",
        f"No single {variant_type.lower()} mechanism dominates the current rule-based evidence.",
    )


def _bed_evidence_lines(ranked: List[Dict]) -> List[str]:
    lines: List[str] = []
    for item in _top_ranked_by_type(ranked, "BED", k=2):
        delta = float(item.get("delta", 0.0))
        label, _ = _bed_mechanism(str(item.get("track_id", "other")))
        strength = _strength_from_delta(delta)
        direction = _direction_word(delta)
        ref_val = _fmt_value(item.get("ref_val"))
        alt_val = _fmt_value(item.get("alt_val"))
        lines.append(
            f"BED `{item['track_id']}` shows a {strength} {direction} in {label} "
            f"(REF {ref_val} → ALT {alt_val}, Δ={delta:+.3f})."
        )
    return lines


def _bw_evidence_lines(ranked: List[Dict]) -> List[str]:
    lines: List[str] = []
    for item in _top_ranked_by_type(ranked, "BigWig", k=2):
        delta = float(item.get("delta", 0.0))
        strength = _strength_from_delta(delta)
        direction = _direction_word(delta)
        ref_val = _fmt_value(item.get("ref_val"))
        alt_val = _fmt_value(item.get("alt_val"))
        context = _humanize_bw_context(
            str(item.get("display_name", item.get("track_id", "track")))
        )
        lines.append(
            f"Context track `{context}` has a {strength} {direction} "
            f"(REF {ref_val} → ALT {alt_val}, Δ={delta:+.3f})."
        )
    return lines


def _sequence_evidence_lines(row: pd.Series, variant_type: str) -> List[str]:
    lines: List[str] = []

    llr = row.get("LLR", np.nan)
    if not _is_missing(llr):
        llr = float(llr)
        if llr <= LLR_STRONG_NEG:
            lines.append(
                f"Sequence-model evidence is strong: LLR {llr:.3f} makes the alternate sequence much less plausible than reference."
            )
        elif llr <= LLR_MODERATE_NEG:
            lines.append(
                f"Sequence-model evidence is supportive: LLR {llr:.3f} penalizes the alternate sequence."
            )
        elif llr >= LLR_POSITIVE:
            lines.append(
                f"LLR {llr:.3f} does not penalize the alternate allele, so sequence-only support is limited."
            )

    kl_mean = row.get("MLM_KL_mean", np.nan)
    if not _is_missing(kl_mean):
        kl_mean = float(kl_mean)
        if kl_mean >= KL_STRONG:
            lines.append(
                f"Mean KL {kl_mean:.3f} indicates a pronounced redistribution of token probabilities around the variant."
            )
        elif kl_mean >= KL_MODERATE:
            lines.append(
                f"Mean KL {kl_mean:.3f} indicates moderate local sequence perturbation."
            )

    if variant_type == "Indel":
        emb_cosine = row.get("EMB_cosine_dist", np.nan)
        emb_l2 = row.get("EMB_l2_dist", np.nan)
        if not _is_missing(emb_cosine) or not _is_missing(emb_l2):
            lines.append(
                "Indel embedding distances are available as supportive context: "
                f"cosine={_fmt_value(emb_cosine)}, L2={_fmt_value(emb_l2)}."
            )

    return lines


def build_signal_interpretation(
    row: pd.Series,
    ranked: List[Dict],
    variant_type: str,
) -> str:
    """Create a short deterministic interpretation panel."""
    gene_name = row.get("gene_name", "Intergenic") or "Intergenic"
    region_class = row.get("region_class", "OTHER")
    if str(gene_name).startswith("N/A (non-human)") or str(region_class) == "NON_HUMAN":
        context_anchor = "non-human BED + sequence outputs"
    else:
        context_anchor = f"{gene_name} / {region_class}"

    if not ranked:
        return (
            "### Rule-Based Signal Interpretation\n\n"
            "No ranked BED or BigWig signals available for rule-based interpretation.\n\n"
            "**Note:** The MAGI score is computed separately from bundled baseline statistics and is shown in the Variant Summary rather than in this panel."
        )

    mechanism, rationale = _primary_mechanism(ranked, row, variant_type)
    evidence_lines = []
    evidence_lines.extend(_bed_evidence_lines(ranked))
    evidence_lines.extend(_bw_evidence_lines(ranked))
    evidence_lines.extend(_sequence_evidence_lines(row, variant_type))

    if not evidence_lines:
        evidence_lines.append(
            "The current ranked outputs are weak, so this should be treated as a low-confidence summary only."
        )

    bullet_block = "\n".join(f"- {line}" for line in evidence_lines[:5])

    return f"""
### Rule-Based Signal Interpretation

**Primary hypothesis:** {mechanism.capitalize()}  
**Context anchor:** {context_anchor}  
**Why this is suggested:** {rationale}

**Top evidence**
{bullet_block}
""".strip()
