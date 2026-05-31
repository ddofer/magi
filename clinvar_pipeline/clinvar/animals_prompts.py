"""System and per-variant LLM prompts for OMIA animal concordance evaluation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_PROMPTS_DIR = Path(__file__).resolve().parent / "prompts"


def _load_prompt(name: str) -> str:
    return (_PROMPTS_DIR / name).read_text(encoding="utf-8")


def get_system_prompt(variant_type: str = "snp") -> str:
    if variant_type == "indel":
        return _load_prompt("system_prompt_animals_indel.txt")
    return _load_prompt("system_prompt_animals_snp.txt")


def _safe_get(row: pd.Series, col: str, default: str = "(not available)") -> str:
    value = row.get(col, default)
    if pd.isna(value) or value == "":
        return default
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def make_variant_prompt_snp(row: pd.Series) -> str:
    variant_id = row.name if row.name is not None else _safe_get(row, "#VariationID", "Unknown")
    rationale = _safe_get(row, "FullRationale")
    if rationale in {"(not available)", ""}:
        rationale = "(No detailed rationale provided)"
    elif len(rationale) > 1500:
        rationale = rationale[:1500] + "... [truncated]"

    bed_abs = _safe_get(row, "BED_Top_Abs", "None")
    bed_gains = _safe_get(row, "BED_Top_Gains", "None")
    bed_losses = _safe_get(row, "BED_Top_Losses", "None")
    mlm_summary = _safe_get(row, "MLM_Summary", "None")

    return f"""## Variant: {variant_id}

### ClinVar Rationale
{rationale}

### BED Feature Signals (Genomic Annotations)
**Top by magnitude:** {bed_abs}
**Top gains (positive delta):** {bed_gains}
**Top losses (negative delta):** {bed_losses}


### MLM Sequence Context
{mlm_summary}

---
Evaluate concordance between the NT signals and the ClinVar rationale. Return JSON only."""


def make_variant_prompt_indel(row: pd.Series) -> str:
    variant_id = row.name if row.name is not None else _safe_get(row, "#VariationID", "Unknown")
    indel_size = row.get("indel_size", "Unknown")
    variant_type = row.get("variant_type", "Unknown")
    rationale = _safe_get(row, "FullRationale")
    if rationale in {"(not available)", ""}:
        rationale = "(No detailed rationale provided)"
    elif len(rationale) > 1500:
        rationale = rationale[:1500] + "... [truncated]"

    bed_abs = _safe_get(row, "BED_Top_Abs", "None")
    bed_gains = _safe_get(row, "BED_Top_Gains", "None")
    bed_losses = _safe_get(row, "BED_Top_Losses", "None")
    mlm_summary = _safe_get(row, "MLM_Summary", "None")
    mlm_logprob_delta = _safe_get(row, "MLM_logprob_delta", "N/A")
    emb_cosine = _safe_get(row, "EMB_cosine_dist", "N/A")
    emb_l2 = _safe_get(row, "EMB_l2_dist", "N/A")
    kl_max = _safe_get(row, "MLM_KL_max", "N/A")
    kl_mean = _safe_get(row, "MLM_KL_mean", "N/A")

    return f"""## Variant: {variant_id}
**Variant Type:** {variant_type}
**Variant Size:** {indel_size}

### ClinVar Rationale
{rationale}

### BED Feature Signals (Genomic Annotations)
**Top by magnitude:** {bed_abs}
**Top gains (positive delta):** {bed_gains}
**Top losses (negative delta):** {bed_losses}


### MLM Sequence Context
**Summary:** {mlm_summary}

**Detailed Indel Metrics:**
- Log-probability delta: {mlm_logprob_delta}
- Embedding cosine distance: {emb_cosine}
- Embedding L2 distance: {emb_l2}
- KL divergence (max): {kl_max}
- KL divergence (mean): {kl_mean}

---
Evaluate concordance between the NT signals and the ClinVar rationale. Return JSON only."""


def make_variant_prompt(row: pd.Series, variant_type: str = "snp") -> str:
    if variant_type == "indel":
        return make_variant_prompt_indel(row)
    return make_variant_prompt_snp(row)


def build_prompts_table(df: pd.DataFrame, variant_type: str = "snp") -> pd.DataFrame:
    out = df.copy()
    out["llm_prompt"] = out.apply(lambda row: make_variant_prompt(row, variant_type), axis=1)
    return out
