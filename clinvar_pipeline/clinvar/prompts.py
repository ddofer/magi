"""System and per-variant LLM prompts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_PROMPTS_DIR = Path(__file__).resolve().parent / "prompts"


def _load_prompt(name: str) -> str:
    return (_PROMPTS_DIR / name).read_text(encoding="utf-8")


def get_system_prompt(variant_type: str = "snp") -> str:
    if variant_type == "indel":
        return _load_prompt("system_prompt_indel.txt")
    return _load_prompt("system_prompt_snp.txt")


def make_variant_prompt_snp(row: pd.Series) -> str:
    def safe_get(col, default="(not available)"):
        v = row.get(col, default)
        if pd.isna(v) or v == "":
            return default
        return str(v)

    variant_id = row.name if row.name else safe_get("#VariationID", "Unknown")
    region_class = safe_get("region_class", "Unknown")

    rationale = safe_get("FullRationale")
    if rationale == "(not available)" or rationale.strip() == "":
        rationale = "(No detailed rationale provided)"
    elif len(rationale) > 1500:
        rationale = rationale[:1500] + "... [truncated]"

    bed_abs = safe_get("BED_Top_Abs", "None")
    bed_gains = safe_get("BED_Top_Gains", "None")
    bed_losses = safe_get("BED_Top_Losses", "None")
    bw_abs = safe_get("BW_Top_Abs", "None")
    bw_gains = safe_get("BW_Top_Gains", "None")
    bw_losses = safe_get("BW_Top_Losses", "None")
    mlm_summary = safe_get("MLM_Summary", "None")

    return f"""## Variant: {variant_id}
**Region:** {region_class}

### ClinVar Rationale
{rationale}

### BED Feature Signals (Genomic Annotations)
**Top by magnitude:** {bed_abs}
**Top gains (positive delta):** {bed_gains}
**Top losses (negative delta):** {bed_losses}

### Epigenomic Track Signals (BW)
**Top by magnitude:** {bw_abs}
**Top gains:** {bw_gains}
**Top losses:** {bw_losses}

### MLM Sequence Context
{mlm_summary}

---
Evaluate concordance between the NT signals and the ClinVar rationale. Return JSON only."""


def make_variant_prompt_indel(row: pd.Series) -> str:
    def safe_get(col, default="(not available)"):
        v = row.get(col, default)
        if pd.isna(v) or v == "":
            return default
        if isinstance(v, float):
            return f"{v:.4f}"
        return str(v)

    variant_id = row.name if row.name else safe_get("#VariationID", "Unknown")
    region = safe_get("region", "Unknown")
    indel_size = row.get("indel_size", "Unknown")
    variant_type = row.get("variant_type", "Unknown")

    rationale = safe_get("FullRationale")
    if rationale == "(not available)" or rationale.strip() == "":
        rationale = "(No detailed rationale provided)"
    elif len(rationale) > 1500:
        rationale = rationale[:1500] + "... [truncated]"

    bed_abs = safe_get("BED_Top_Abs", "None")
    bed_gains = safe_get("BED_Top_Gains", "None")
    bed_losses = safe_get("BED_Top_Losses", "None")
    bw_abs = safe_get("BW_Top_Abs", "None")
    bw_gains = safe_get("BW_Top_Gains", "None")
    bw_losses = safe_get("BW_Top_Losses", "None")
    mlm_summary = safe_get("MLM_Summary", "None")
    mlm_logprob_delta = safe_get("MLM_logprob_delta", "N/A")
    emb_cosine = safe_get("EMB_cosine_dist", "N/A")
    emb_l2 = safe_get("EMB_l2_dist", "N/A")
    kl_max = safe_get("MLM_KL_max", "N/A")
    kl_mean = safe_get("MLM_KL_mean", "N/A")

    return f"""## Variant: {variant_id}
**Variant Type:** {variant_type}
**Variant Size:** {indel_size}
**Region:** {region}

### ClinVar Rationale
{rationale}

### BED Feature Signals (Genomic Annotations)
**Top by magnitude:** {bed_abs}
**Top gains (positive delta):** {bed_gains}
**Top losses (negative delta):** {bed_losses}

### Epigenomic Track Signals (BW)
**Top by magnitude:** {bw_abs}
**Top gains:** {bw_gains}
**Top losses:** {bw_losses}

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
    out["llm_prompt"] = out.apply(
        lambda row: make_variant_prompt(row, variant_type),
        axis=1,
    )
    return out
