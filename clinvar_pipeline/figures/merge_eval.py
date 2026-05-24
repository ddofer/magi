"""Merge gathered LLM evaluations with signaled cohort + impact scores."""

from __future__ import annotations

import pandas as pd

from config import IMPACT_SCORE_COLS, PATHS
from clinvar.impact_scoring import attach_impact_to_signaled
from figures.mechanism_canonicalize import (
    build_primary_mechanism_canonicalizer,
    build_rationale_mechanism_canonicalizer,
)


def load_merged_eval_frame() -> pd.DataFrame:
    snps_d = attach_impact_to_signaled(PATHS["snp_annotated_signaled"], PATHS["snp_deltas_impact"])
    indels_d = attach_impact_to_signaled(PATHS["indel_annotated_signaled"], PATHS["indel_deltas_impact"])
    if "variant_id" in indels_d.columns and "#VariationID" not in indels_d.columns:
        indels_d = indels_d.rename(columns={"variant_id": "#VariationID"})

    snps = pd.read_parquet(PATHS["snp_llm_results"])
    indels = pd.read_parquet(PATHS["indel_llm_results"])
    snps["variant_type"] = "snp"
    indels["variant_type"] = "indel"

    impact_cols = list(dict.fromkeys(["#VariationID", "label", *IMPACT_SCORE_COLS, "LLR"]))
    delta_cols = [c for c in impact_cols if c in snps_d.columns]

    merged = snps.merge(snps_d[delta_cols], on="#VariationID", how="left")

    indel_llr = indels_d[["#VariationID", "LLR"]].rename(columns={"LLR": "LLR_indel"})
    merged = merged.merge(indel_llr, on="#VariationID", how="left")
    merged["LLR"] = merged["LLR"].fillna(merged["LLR_indel"])
    merged = merged.drop(columns=["LLR_indel"])

    impact_no_llr = [c for c in IMPACT_SCORE_COLS if c != "LLR" and c in indels_d.columns]
    double_cols = ["label", *impact_no_llr]
    if double_cols:
        ind_tmp = indels_d[["#VariationID", *double_cols]].rename(
            columns={c: f"{c}_indel" for c in double_cols}
        )
        merged = merged.merge(ind_tmp, on="#VariationID", how="left")
        for c in double_cols:
            if c in merged.columns:
                merged[c] = merged[c].fillna(merged[f"{c}_indel"])
        merged = merged.drop(columns=[f"{c}_indel" for c in double_cols])

    indel_part = indels.merge(
        indels_d[[c for c in delta_cols if c in indels_d.columns]],
        on="#VariationID",
        how="left",
    )
    merged = pd.concat([merged, indel_part], ignore_index=True)

    merged["concordance"] = merged["concordance"].astype(str).str.strip().str.upper()

    primary_canon = build_primary_mechanism_canonicalizer()
    rationale_canon = build_rationale_mechanism_canonicalizer()
    if "primary_signal_mechanism" in merged.columns:
        merged["primary_signal_mechanism"] = merged["primary_signal_mechanism"].map(primary_canon)
    if "rationale_mechanism" in merged.columns:
        merged["rationale_mechanism"] = merged["rationale_mechanism"].map(rationale_canon)
    return merged
