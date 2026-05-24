"""ClinVar rationale merge and strict filtering."""

from __future__ import annotations

from typing import Optional

import pandas as pd

from clinvar.constants import (
    ALLOWED_OTHER_STATUS,
    CRITERIA_SINGLE,
    PLACEHOLDER_RATIONALES,
    SUBMISSION_COLS,
)


def _clean_text(text_series: pd.Series) -> str:
    valid = [
        str(t)
        for t in text_series.dropna()
        if str(t).strip() not in ("-", "")
    ]
    seen: set[str] = set()
    uniq: list[str] = []
    for v in valid:
        if v not in seen:
            uniq.append(v)
            seen.add(v)
    return " | ".join(uniq)


def _combine_rationale(row: pd.Series) -> str:
    parts = []
    if row.get("Description"):
        parts.append(f"DESC: {row['Description']}")
    if row.get("ExplanationOfInterpretation"):
        parts.append(f"EXPL: {row['ExplanationOfInterpretation']}")
    return "\n".join(parts) if parts else "No detailed rationale provided."


def aggregate_review_status(sub_df: pd.DataFrame) -> pd.Series:
    """Per-VariationID review-status aggregation (legacy; not used in strict filter)."""
    tmp = (
        sub_df[["#VariationID", "ReviewStatus"]]
        .dropna(subset=["#VariationID"])
        .copy()
    )
    tmp["#VariationID"] = tmp["#VariationID"].astype("Int64")
    tmp["ReviewStatus"] = tmp["ReviewStatus"].astype(str).str.strip().str.lower()

    allowed_all = {CRITERIA_SINGLE} | ALLOWED_OTHER_STATUS

    def _agg(s: pd.Series):
        vals = s.dropna().tolist()
        n = len(vals)
        if n == 0:
            return pd.NA
        if n == 1 and vals[0] == "reviewed by expert panel":
            return "reviewed_by_expert_panel"
        if n <= 1:
            return pd.NA
        if any(v not in allowed_all for v in vals):
            return pd.NA
        n_criteria = sum(v == CRITERIA_SINGLE for v in vals)
        if n_criteria < 1:
            return pd.NA
        n_other = n - n_criteria
        if n_criteria >= 2 or n_other >= 1:
            return "criteria_provided,_multiple_submitters,_no_conflicts"
        return pd.NA

    return tmp.groupby("#VariationID")["ReviewStatus"].apply(_agg)


def load_submission_summary(path: str) -> pd.DataFrame:
    print(f"  Loading submission summary from {path} ...")
    return pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=SUBMISSION_COLS,
        skiprows=18,
        low_memory=False,
    )


def fetch_rationales_and_filter(
    mapped_df: pd.DataFrame,
    submission_summary_path: str,
    min_gold_stars: int = 2,
    sub_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    if sub_df is None:
        sub_df = load_submission_summary(submission_summary_path)

    sub_df = sub_df.dropna(subset=["#VariationID"]).copy()
    mapped_df = mapped_df[mapped_df["gold_stars"] >= min_gold_stars].copy()

    print("  Aggregating rationales ...")
    grouped = (
        sub_df.groupby("#VariationID", as_index=False)
        .agg(
            {
                "Description": _clean_text,
                "ExplanationOfInterpretation": _clean_text,
            }
        )
    )
    grouped["FullRationale"] = grouped.apply(_combine_rationale, axis=1)
    grouped = grouped[["#VariationID", "FullRationale"]]

    print("  Merging rationales ...")
    final_df = mapped_df.copy()
    final_df["FullRationale"] = "No rationale provided."
    final_df["#VariationID"] = final_df["#VariationID"].astype("Int64")
    grouped["#VariationID"] = grouped["#VariationID"].astype("Int64")

    rationale_map = grouped.set_index("#VariationID")["FullRationale"]
    mask_has_id = final_df["#VariationID"].notna()
    final_df.loc[mask_has_id, "FullRationale"] = (
        final_df.loc[mask_has_id, "#VariationID"]
        .map(rationale_map)
        .fillna("No rationale provided.")
    )

    dup_mask = final_df.duplicated(subset=["chrom", "pos", "ref", "alt"], keep=False)
    dedup_df = final_df[~dup_mask].copy()
    print(
        f"  After dedup: {len(dedup_df):,}  "
        f"(dropped {dup_mask.sum():,} duplicate-locus rows)"
    )

    strict_df = dedup_df[
        ~dedup_df["FullRationale"].isin(PLACEHOLDER_RATIONALES)
    ].copy()
    print(f"  After strict filter: {len(strict_df):,}")
    funnel = {"quality": len(dedup_df), "rationaled": len(strict_df)}
    return strict_df, funnel
