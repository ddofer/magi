"""ClinVar rationale merge and strict filtering."""

from __future__ import annotations

from typing import Optional

import pandas as pd

from clinvar.constants import (
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
