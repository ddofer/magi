#!/usr/bin/env python3
"""Build two canonical LLM result files (SNP + indel) for figure generation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from clinvar.llm_results import INDEL_LLM_CORE, SNP_LLM_CORE
from config import DATA_ROOT, PATHS, PROJECT_ROOT


def _pick_columns(df: pd.DataFrame, wanted: list[str]) -> pd.DataFrame:
    cols = [c for c in wanted if c in df.columns]
    use = cols if cols else [c for c in df.columns if c not in ("raw_response", "parse_success")]
    return df[use].copy()


def _normalize_snp(df: pd.DataFrame) -> pd.DataFrame:
    if "primary_mechanism" in df.columns and "primary_signal_mechanism" not in df.columns:
        df["primary_signal_mechanism"] = df["primary_mechanism"]
    if "label" in df.columns and df["label"].dtype == bool:
        df["label"] = df["label"].map({True: "Pathogenic", False: "Benign"})
    return df


def build_snp_llm() -> pd.DataFrame:
    eval_path = PROJECT_ROOT / "parquet/snp_evaluation_results_v2.csv"
    search_path = PROJECT_ROOT / "snp_for_search.csv"
    if not eval_path.exists():
        raise FileNotFoundError(eval_path)

    df = pd.read_csv(eval_path)
    df = _normalize_snp(df)

    if search_path.exists():
        search = pd.read_csv(search_path)
        search = _normalize_snp(search)
        extras = [c for c in ("label", "primary_mechanism", "GeneSymbol") if c in search.columns]
        if extras:
            df = df.merge(
                search[["#VariationID", *extras]].drop_duplicates("#VariationID"),
                on="#VariationID",
                how="left",
                suffixes=("", "_search"),
            )
            for col in extras:
                sc = f"{col}_search"
                if sc in df.columns:
                    df[col] = df[col].fillna(df[sc])
                    df.drop(columns=[sc], inplace=True)

    return _pick_columns(df, SNP_LLM_CORE)


def build_indel_llm() -> pd.DataFrame:
    # v2 CSV is the canonical full-cohort export (12,301 variants; ~43% C+P).
    # Prefer it over llm_results/indel_results.parquet, which may be an older partial run.
    csv_path = PROJECT_ROOT / "parquet/indel_evaluation_results_v2.csv"
    parquet_path = PROJECT_ROOT / "llm_results/indel_results.parquet"

    if csv_path.exists():
        df = pd.read_csv(csv_path)
    elif parquet_path.exists():
        df = pd.read_parquet(parquet_path)
    else:
        raise FileNotFoundError(f"Need {csv_path} or {parquet_path}")

    if "variant_id" in df.columns and "#VariationID" not in df.columns:
        df = df.rename(columns={"variant_id": "#VariationID"})
    if "size" in df.columns and "indel_size" not in df.columns:
        df["indel_size"] = df["size"]
    if "indel_size" not in df.columns:
        sig_path = Path(PATHS["indel_annotated_signaled"])
        if sig_path.exists():
            sig = pd.read_parquet(sig_path, columns=["#VariationID", "indel_size"])
            df = df.merge(sig.drop_duplicates("#VariationID"), on="#VariationID", how="left")
    return _pick_columns(df, INDEL_LLM_CORE)


def gather(*, force: bool = False) -> None:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    (DATA_ROOT / "llm").mkdir(parents=True, exist_ok=True)
    snp_out = Path(PATHS["snp_llm_results"])
    indel_out = Path(PATHS["indel_llm_results"])

    if snp_out.exists() and not force:
        print(f"Skip existing {snp_out}")
    else:
        snp = build_snp_llm()
        snp.to_parquet(snp_out, index=False)
        print(f"  SNP  → {snp_out}  ({len(snp):,} rows, {len(snp.columns)} cols)")

    if indel_out.exists() and not force:
        print(f"Skip existing {indel_out}")
    else:
        indel = build_indel_llm()
        indel.to_parquet(indel_out, index=False)
        print(f"  INDEL → {indel_out}  ({len(indel):,} rows, {len(indel.columns)} cols)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    gather(force=args.force)


if __name__ == "__main__":
    main()
