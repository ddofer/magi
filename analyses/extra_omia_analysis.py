"""Summarize OMIA post-inference outputs for manuscript tables.

Usage:
  python analyses/extra_omia_analysis.py
"""

from __future__ import annotations

from pathlib import Path
import json
import pandas as pd


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _find_species_col(df: pd.DataFrame) -> str | None:
    for c in df.columns:
        if c.lower() in {"species", "organism", "taxon"}:
            return c
    return None


def summarize_from_report(report: dict) -> pd.DataFrame:
    rows: list[dict] = []
    for item in report.get("load_report", []):
        rows.append(
            {
                "species": item.get("species"),
                "rows_loaded": item.get("rows_loaded"),
                "cols_loaded": item.get("cols_loaded"),
                "selected_file": item.get("selected_file"),
            }
        )
    return pd.DataFrame(rows)


def summarize_from_merged_parquet(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_parquet(path)
    species_col = _find_species_col(df)
    if species_col is None:
        return pd.DataFrame(
            [
                {
                    "species": "unknown",
                    "n_variants": len(df),
                    "mean_abs_D_BED": pd.NA,
                    "max_abs_D_BED": pd.NA,
                }
            ]
        )

    bed_cols = [c for c in df.columns if c.startswith("D_BED_")]
    out_rows: list[dict] = []
    for sp, g in df.groupby(species_col, dropna=False):
        row: dict[str, object] = {
            "species": str(sp),
            "n_variants": len(g),
        }
        if bed_cols:
            abs_vals = g[bed_cols].abs()
            row["mean_abs_D_BED"] = float(abs_vals.to_numpy().mean())
            row["max_abs_D_BED"] = float(abs_vals.to_numpy().max())
        else:
            row["mean_abs_D_BED"] = pd.NA
            row["max_abs_D_BED"] = pd.NA
        out_rows.append(row)

    return pd.DataFrame(out_rows).sort_values("species").reset_index(drop=True)


def main() -> None:
    repo = _repo_root()
    out_dir = repo / "analyses" / "generated_tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    report_path = repo / "data" / "results" / "omia_post_inference_report.json"
    summary_csv_path = repo / "data" / "results" / "omia_post_inference_summary.csv"
    merged_parquet_path = (
        repo / "data" / "results" / "omia_all_species_inferred_relevant.parquet"
    )

    report = _load_json(report_path)
    if report is not None:
        report_df = summarize_from_report(report)
        report_df.to_csv(out_dir / "omia_load_report_summary.csv", index=False)
        coverage = report.get("coverage", {})
        if isinstance(coverage, dict):
            coverage_rows = [{"species": k, "coverage": v} for k, v in coverage.items()]
        elif isinstance(coverage, list):
            coverage_rows = []
            for item in coverage:
                if isinstance(item, dict):
                    coverage_rows.append(
                        {
                            "species": item.get("species", "unknown"),
                            "coverage": item.get("coverage", item),
                        }
                    )
                else:
                    coverage_rows.append({"species": "unknown", "coverage": item})
        else:
            coverage_rows = [{"species": "unknown", "coverage": coverage}]

        coverage_df = pd.DataFrame(coverage_rows)
        coverage_df.to_csv(out_dir / "omia_coverage_summary.csv", index=False)

    if summary_csv_path.exists():
        summary_df = pd.read_csv(summary_csv_path)
        summary_df.to_csv(out_dir / "omia_post_inference_summary_copy.csv", index=False)

    merged_df = summarize_from_merged_parquet(merged_parquet_path)
    if not merged_df.empty:
        merged_df.to_csv(out_dir / "omia_species_feature_summary.csv", index=False)

    print("Wrote OMIA summaries to analyses/generated_tables/")


if __name__ == "__main__":
    main()
