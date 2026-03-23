"""Generate manuscript-ready summary tables for MAGI.

Usage:
  python analyses/figures_tables.py
"""

from __future__ import annotations

from pathlib import Path
import json
import pandas as pd


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _safe_read_csv(path: Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path)
    return None


def _safe_read_json(path: Path) -> dict | None:
    if path.exists():
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return None


def build_table1_dataset_summary(repo: Path) -> pd.DataFrame:
    rows: list[dict] = []

    # Concordance evaluation (human curated subset)
    eval_path = repo / "analyses" / "stav-llm_evaluation_results.csv"
    eval_df = _safe_read_csv(eval_path)
    if eval_df is not None and not eval_df.empty:
        label_col = None
        for c in eval_df.columns:
            if c.lower() == "label":
                label_col = c
                break
        n_total = len(eval_df)
        label_summary = "N/A"
        if label_col is not None:
            counts = eval_df[label_col].value_counts(dropna=False)
            label_summary = ", ".join([f"{k}: {v}" for k, v in counts.items()])

        rows.append(
            {
                "dataset": "ClinVar rationale evaluation",
                "species": "human",
                "n_variants": n_total,
                "labels": label_summary,
                "source_file": str(eval_path.relative_to(repo)),
            }
        )

    # OMIA summary
    omia_summary_path = repo / "data" / "results" / "omia_post_inference_summary.csv"
    omia_df = _safe_read_csv(omia_summary_path)
    if omia_df is not None and not omia_df.empty:
        cols_lower = {c.lower(): c for c in omia_df.columns}
        species_col = cols_lower.get("species")
        rows_col = cols_lower.get("rows")
        for _, r in omia_df.iterrows():
            species = str(r[species_col]) if species_col else "unknown"
            n_rows = int(r[rows_col]) if rows_col and pd.notna(r[rows_col]) else None
            rows.append(
                {
                    "dataset": "OMIA inferred",
                    "species": species,
                    "n_variants": n_rows,
                    "labels": "curated pathogenic (OMIA)",
                    "source_file": str(omia_summary_path.relative_to(repo)),
                }
            )

    if not rows:
        rows.append(
            {
                "dataset": "No dataset metadata found",
                "species": "N/A",
                "n_variants": "N/A",
                "labels": "N/A",
                "source_file": "N/A",
            }
        )

    return pd.DataFrame(rows)


def build_table2_feature_summary(repo: Path) -> pd.DataFrame:
    report_path = repo / "data" / "results" / "omia_post_inference_report.json"
    report = _safe_read_json(report_path)

    # Defaults from validated project notes
    bigwig_total = 7362
    bigwig_filtered = 7191
    kai = 222
    cage = 1276
    gtex = 84
    encode = 5609

    if report:
        coltax = report.get("column_taxonomy", {})
        # OMIA files do not always include BigWig, so keep defaults unless useful info exists.
        _ = coltax

    rows = [
        {
            "feature_group": "BED elements",
            "count": 21,
            "description": "Structural genomic annotations",
        },
        {
            "feature_group": "BigWig total (pre-filter)",
            "count": bigwig_total,
            "description": "All available functional tracks before filtering",
        },
        {
            "feature_group": "BigWig selected (post-filter)",
            "count": bigwig_filtered,
            "description": "Tracks retained for inference",
        },
        {
            "feature_group": "Catlas (kai*)",
            "count": kai,
            "description": "Single-cell accessibility tracks",
        },
        {
            "feature_group": "CAGE (CNhs*)",
            "count": cage,
            "description": "Transcription initiation tracks",
        },
        {
            "feature_group": "GTEx",
            "count": gtex,
            "description": "Tissue expression-derived tracks",
        },
        {
            "feature_group": "ENCODE (ENCSR*)",
            "count": encode,
            "description": "Histone/chromatin tracks from ENCODE",
        },
        {
            "feature_group": "MLM core metrics",
            "count": 8,
            "description": "LLR, priors, deltas, KL and log-prob features",
        },
    ]
    return pd.DataFrame(rows)


def build_table3_concordance(repo: Path) -> pd.DataFrame:
    eval_path = repo / "analyses" / "stav-llm_evaluation_results.csv"
    eval_df = _safe_read_csv(eval_path)
    if eval_df is None or eval_df.empty:
        return pd.DataFrame(
            [
                {
                    "concordance": "N/A",
                    "count": "N/A",
                    "percentage": "N/A",
                }
            ]
        )

    conc_col = None
    for c in eval_df.columns:
        if c.lower() == "concordance":
            conc_col = c
            break
    if conc_col is None:
        return pd.DataFrame(
            [
                {
                    "concordance": "Column not found",
                    "count": "N/A",
                    "percentage": "N/A",
                }
            ]
        )

    counts = eval_df[conc_col].value_counts(dropna=False)
    total = len(eval_df)
    out = pd.DataFrame(
        {
            "concordance": counts.index.astype(str),
            "count": counts.values,
            "percentage": (counts.values / total * 100.0).round(2),
        }
    )
    return out.sort_values("count", ascending=False).reset_index(drop=True)


def main() -> None:
    repo = _repo_root()
    out_dir = repo / "analyses" / "generated_tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    t1 = build_table1_dataset_summary(repo)
    t2 = build_table2_feature_summary(repo)
    t3 = build_table3_concordance(repo)

    t1.to_csv(out_dir / "table1_dataset_summary.csv", index=False)
    t2.to_csv(out_dir / "table2_feature_summary.csv", index=False)
    t3.to_csv(out_dir / "table3_concordance_summary.csv", index=False)

    print("Wrote:")
    print(f"- {(out_dir / 'table1_dataset_summary.csv').relative_to(repo)}")
    print(f"- {(out_dir / 'table2_feature_summary.csv').relative_to(repo)}")
    print(f"- {(out_dir / 'table3_concordance_summary.csv').relative_to(repo)}")


if __name__ == "__main__":
    main()
