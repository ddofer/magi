"""Compute concordance summaries, including signal-stratified rates.

Usage:
  python analyses/extra_concordance_analysis.py
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _detect_column(df: pd.DataFrame, target: str) -> str | None:
    for c in df.columns:
        if c.lower() == target.lower():
            return c
    return None


def _summary_by_concordance(df: pd.DataFrame, conc_col: str) -> pd.DataFrame:
    counts = df[conc_col].value_counts(dropna=False)
    total = len(df)
    out = pd.DataFrame(
        {
            "concordance": counts.index.astype(str),
            "count": counts.values,
            "percentage": (counts.values / total * 100.0).round(2),
        }
    )
    return out.sort_values("count", ascending=False).reset_index(drop=True)


def _actionable_rate(df: pd.DataFrame, conc_col: str) -> pd.DataFrame:
    actionable = df[conc_col].astype(str).isin(["CONCORDANT", "PARTIAL"]).sum()
    discordant = df[conc_col].astype(str).eq("DISCORDANT").sum()
    total = len(df)
    return pd.DataFrame(
        [
            {
                "n_total": total,
                "n_actionable": actionable,
                "actionable_pct": round(actionable / total * 100.0, 2)
                if total
                else 0.0,
                "n_discordant": discordant,
                "discordant_pct": round(discordant / total * 100.0, 2)
                if total
                else 0.0,
            }
        ]
    )


def _signal_stratified(
    df: pd.DataFrame, conc_col: str, signal_col: str
) -> pd.DataFrame:
    rows: list[dict] = []
    for sig, group in df.groupby(signal_col, dropna=False):
        total = len(group)
        conc = group[conc_col].astype(str)
        n_conc = conc.eq("CONCORDANT").sum()
        n_part = conc.eq("PARTIAL").sum()
        n_disc = conc.eq("DISCORDANT").sum()
        n_na = conc.eq("NOT_APPLICABLE").sum()
        rows.append(
            {
                "signal_category": str(sig),
                "n_total": total,
                "n_concordant": n_conc,
                "concordant_pct": round(n_conc / total * 100.0, 2) if total else 0.0,
                "n_partial": n_part,
                "partial_pct": round(n_part / total * 100.0, 2) if total else 0.0,
                "n_actionable": n_conc + n_part,
                "actionable_pct": round((n_conc + n_part) / total * 100.0, 2)
                if total
                else 0.0,
                "n_discordant": n_disc,
                "discordant_pct": round(n_disc / total * 100.0, 2) if total else 0.0,
                "n_not_applicable": n_na,
            }
        )
    return pd.DataFrame(rows).sort_values("signal_category").reset_index(drop=True)


def analyze_file(path: Path, out_dir: Path) -> None:
    if not path.exists():
        print(f"Skip missing file: {path}")
        return

    df = pd.read_csv(path)
    conc_col = _detect_column(df, "concordance")
    if conc_col is None:
        print(f"No concordance column in {path}")
        return

    stem = path.stem
    conc_summary = _summary_by_concordance(df, conc_col)
    conc_summary.to_csv(out_dir / f"{stem}__concordance_distribution.csv", index=False)

    actionable = _actionable_rate(df, conc_col)
    actionable.to_csv(out_dir / f"{stem}__actionable_rate.csv", index=False)

    sig_col = _detect_column(df, "signal_category")
    if sig_col is not None:
        signal_summary = _signal_stratified(df, conc_col, sig_col)
        signal_summary.to_csv(
            out_dir / f"{stem}__signal_stratified_concordance.csv", index=False
        )

    print(f"Analyzed {path.name} -> outputs written")


def main() -> None:
    repo = _repo_root()
    out_dir = repo / "analyses" / "generated_tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = [
        repo / "analyses" / "stav-llm_evaluation_results.csv",
        repo / "analyses" / "pilot_llm_concordance.csv",
        repo / "stav_data" / "stav data" / "snp_LLM_results.csv",
        repo / "stav_data" / "stav data" / "indel_LLM_results.csv",
    ]

    for p in candidates:
        analyze_file(p, out_dir)

    print(f"Done. See outputs in {out_dir.relative_to(repo)}")


if __name__ == "__main__":
    main()
