"""Supplementary Fig S1 — per-species OMIA concordance (SNP vs indel) → figs1.png."""

from __future__ import annotations

from pathlib import Path

from config import PATHS
from clinvar.animals_data import load_animals_eval_res, per_species_accuracy_table
from figures.common import ensure_dirs
from figures.plotting import plot_per_species_accuracy_snp_vs_indel


def run() -> None:
    ensure_dirs()
    res = load_animals_eval_res()
    summary = per_species_accuracy_table(res)

    out = Path(PATHS["figs1"])
    out.parent.mkdir(parents=True, exist_ok=True)
    plot_per_species_accuracy_snp_vs_indel(res, out)

    sidecar = out.with_suffix(".csv")
    summary.to_csv(sidecar, index=False)
    print(f"  FigS1 → {out}")
    print(f"  FigS1 summary → {sidecar}  (n={len(res):,})")
