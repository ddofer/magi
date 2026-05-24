#!/usr/bin/env python3
"""Compute global z-score impact columns for the full delta cohort (SNP + indel)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import PATHS
from clinvar.impact_scoring import compute_impact_frame, prepare_combined_deltas, slim_impact_columns


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    print("Loading delta parquets (track columns only) …")
    df = prepare_combined_deltas(
        PATHS["snp_deltas"],
        PATHS["indel_deltas"],
        [PATHS["af_snps"], PATHS["af_indels"]],
    )
    print(f"  Combined: {len(df):,} variants")

    print("Computing impact scores (z-score background: Benign, per variant_type) …")
    df = compute_impact_frame(df)

    snp_impact = slim_impact_columns(df[df["variant_type"] == "SNP"])
    indel_impact = slim_impact_columns(df[df["variant_type"] != "SNP"])

    for path in (PATHS["snp_deltas_impact"], PATHS["indel_deltas_impact"]):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    snp_impact.to_parquet(PATHS["snp_deltas_impact"], index=False)
    indel_impact.to_parquet(PATHS["indel_deltas_impact"], index=False)
    print(f"  Saved → {PATHS['snp_deltas_impact']}")
    print(f"  Saved → {PATHS['indel_deltas_impact']}")
    print("\n✓ Impact scoring complete.")


if __name__ == "__main__":
    main()
