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
    parser.add_argument("--background", help="frozen benign background CSV")
    args = parser.parse_args()

    print("Loading delta parquets (track columns only) …")
    df = prepare_combined_deltas(
        PATHS["snp_deltas"],
        PATHS["indel_deltas"],
        [PATHS["af_snps"], PATHS["af_indels"]],
    )
    print(f"  Combined: {len(df):,} variants")

    print(f"Computing impact scores (z-score background: {args.background or 'Benign, per variant_type'}) …")
    df = compute_impact_frame(df, background=args.background)

    snp_impact = slim_impact_columns(df[df["variant_type"] == "SNP"])
    indel_impact = slim_impact_columns(df[df["variant_type"] != "SNP"])

    # A frozen-background run gets its own files, so the paper figures keep reading the default ones.
    for impact, key in ((snp_impact, "snp_deltas_impact"), (indel_impact, "indel_deltas_impact")):
        path = Path(PATHS[key])
        if args.background:
            path = path.with_name(f"{path.stem}_frozen_bg{path.suffix}")
        path.parent.mkdir(parents=True, exist_ok=True)
        impact.to_parquet(path, index=False)
        print(f"  Saved → {path}")
    print("\n✓ Impact scoring complete.")


if __name__ == "__main__":
    main()
