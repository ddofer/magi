#!/usr/bin/env python3
"""Write the rule layer (SNV and indel rules) as one CSV for the supplement."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "vendor"))

from config import OUTPUT_ROOT
from nt_mechanism_assignment import rule_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(OUTPUT_ROOT / "mechanisms" / "rule_table.csv"))
    args = parser.parse_args()

    table = pd.concat([rule_table("snp"), rule_table("indel")], ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, index=False)
    print(f"Saved {len(table)} conditions to {args.out}")
    blocked = table.loc[~table["can_be_primary"], ["variant_type", "mechanism"]].drop_duplicates()
    if len(blocked):
        print("Rules that can never be the assigned mechanism:")
        print(blocked.to_string(index=False))


if __name__ == "__main__":
    main()
