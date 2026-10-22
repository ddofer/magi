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

RULE_TABLE = OUTPUT_ROOT / "mechanisms" / "rule_table.csv"


def write_rule_table(out=RULE_TABLE) -> None:
    table = pd.concat([rule_table("snp"), rule_table("indel")], ignore_index=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False)
    print(f"Saved {len(table)} conditions to {out}")
    blocked = table.loc[~table["can_be_primary"], ["variant_type", "mechanism"]].drop_duplicates()
    if len(blocked):
        print("Rules that can never be the assigned mechanism:")
        print(blocked.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(RULE_TABLE))
    write_rule_table(parser.parse_args().out)


if __name__ == "__main__":
    main()
