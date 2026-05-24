#!/usr/bin/env python3
"""Stage 3: build llm_prompt column on signaled parquets."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import PATHS
from clinvar.prompts import build_prompts_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant-type",
        choices=["snp", "indel"],
        required=True,
    )
    args = parser.parse_args()

    if args.variant_type == "snp":
        input_path = PATHS["snp_annotated_signaled"]
        output_path = PATHS["snp_with_prompts"]
    else:
        input_path = PATHS["indel_annotated_signaled"]
        output_path = PATHS["indel_with_prompts"]

    df = pq.read_table(input_path).to_pandas()
    print(f"Loaded {len(df)} variants from {input_path}")

    out = build_prompts_table(df, variant_type=args.variant_type)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(output_path, index=False)
    print(f"Saved {len(out)} prompts → {output_path}")


if __name__ == "__main__":
    main()
