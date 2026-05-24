#!/usr/bin/env python3
"""Stage 4a: threaded Gemini evaluation with checkpointing."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import DEFAULT_LLM_MODEL, PATHS
from clinvar.llm.evaluate_parallel import run_evaluation_pipeline
from clinvar.prompts import build_prompts_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant-type", choices=["snp", "indel"], required=True)
    parser.add_argument("--model", default=DEFAULT_LLM_MODEL)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--sample-size", type=int, default=None)
    parser.add_argument("--checkpoint", default=None, help="Checkpoint parquet path")
    parser.add_argument(
        "--skip-prompt-build",
        action="store_true",
        help="Use existing *_with_prompts parquet",
    )
    args = parser.parse_args()

    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        raise SystemExit("Set GEMINI_API_KEY environment variable")

    if args.variant_type == "snp":
        prompted_path = PATHS["snp_with_prompts"]
        output_path = PATHS["snp_evaluation_results"]
    else:
        prompted_path = PATHS["indel_with_prompts"]
        output_path = PATHS["indel_evaluation_results"]

    if args.skip_prompt_build and Path(prompted_path).exists():
        df = pq.read_table(prompted_path).to_pandas()
    else:
        signaled = (
            PATHS["snp_annotated_signaled"]
            if args.variant_type == "snp"
            else PATHS["indel_annotated_signaled"]
        )
        df = pq.read_table(signaled).to_pandas()
        df = build_prompts_table(df, variant_type=args.variant_type)

    checkpoint = args.checkpoint or output_path.replace(".parquet", "_checkpoint.parquet")

    run_evaluation_pipeline(
        df_with_prompts=df,
        variant_type=args.variant_type,
        output_path=output_path,
        checkpoint_path=checkpoint,
        model_name=args.model,
        api_key=api_key,
        max_workers=args.workers,
        sample_size=args.sample_size,
    )


if __name__ == "__main__":
    main()
