#!/usr/bin/env python3
"""Stage 10: threaded Gemini concordance evaluation for OMIA animal variants."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import DEFAULT_LLM_MODEL, PATHS
from clinvar.animals_prompts import build_prompts_table, get_system_prompt
from clinvar.llm.evaluate_parallel import run_evaluation_pipeline


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
        help="Use existing animals_*_with_prompts parquet",
    )
    args = parser.parse_args()

    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        raise SystemExit("Set GEMINI_API_KEY environment variable")

    if args.variant_type == "snp":
        prompted_path = PATHS["animals_snp_with_prompts"]
        signaled_path = PATHS["animals_snp_signaled"]
        output_path = PATHS["animals_snp_eval"]
    else:
        prompted_path = PATHS["animals_indel_with_prompts"]
        signaled_path = PATHS["animals_indel_signaled"]
        output_path = PATHS["animals_indel_eval"]

    if args.skip_prompt_build and Path(prompted_path).exists():
        df = pq.read_table(prompted_path).to_pandas()
    else:
        df = pq.read_table(signaled_path).to_pandas()
        df = build_prompts_table(df, variant_type=args.variant_type)

    checkpoint = args.checkpoint or output_path.replace(".parquet", "_checkpoint.parquet")
    system_prompt = get_system_prompt(args.variant_type)

    run_evaluation_pipeline(
        df_with_prompts=df,
        variant_type=args.variant_type,
        output_path=output_path,
        checkpoint_path=checkpoint,
        model_name=args.model,
        api_key=api_key,
        max_workers=args.workers,
        sample_size=args.sample_size,
        system_prompt=system_prompt,
    )


if __name__ == "__main__":
    main()
