#!/usr/bin/env python3
"""Stage 4b: Gemini batch API evaluation."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import DEFAULT_LLM_MODEL, PATHS
from clinvar.llm.evaluate_batch import run_batch_evaluation_pipeline
from clinvar.prompts import build_prompts_table, get_system_prompt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant-type", choices=["snp", "indel"], required=True)
    parser.add_argument("--model", default=DEFAULT_LLM_MODEL)
    parser.add_argument("--batch-name", default=None)
    parser.add_argument("--submit-only", action="store_true")
    parser.add_argument("--sample-size", type=int, default=None)
    args = parser.parse_args()

    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        raise SystemExit("Set GEMINI_API_KEY environment variable")

    signaled = (
        PATHS["snp_annotated_signaled"]
        if args.variant_type == "snp"
        else PATHS["indel_annotated_signaled"]
    )
    df = pq.read_table(signaled).to_pandas()
    if args.sample_size:
        df = df.sample(n=args.sample_size, random_state=42)

    df = build_prompts_table(df, variant_type=args.variant_type)
    system_prompt = get_system_prompt(args.variant_type)
    batch_name = args.batch_name or f"{args.variant_type}_concordance_eval"
    work_dir = str(Path(PATHS["batch_work_dir"]) / f"{args.variant_type}s")

    run_batch_evaluation_pipeline(
        df_with_prompts=df,
        system_prompt=system_prompt,
        variant_type=args.variant_type,
        batch_name=batch_name,
        work_dir=work_dir,
        api_key=api_key,
        model=args.model,
        wait_for_completion=not args.submit_only,
    )


if __name__ == "__main__":
    main()
