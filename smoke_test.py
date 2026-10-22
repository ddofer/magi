#!/usr/bin/env python3
"""Run two variants through NTv3 and check track outputs. Files go in scratch/."""

import argparse
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd

import inference as base


def run_smoke_test(mode="snp", model_size="100M"):
    """Run a tiny inference check with an 8K window; embeddings on for indels, as in inference_indel.py."""
    work_dir = Path.cwd() / "scratch" / f"smoke_{mode}"
    work_dir.mkdir(parents=True, exist_ok=True)
    base.CONTEXT_LEN = 8 * 1024
    base.USE_BED = True
    base.USE_BIGWIGS = True
    base.USE_KL_DIVERGENCE = True
    base.USE_EMBEDDINGS = mode == "indel"
    base.MODEL_NAME = f"InstaDeepAI/NTv3_{model_size}_post"
    base.OUTPUT_RESULTS_FILE = str(work_dir / "deltas.parquet")
    input_file = "clinvar_input.parquet" if mode == "snp" else "clinvar_indel.parquet"

    try:
        df = pd.read_parquet(input_file).head(2)
        if df.empty:
            raise ValueError("Input contains no variants")
        base.INPUT_DATA_FILE = str(work_dir / "input.parquet")
        df.to_parquet(base.INPUT_DATA_FILE, index=False)
        print(f"Running {mode} smoke test: {len(df)} variants, {model_size}, 8K window")
        result = base.main()
        if result is None or result.empty:
            raise ValueError("Inference returned no variants")

        for prefix in ("D_BED_", "REF_BED_", "D_BW_", "REF_BW_", "LFC_", "WIN_"):
            columns = [c for c in result.columns if c.startswith(prefix)]
            if not columns:
                if prefix in ("LFC_", "WIN_"):  # absent from older outputs
                    continue
                raise ValueError(f"Missing {prefix} output columns")
            frame = result[columns].select_dtypes("number") if prefix == "LFC_" else result[columns]  # LFC_top_track
            values = frame.to_numpy().ravel()
            if np.isinf(values).any():
                raise ValueError(f"Infinite values in {prefix} outputs")
            values = values[~np.isnan(values)]  # NaN marks a variant outside the track window
            if values.size == 0:
                raise ValueError(f"All {prefix} outputs are NaN")
            if prefix.startswith("REF_") and not ((values >= 0) & (values <= 1)).all():
                raise ValueError(f"{prefix} probabilities outside [0, 1]")
            if prefix.startswith("D_") and not values.any():
                raise ValueError(f"All {prefix} deltas are zero: ref and alt windows identical")
            print(f"  {prefix}: range [{values.min():.6f}, {values.max():.6f}]")

        print(f"Smoke test passed: {base.OUTPUT_RESULTS_FILE}")
        return True
    except Exception as exc:
        traceback.print_exc()
        print(f"Smoke test failed: {exc}")
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["snp", "indel"], default="snp")
    parser.add_argument("--model-size", choices=["100M", "650M"], default="100M")
    args = parser.parse_args()
    sys.exit(0 if run_smoke_test(args.mode, args.model_size) else 1)
