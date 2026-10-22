#!/usr/bin/env python3
"""Run two variants through NTv3 and check track outputs. Files go in scratch/."""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

import inference as base


def run_smoke_test(mode="snp", model_size="100M"):
    """Run a tiny inference check with an 8K window and no embeddings."""
    work_dir = Path(__file__).resolve().parent / "scratch" / f"smoke_{mode}"
    work_dir.mkdir(parents=True, exist_ok=True)
    base.DEBUG_MODE = False
    base.CONTEXT_LEN = 8 * 1024
    base.BATCH_SIZE = 2
    base.USE_BED = True
    base.USE_BIGWIGS = True
    base.USE_KL_DIVERGENCE = True
    base.USE_EMBEDDINGS = False
    base.MODEL_NAME = f"InstaDeepAI/NTv3_{model_size}_post"
    base.OUTPUT_RESULTS_FILE = str(work_dir / "deltas.parquet")
    input_file = "clinvar_input.parquet" if mode == "snp" else "clinvar_indel.parquet"

    try:
        df = pd.read_parquet(input_file).head(2).copy()
        if df.empty:
            raise ValueError("Input contains no variants")
        base.INPUT_DATA_FILE = str(work_dir / "input.parquet")
        df.to_parquet(base.INPUT_DATA_FILE, index=False)
        print(f"Running {mode} smoke test: {len(df)} variants, {model_size}, 8K window")
        result = base.main()
        if result is None or result.empty:
            raise ValueError("Inference returned no variants")

        for kind in ("BED", "BW"):
            for prefix in (f"D_{kind}_", f"REF_{kind}_"):
                columns = [c for c in result.columns if c.startswith(prefix)]
                if not columns:
                    raise ValueError(f"Missing {prefix} output columns")
                values = result[columns].to_numpy().ravel()
                if not np.isfinite(values).all():
                    raise ValueError(f"Non-finite values in {prefix} outputs")
                if prefix.startswith("REF_") and not ((values >= 0) & (values <= 1)).all():
                    raise ValueError(f"{prefix} probabilities outside [0, 1]")
                print(f"  {prefix}: range [{values.min():.6f}, {values.max():.6f}]")

        print(f"Smoke test passed: {base.OUTPUT_RESULTS_FILE}")
        return True
    except Exception as exc:
        print(f"Smoke test failed: {exc}")
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["snp", "indel"], default="snp")
    parser.add_argument("--model-size", choices=["100M", "650M"], default="100M")
    args = parser.parse_args()
    sys.exit(0 if run_smoke_test(args.mode, args.model_size) else 1)
