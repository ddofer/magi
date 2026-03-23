#!/usr/bin/env python3
"""
Smoke test: verify center-focused inference with smaller config.

Usage:
  python smoke_test.py --mode snp|indel [--model-size 100M|650M]
"""

import sys
import argparse
import numpy as np
import nt3_inference as base


def run_smoke_test(mode="snp", model_size="100M"):
    """Configure and run inference with reduced window + sample count."""

    # Override config for faster iteration
    base.DEBUG_MODE = False
    base.CONTEXT_LEN = 8 * 1024  # 8K window instead of 64K
    base.BATCH_SIZE = 2
    base.USE_BED = True
    base.USE_BIGWIGS = True
    base.USE_KL_DIVERGENCE = True
    base.USE_EMBEDDINGS = False  # Disable for speed

    # Use smaller model
    model_key = "100M" if model_size == "100M" else "650M"
    base.MODEL_NAME = f"InstaDeepAI/NTv3_{model_key}_post"

    if mode == "snp":
        print(f"🧪 SNP Smoke Test (8K window, {model_size} model)")
        base.INPUT_DATA_FILE = "clinvar_input.parquet"
        base.OUTPUT_RESULTS_FILE = "smoke_test_snp_deltas.parquet"
    else:  # indel
        print(f"🧪 INDEL Smoke Test (8K window, {model_size} model)")
        base.INPUT_DATA_FILE = "clinvar_indel.parquet"
        base.OUTPUT_RESULTS_FILE = "smoke_test_indel_deltas.parquet"
        base.USE_EMBEDDINGS = True  # INDELs benefit from embeddings

    try:
        # Load data and limit to 150 samples
        import pandas as pd

        df = pd.read_parquet(base.INPUT_DATA_FILE)
        df_small = df.head(150).copy()

        print(f"   Loaded {len(df_small)} samples from {base.INPUT_DATA_FILE}")
        print(f"   Label distribution: {dict(df_small['label'].value_counts())}")

        # Save temp input
        temp_input = f"smoke_test_{mode}_input.parquet"
        df_small.to_parquet(temp_input, index=False)
        base.INPUT_DATA_FILE = temp_input

        # Run
        print("   Starting inference...")
        rdf = base.main()

        # Quick validation
        if rdf is not None:
            print("\n✅ Smoke test completed!")
            print(f"   Output shape: {rdf.shape}")

            # Check delta ranges
            delta_bed_cols = [c for c in rdf.columns if c.startswith("D_BED_")]
            delta_bw_cols = [c for c in rdf.columns if c.startswith("D_BW_")]

            if delta_bed_cols:
                bed_deltas = rdf[delta_bed_cols].values.flatten()
                bed_deltas = bed_deltas[~np.isnan(bed_deltas)]
                print(f"\n   BED Deltas ({len(bed_deltas)} values):")
                print(f"      Range: [{bed_deltas.min():.6f}, {bed_deltas.max():.6f}]")
                print(f"      Mean |Δ|: {np.abs(bed_deltas).mean():.6f}")
                print(f"      Median |Δ|: {np.median(np.abs(bed_deltas)):.6f}")

                # Check for start_codon
                if "D_BED_start_codon" in rdf.columns:
                    sc_vals = rdf["D_BED_start_codon"].dropna()
                    print(f"\n   START_CODON Deltas (n={len(sc_vals)}):")
                    print(f"      Range: [{sc_vals.min():.6f}, {sc_vals.max():.6f}]")
                    print(f"      Mean: {sc_vals.mean():.6f}")
                    print(f"      StdDev: {sc_vals.std():.6f}")

            if delta_bw_cols:
                bw_deltas = rdf[delta_bw_cols].values.flatten()
                bw_deltas = bw_deltas[~np.isnan(bw_deltas)]
                print(f"\n   BigWig Deltas ({len(bw_deltas)} values):")
                print(f"      Range: [{bw_deltas.min():.6f}, {bw_deltas.max():.6f}]")
                print(f"      Mean |Δ|: {np.abs(bw_deltas).mean():.6f}")
                print(f"      Median |Δ|: {np.median(np.abs(bw_deltas)):.6f}")

            # Check REF bounds (should be probabilities in [0, 1])
            ref_bed_cols = [c for c in rdf.columns if c.startswith("REF_BED_")]
            ref_bw_cols = [c for c in rdf.columns if c.startswith("REF_BW_")]

            if ref_bed_cols:
                ref_bed_vals = rdf[ref_bed_cols].values.flatten()
                ref_bed_vals = ref_bed_vals[~np.isnan(ref_bed_vals)]
                in_bounds = ((ref_bed_vals >= 0) & (ref_bed_vals <= 1)).sum()
                print(f"\n   REF_BED Bounds: {in_bounds}/{len(ref_bed_vals)} in [0, 1]")
                if in_bounds != len(ref_bed_vals):
                    print("      ⚠️  Out-of-bounds values detected!")
                    print(
                        f"      Range: [{ref_bed_vals.min():.6f}, {ref_bed_vals.max():.6f}]"
                    )

            if ref_bw_cols:
                ref_bw_vals = rdf[ref_bw_cols].values.flatten()
                ref_bw_vals = ref_bw_vals[~np.isnan(ref_bw_vals)]
                in_bounds = ((ref_bw_vals >= 0) & (ref_bw_vals <= 1)).sum()
                print(f"   REF_BW Bounds: {in_bounds}/{len(ref_bw_vals)} in [0, 1]")
                if in_bounds != len(ref_bw_vals):
                    print("      ⚠️  Out-of-bounds values detected!")
                    print(
                        f"      Range: [{ref_bw_vals.min():.6f}, {ref_bw_vals.max():.6f}]"
                    )

            print(f"\n   Output: {base.OUTPUT_RESULTS_FILE}")
            return True
        else:
            print("❌ Inference failed")
            return False
    except Exception as e:
        print(f"❌ Smoke test failed: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Smoke test center-focused inference")
    parser.add_argument(
        "--mode",
        choices=["snp", "indel"],
        default="snp",
        help="SNP or INDEL inference mode",
    )
    parser.add_argument(
        "--model-size",
        choices=["100M", "650M"],
        default="100M",
        help="Model size to use",
    )
    args = parser.parse_args()

    success = run_smoke_test(mode=args.mode, model_size=args.model_size)
    sys.exit(0 if success else 1)
