#!/usr/bin/env python3
"""
NTv3 INDEL Inference (thin wrapper around nt3_inference.py)

Purpose
  Configure the shared pipeline for INDEL inputs and run the same inference
  stack as the SNP pipeline, with embeddings enabled by default.

Usage
  python nt3_inference_indel.py

Notes
  - This wrapper only sets config values; logic lives in nt3_inference.py
  - To run a small test, uncomment DEBUG_MODE and CONTEXT_LEN below
  - Output file is separate from SNP outputs to avoid overwrites
"""

import nt3_inference as base

# ── Config overrides ────────────────────────────────────────────────────────
base.INPUT_DATA_FILE = "clinvar_indel.parquet"
base.OUTPUT_RESULTS_FILE = "clinvar_indel_deltas.parquet"
base.USE_KL_DIVERGENCE = True
base.USE_EMBEDDINGS = True  # INDEL benefits from embedding distances
# base.DEBUG_MODE        = True   # uncomment for quick test (33 variants)
# base.CONTEXT_LEN       = 1 * 1024

if __name__ == "__main__":
    base.main()
