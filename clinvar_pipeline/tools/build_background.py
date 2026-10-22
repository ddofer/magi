#!/usr/bin/env python3
"""Rebuild the frozen benign background (gradio_app/data/magi_baseline_stats.csv format).

For each variant type: mean and SD of every D_BED_/D_BW_ delta over ClinVar benign rows.
SD is the population SD (np.nanstd: ddof=0, NaN ignored), the same as the recomputed
background in clinvar.impact_scoring.z_score_matrix.
Blocks, in order: snp, indel (all indels pooled), deletion, insertion, indel (ClinVar
"Indel", i.e. delins). Rows follow the parquet column order. A parquet without a
variant_type column is read as SNVs.

    python tools/build_background.py SNV_DELTAS.parquet INDEL_DELTAS.parquet -o background.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from clinvar.impact_scoring import TRACK_PREFIXES, _harmonize_labels

INDEL_TYPES = ("Deletion", "Insertion", "Indel")


def read_benign(path: str) -> pd.DataFrame:
    pf = pq.ParquetFile(path)
    names = pf.schema.names
    keep = ["label", *(["variant_type"] if "variant_type" in names else [])]
    df = _harmonize_labels(pf.read(columns=keep + [c for c in names if c.startswith(TRACK_PREFIXES)]).to_pandas())
    if "variant_type" not in df.columns:
        df["variant_type"] = "SNP"
    return df[df["label"] == "Benign"]


def build_background(benign: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in benign.columns if c.startswith(TRACK_PREFIXES)]
    is_snp = benign["variant_type"] == "SNP"
    blocks = [("snp", is_snp), ("indel", ~is_snp)]
    blocks += [(t.lower(), benign["variant_type"] == t) for t in INDEL_TYPES]
    frames = []
    for name, mask in blocks:
        if not mask.any():
            continue
        m = benign.loc[mask, cols].to_numpy(np.float64)
        frames.append(pd.DataFrame({
            "variant_type": name,
            "delta_col": cols,
            "mean": np.nanmean(m, axis=0),
            "std": np.nanstd(m, axis=0),
        }))
        print(f"  {name}: {mask.sum():,} benign rows")
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("parquets", nargs="+", help="scored delta parquets (label, variant_type, D_BED_*, D_BW_*)")
    parser.add_argument("-o", "--out", required=True, help="output CSV")
    args = parser.parse_args()
    benign = pd.concat([read_benign(p) for p in args.parquets], ignore_index=True)
    build_background(benign).to_csv(args.out, index=False)
    print(f"  Saved → {args.out}")


if __name__ == "__main__":
    main()
