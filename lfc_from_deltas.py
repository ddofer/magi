#!/usr/bin/env python3
"""Add inference.py's LFC_ columns to an existing deltas parquet, from its REF_BW_ and D_BW_ columns.

MAGI stored sigmoid(raw signal), so the raw signal is logit(REF_BW) for Ref and logit(REF_BW + D_BW)
for Alt. This is exact for SNVs only: for longer variants the stored value is a mean of sigmoids over
the variant letters, not the sigmoid of a mean. Values stored as 1.0 (raw signal above ~17) cannot be
recovered; they become NaN and the LFC_ summaries skip them.

    python lfc_from_deltas.py clinvar_new_deltas.parquet clinvar_new_deltas_lfc.parquet

Reads the file once, converts a batch of tracks at a time, then writes the input plus the LFC_
columns. The 2 GB SNV file needs a few GB of memory: run it with srun.
"""
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

import inference

BATCH = 256  # tracks converted at a time


def logit(p):
    """Raw BigWig signal from a stored sigmoid value; NaN where it was stored as 1.0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.log(p) - np.log1p(-p)
    return np.where(np.isfinite(x), x, np.nan)


def main(src, dst):
    tab = pq.read_table(src)
    tracks = [c[len("REF_BW_") :] for c in tab.column_names if c.startswith("REF_BW_")]
    lfc = np.empty((tab.num_rows, len(tracks)), np.float32)
    for j in range(0, len(tracks), BATCH):
        t = tracks[j : j + BATCH]
        ref = np.column_stack([tab[f"REF_BW_{x}"].to_numpy() for x in t])
        d = np.column_stack([tab[f"D_BW_{x}"].to_numpy() for x in t])
        lfc[:, j : j + len(t)] = inference.log2_fold_change(logit(ref), logit(ref + d))
    cols = inference.lfc_summary(lfc, tracks, inference.track_groups(tracks, inference.METADATA_FILE))
    del lfc
    for name, v in cols.items():
        tab = tab.append_column(name, pa.array(v))
    pq.write_table(tab, dst)
    print(f"{dst}: {tab.num_rows} rows, +{len(cols)} LFC_ columns")


if __name__ == "__main__":
    main(*sys.argv[1:])
