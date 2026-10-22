"""A frozen benign background must reproduce the recomputed z-scores and work without labels."""

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "clinvar_pipeline"))
from clinvar.impact_scoring import compute_impact_frame  # noqa: E402
from tools.build_background import build_background  # noqa: E402

SCORES = ["BW_z_mean_abs", "BW_z_max_abs", "BW_z_sum_log", "BW_z_top10_mean",
          "Global_z_mean_abs", "Global_z_max_abs", "Global_z_sum_log"]


def toy_cohort():
    rng = np.random.default_rng(0)
    types = np.repeat(["SNP", "Deletion", "Insertion", "Indel"], 30)
    df = pd.DataFrame(rng.normal(size=(len(types), 5)),
                      columns=["D_BED_exon", "D_BED_ORF", "D_BW_a", "D_BW_b", "D_BW_c"])
    df.iloc[::4] *= 5  # per-type benign moments differ from the pooled ones
    df["variant_type"] = types
    df["label"] = np.tile([True, False], len(types) // 2)
    return df


class FrozenBackgroundTest(unittest.TestCase):
    def test_own_background_matches_recomputed(self):
        df = toy_cohort()
        background = build_background(df[~df["label"]])
        self.assertEqual(list(background["variant_type"].unique()), ["snp", "indel", "deletion", "insertion"])
        recomputed = compute_impact_frame(df)
        frozen = compute_impact_frame(df, background=background)
        np.testing.assert_allclose(frozen[SCORES], recomputed[SCORES], rtol=1e-12)

    def test_unlabelled_rows_get_nonzero_scores(self):
        df = toy_cohort()
        background = build_background(df[~df["label"]])
        unlabelled = compute_impact_frame(df.drop(columns="label"), background=background)
        self.assertTrue((unlabelled["Global_z_sum_log"] > 0).all())


if __name__ == "__main__":
    unittest.main()
