"""Check smoke-test failure reporting without loading a model."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import os
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

import smoke_test


ROOT = Path(__file__).resolve().parents[1]


class SmokeTest(unittest.TestCase):
    def run_smoke(self, result, **kwargs):
        scratch = ROOT / "scratch"
        scratch.mkdir(exist_ok=True)
        original_cwd = Path.cwd()
        original_config = vars(smoke_test.base).copy()
        try:
            with tempfile.TemporaryDirectory(dir=scratch) as work:
                os.chdir(work)
                input_df = pd.DataFrame({"label": ["Benign"] * 200})

                def fake_main():
                    self.use_embeddings = smoke_test.base.USE_EMBEDDINGS
                    return result

                with (
                    patch.object(pd, "read_parquet", return_value=input_df),
                    patch.object(pd.DataFrame, "to_parquet"),
                    patch.object(smoke_test.base, "main", fake_main),
                    redirect_stdout(StringIO()),
                ):
                    return smoke_test.run_smoke_test(**kwargs)
        finally:
            os.chdir(original_cwd)
            for name, value in original_config.items():
                setattr(smoke_test.base, name, value)

    def test_valid_probabilities(self):
        result = pd.DataFrame({
            "D_BED_exon": [0.1], "REF_BED_exon": [0.4],
            "D_BW_track": [-0.2], "REF_BW_track": [0.7],
        })
        self.assertTrue(self.run_smoke(result))

    def test_out_of_bounds_probabilities_fail(self):
        result = pd.DataFrame({
            "D_BED_exon": [0.1], "REF_BED_exon": [1.4],
            "D_BW_track": [-0.2], "REF_BW_track": [0.7],
        })
        self.assertFalse(self.run_smoke(result))

    def test_missing_track_outputs_fail(self):
        self.assertFalse(self.run_smoke(pd.DataFrame({"LLR": [0.1]})))

    def test_partial_nan_track_output_passes(self):
        result = pd.DataFrame({
            "D_BED_exon": [float("nan"), 0.1], "REF_BED_exon": [float("nan"), 0.4],
            "D_BW_track": [-0.2, 0.0], "REF_BW_track": [0.7, 0.2],
        })
        self.assertTrue(self.run_smoke(result))

    def test_all_zero_delta_column_fails(self):
        # Every delta zero means ref and alt windows were identical (e.g. all-N sequence).
        result = pd.DataFrame({
            "D_BED_exon": [0.0, 0.0], "REF_BED_exon": [0.4, 0.3],
            "D_BW_track": [-0.2, 0.1], "REF_BW_track": [0.7, 0.2],
        })
        self.assertFalse(self.run_smoke(result))

    def test_infinite_output_fails(self):
        result = pd.DataFrame({
            "D_BED_exon": [float("inf")], "REF_BED_exon": [0.4],
            "D_BW_track": [-0.2], "REF_BW_track": [0.7],
        })
        self.assertFalse(self.run_smoke(result))

    def test_indel_mode_enables_embeddings(self):
        result = pd.DataFrame({
            "D_BED_exon": [0.1], "REF_BED_exon": [0.4],
            "D_BW_track": [-0.2], "REF_BW_track": [0.7],
        })
        self.run_smoke(result, mode="indel")
        self.assertTrue(self.use_embeddings)
        self.run_smoke(result, mode="snp")
        self.assertFalse(self.use_embeddings)

    def test_empty_output_fails(self):
        self.assertFalse(self.run_smoke(pd.DataFrame()))

    def test_all_nan_track_output_fails(self):
        result = pd.DataFrame({
            "D_BED_exon": [float("nan")], "REF_BED_exon": [0.4],
            "D_BW_track": [-0.2], "REF_BW_track": [0.7],
        })
        self.assertFalse(self.run_smoke(result))


if __name__ == "__main__":
    unittest.main()
