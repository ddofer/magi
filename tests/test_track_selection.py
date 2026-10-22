"""Track filtering must preserve model order and use the first metadata match."""

from contextlib import redirect_stdout
import importlib.util
from io import StringIO
from pathlib import Path
import unittest
from unittest.mock import patch

import pandas as pd

import inference


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("app_inference", ROOT / "gradio_app/inference.py")
app_inference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(app_inference)
SELECTORS = (inference.get_track_indices, app_inference.get_track_indices)


class TrackSelectionTest(unittest.TestCase):
    def test_filter_and_order_with_duplicate_metadata(self):
        metadata = pd.DataFrame([
            ("histone", "ChIP-seq", "H3K27ac", "encode_v3"),
            ("other_chip", "ChIP-seq", "CTCF", "encode_v3"),
            ("accessibility", "ATAC-seq", "", "encode_v3"),
            ("wrong_source", "ATAC-seq", "", "geo"),
            ("cage", "CAGE", "", "other"),
            ("fantom", "", "", "fantom5"),
            ("expression", "", "", "gtex"),
            ("duplicate", "ChIP-seq", "CTCF", "geo"),
            ("duplicate", "ChIP-seq", "H3K4me3", "geo"),
        ], columns=["file_id", "assay", "experiment_target", "dataset"])
        names = ["expression", "unknown", "other_chip", "histone", "duplicate",
                 "accessibility", "wrong_source", "cage", "fantom", "histone"]
        expected_indices = [0, 3, 5, 7, 8, 9]
        expected = (expected_indices, [names[i] for i in expected_indices])
        for select in SELECTORS:
            with (
                self.subTest(module=select.__module__),
                patch.object(pd, "read_csv", return_value=metadata),
                redirect_stdout(StringIO()),
            ):
                self.assertEqual(select(names, Path(__file__)), expected)

    def test_missing_metadata_keeps_all_tracks(self):
        names = ["track_b", "track_a"]
        for select in SELECTORS:
            with self.subTest(module=select.__module__), redirect_stdout(StringIO()):
                self.assertEqual(
                    select(names, ROOT / "scratch/nonexistent.csv"), ([0, 1], names)
                )


if __name__ == "__main__":
    unittest.main()
