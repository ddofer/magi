"""Behaviour of the gradio_app modules that does not need a model."""

from contextlib import redirect_stdout
import importlib.util
from io import StringIO
import math
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
APP_DIR = ROOT / "gradio_app"


def load(name):
    spec = importlib.util.spec_from_file_location(f"app_{name}", APP_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


inference = load("inference")
analysis = load("analysis")
tracks = load("tracks")

GENOME = {"chr1": "ACGT" * 25}  # 100 bp; pyfaidx-like: `in` and [chrom][start:end]


class SequenceTest(unittest.TestCase):
    def test_ref_mismatch_is_an_error_not_a_warning(self):
        with redirect_stdout(StringIO()), self.assertRaisesRegex(ValueError, "REF mismatch"):
            inference.get_genomic_sequence(GENOME, "chr1", 50, "G", "T", context_len=20)

    def test_matching_ref_builds_alt_window(self):
        # 1-based 50 is 0-based 49 -> "C" in ACGT repeats
        ref_seq, alt_seq, center = inference.get_genomic_sequence(GENOME, "chr1", 50, "C", "T", context_len=20)
        self.assertEqual((ref_seq[center], alt_seq[center], center), ("C", "T", 10))

    def test_empty_slice_from_truncated_contig_is_none(self):
        self.assertIsNone(inference._fetch_sequence_from_local({"chr15": "ACGT"}, "chr15", 80_000_000, 80_000_100))


class BaselineTest(unittest.TestCase):
    def test_indel_uses_pooled_benign_background(self):
        analysis._MAGI_BASELINE_CACHE = None
        baseline = analysis.load_magi_baseline()
        # pooled benign indel moments from clinvar_indel_deltas.parquet (first block of the CSV)
        row = pd.Series({"ref": "AG", "alt": "TC", "D_BED_protein_coding_gene": -0.000676 + 0.020592})
        z_sum = analysis.compute_global_z_sum_log(row, baseline, variant_type="indel")
        self.assertAlmostEqual(z_sum, math.log(2), places=3)


def fake_profiles(chrom, pos, ref, alt):
    length = 40
    return {
        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
        "variant_center": 20, "input_len": length,
        "bed_names": ["exon"], "bigwig_names": [],
        "bed_ref": np.full((length, 1), 0.2), "bed_alt": np.full((length, 1), 0.6),
        "bed_track_len": length, "bed_track_start": 0,
        "bw_ref": None, "bw_alt": None,
    }


RANKED = [{"track_id": "exon", "track_type": "BED", "display_name": "exon", "delta": 0.4}]


class RegionPlotTest(unittest.TestCase):
    def test_plot_uses_given_profiles_not_module_global(self):
        inference._LAST_TRACK_PROFILES = fake_profiles("chr17", 43092418, "T", "C")
        mine = fake_profiles("chr7", 117559590, "ATCT", "A")
        fig = tracks.generate_region_tracks_plot(ranked_tracks=RANKED, visible_radius_bp=10, profiles=mine)
        self.assertIn("chr7:117,559,590 ATCT>A", fig._suptitle.get_text())
        self.assertEqual(tracks.get_track_view_bounds(mine)["max_radius"], 19)


if __name__ == "__main__":
    unittest.main()
