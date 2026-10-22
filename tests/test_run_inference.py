"""run_inference on a tiny fake NTv3 (CPU, no model load): order, failures, shared windows, resume."""

import tempfile
import types
import unittest
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import inference  # noqa: E402

C, CTX = 8, 1024
VOCAB = {"<unk>": 0, "<pad>": 1, "<mask>": 2, "<cls>": 3, "<eos>": 4, "<bos>": 5, "A": 6, "T": 7, "C": 8, "G": 9, "N": 10}


class Head(nn.Module):
    def __init__(self, n):
        super().__init__()
        self.head = nn.Linear(C, n)

    def forward(self, x):
        return F.softplus(self.head(x))


class Bed(nn.Module):
    def __init__(self):
        super().__init__()
        self.head = nn.Linear(C, 4)

    def forward(self, x):
        return self.head(x).reshape(*x.shape[:-1], 2, 2)


class FakeNTv3(nn.Module):
    """Same interface as NTv3PostTrained; a local conv, so letters far from an edit do not change."""

    def __init__(self):
        super().__init__()
        torch.manual_seed(0)
        self.emb, self.conv, self.lm = nn.Embedding(11, C), nn.Conv1d(C, C, 5, padding=2), nn.Linear(C, 11)
        self.core = types.SimpleNamespace(bed_head=Bed(), bigwig_head=types.SimpleNamespace(species_heads=[Head(3)]))
        self.config = types.SimpleNamespace(
            keep_target_center_fraction=0.375, num_species_special_tokens=6, species_to_token_id={"human": 6}
        )
        self.calls = 0

    def encode_species(self, names):
        return torch.tensor([self.config.species_to_token_id[n] for n in names])

    def forward(self, input_ids, species_ids):
        self.calls += len(input_ids)
        h = torch.tanh(self.conv(self.emb(input_ids).transpose(1, 2)).transpose(1, 2))
        return types.SimpleNamespace(logits=self.lm(h), embedding=h)


class Tok:
    pad_token_id = 1

    def get_vocab(self):
        return VOCAB

    def __call__(self, s, add_special_tokens=False):
        return {"input_ids": [VOCAB[ch] for ch in s]}


class Rec(str):
    def __getitem__(self, sl):
        return types.SimpleNamespace(seq=str.__getitem__(self, sl))


rng = np.random.default_rng(0)
GENOME = {"chr1": Rec("".join(rng.choice(list("ACGT"), 8000)))}


def snv(pos):
    ref = GENOME["chr1"][pos - 1 : pos].seq
    return ("chr1", pos, ref, "A" if ref != "A" else "C")


def run(df, **kw):
    model = FakeNTv3()
    out = inference.run_inference(
        df, model, Tok(), GENOME, ["e0", "e1"], ["t0", "t1", "t2"], [0, 2], "cpu",
        context_len=CTX, batch_size=1, use_embeddings=True, mlm_window=4, **kw,
    )
    return out, model.calls


class RunInferenceTest(unittest.TestCase):
    def setUp(self):
        rows = [snv(4000), snv(4128), ("chr9", 50, "A", "G"), snv(4005), snv(2000)]
        self.df = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt"])

    def test_order_failures_and_columns(self):
        out, _ = run(self.df)
        self.assertEqual(out["pos"].tolist(), self.df["pos"].tolist())
        self.assertTrue(out.iloc[2][["D_BED_e0", "LLR", "EMB_l2_dist"]].isna().all())  # missing contig
        self.assertFalse(out.drop(index=2)[["D_BED_e0", "D_BW_t2", "LLR", "EMB_l2_dist"]].isna().any().any())
        self.assertEqual([c for c in out.columns if c.startswith("D_BW_")], ["D_BW_t0", "D_BW_t2"])
        self.assertNotIn("window_offset", out.columns)

    def test_shared_windows(self):
        centred, n_centred = run(self.df)
        shared, n_shared = run(self.df, share_windows=True)
        self.assertEqual(n_centred, 8)  # 4 scorable variants x (ref + alt)
        self.assertEqual(n_shared, 7)  # 4000 and 4128 (same place in a 128-letter token) share a ref
        self.assertEqual(shared["window_offset"].tolist()[:2], [0, 128])
        # The fake model is local, so moving the window does not change the features.
        cols = [c for c in centred.columns if c.startswith(("D_", "REF_BW", "MLM_", "LLR", "EMB_"))]
        np.testing.assert_allclose(shared[cols].to_numpy(float), centred[cols].to_numpy(float), atol=1e-5)

    def test_shared_window_at_contig_edges_matches_centred(self):
        df = pd.DataFrame([snv(100), snv(7950)], columns=["chrom", "pos", "ref", "alt"])
        centred, _ = run(df)
        shared, _ = run(df, share_windows=True)
        self.assertEqual(shared["window_offset"].tolist(), [-413, 0])  # start clipped at 0; short at the end
        np.testing.assert_allclose(shared["MLM_KL_mean"], centred["MLM_KL_mean"], atol=1e-6)

    def test_shared_groups_split_by_species(self):
        rows = [dict(chrom="chr1", pos=4000, ref="A", _sid=6), dict(chrom="chr1", pos=100, ref="A", _sid=7)]
        self.assertEqual([ks for *_, ks in inference._shared_windows(rows, CTX, 192)], [[0], [1]])

    def test_shared_deletion_keeps_full_length_alt(self):
        pos = 4000
        ref = GENOME["chr1"][pos - 1 : pos + 3].seq
        rows = list(inference._variant_pairs([{"chrom": "chr1", "pos": pos, "ref": ref, "alt": ref[0]}],
                                             GENOME, CTX, True, 192))
        (ref_seq, [(_, alt_seq, vc)]), = rows
        self.assertEqual(len(alt_seq), CTX)
        self.assertEqual(alt_seq[: vc + 1] + GENOME["chr1"][pos + 3 : pos + 3 + 10].seq, alt_seq[: vc + 11])

    def test_parts_resume(self):
        with tempfile.TemporaryDirectory() as d:
            first, _ = run(self.df, parts_dir=d, chunk_rows=2)
            self.assertEqual(len(list(Path(d).glob("part-*.parquet"))), 3)
            again, calls = run(self.df, parts_dir=d, chunk_rows=2)
            self.assertEqual(calls, 0)
            pd.testing.assert_frame_equal(first, again)
            with self.assertRaises(ValueError):  # other settings never reuse these parts
                run(self.df, parts_dir=d, chunk_rows=2, share_windows=True)

    def test_token_lut_matches_tokenizer(self):
        lut = inference.build_token_lut(Tok())
        ids = inference.tokenize(["ACGTNacgtn-"], lut, 16, 1)[0].tolist()
        self.assertEqual(ids, [6, 8, 9, 7, 10, 6, 8, 9, 7, 10, 0, 1, 1, 1, 1, 1])


if __name__ == "__main__":
    unittest.main()
