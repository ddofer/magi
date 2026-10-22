#!/usr/bin/env python3
"""
NTv3 Genomic Variant Inference Pipeline (shared module)

Purpose
  Predict functional track changes for genomic variants using NTv3 and compare
  ref vs alt signals. Supports SNPs and INDELs using a shared pipeline.

Inputs
  - Parquet with required columns: chrom, pos (1-based), ref, alt
  - label and other extra columns are preserved in the output
  - Reference genome FASTA (hg38.fa) or UCSC fallback
  - Track metadata CSV for BigWig filtering (functional_tracks_metadata_human.csv)

Outputs
  - Parquet with ref/alt deltas for BED and BigWig tracks
  - MLM features:
      LLR, MLM_Prior, MLM_Delta (SNPs only)
      MLM_KL_mean, MLM_KL_max (KL(alt || ref) across window)
      MLM_logprob_ref/alt/delta (sequence log-prob around variant)
  - Embedding features (optional): EMB_* distances from the final embedding
  - LFC_*: log2 fold change of the raw BigWig signal (before MAGI's sigmoid) over the selected
    tracks: LFC_max_abs, LFC_mean_abs, LFC_top_track, LFC_max_abs_<assay group> (lfc_summary)
  - WIN_<element>_gain/loss(_pos): largest Alt-Ref rise (>= 0) / drop (<= 0) in the softmax
    probability of splice donor, splice acceptor, start and stop codon within WIN_LETTERS of the
    variant, and its offset in reference letters (NaN when there is no rise / drop)
  - Metadata: indel_size plus any extra input columns (+ window_offset with SHARE_WINDOWS)
  - Variants whose sequence cannot be fetched keep their row, with NaN features

Method notes
  - Ref/alt sequences are centered at pos and cropped to equal length
    (SHARE_WINDOWS: nearby variants share one window; see run_inference)
  - All features (BED, BigWig, MLM/KL) use variant-span extraction:
    single position for SNPs, mean over max(len(ref), len(alt)) for INDELs
  - Track heads cover the central 37.5% of the input at one letter per position;
    they run only on the letters read (forward_trunk)
  - KL features capture surprise/novelty in prediction distributions

Quick start
    - SNPs:   python inference.py
    - INDELs: python inference_indel.py
    - Any set: python inference.py --input x.parquet --output y.parquet [--compile] [--share-windows]
      (SLURM: sbatch inference.sbatch --input ...)

Examples
  - Override config in Python:
            import inference as base
      base.DEBUG_MODE = True
      base.OUTPUT_RESULTS_FILE = "tmp.parquet"
      base.main()
  - INDEL run with embeddings:
            import inference_indel  # sets base.USE_EMBEDDINGS = True
            import inference as base
      base.DEBUG_MODE = True
      base.main()

Config
  - CONTEXT_LEN, BATCH_SIZE, USE_BED, USE_BIGWIGS, SPECIES
  - USE_KL_DIVERGENCE, USE_EMBEDDINGS, MLM_WINDOW
  - SHARE_WINDOWS, COMPILE, CHUNK_ROWS

Dependencies
  torch, transformers, pandas, numpy, pyfaidx, tqdm, requests

Gotchas
  1) Species ID: must use model.encode_species(["human"])
  2) BED tracks can be missing: use getattr(output, "bed_tracks_logits", None)
  3) Use AutoModel (not AutoModelForMaskedLM) for NTv3
  4) BigWig/BED outputs are converted to probabilities via sigmoid
     BED tracks have shape (B, L', num_elements, 2): extract positive-class
     logit [..., 1] then apply sigmoid.  BigWig tracks have shape
     (B, L', num_tracks): apply sigmoid directly.
     This is MAGI's convention, kept so values stay comparable with the paper.
     The stock NTv3 pipeline differs: BigWig head output is already softplus
     (predicted signal, used as is), and BED uses softmax over the 2 classes.
     So REF_BW_* lies in [0.5, 1) and is not a probability.

Authors: Dan Ofer / Michal Linial Lab


Example/Note on NT V3 model example usage and outputs format (without our usage/additions):
https://huggingface.co/spaces/InstaDeepAI/ntv3

```
💻 Use a post-trained model
Here is a quick example of how to use the post-trained NTv3 650M model to predict tracks for a human genomic window.

from transformers import pipeline
import torch

model_name = "InstaDeepAI/NTv3_650M_pos"

ntv3_tracks = pipeline(
    "ntv3-tracks",
    model=model_name,
    trust_remote_code=True
)

# Run track prediction
out = ntv3_tracks(
  {
    "chrom": "chr19",
    "start": 6_700_000,
    "end": 6_831_072,
    "species": "human"
  }
)

# Print output shapes
# 7k human tracks over 37.5 % center region of the input sequence
print("bigwig_tracks_logits:", tuple(out.bigwig_tracks_logits.shape))
# Location of 21 genomic elements over 37.5 % center region of the input sequence
print("bed_tracks_logits:", tuple(out.bed_tracks_logits.shape))
# Language model logits for whole sequence over vocabulary
print("language model logits:", tuple(out.logits.shape))
```
"""

import hashlib
import json
import os
import shutil
import time
import warnings
import torch
import torch.nn.functional as F
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from pyfaidx import Fasta
import requests
from transformers import AutoModel, AutoTokenizer

# ============================================================================
# CONFIGURATION  (override any of these before calling main / run_inference)
# ============================================================================
DEBUG_MODE = False
CONTEXT_LEN = 64 * 1024
BATCH_SIZE = 2
USE_BED = True
USE_BIGWIGS = True
USE_KL_DIVERGENCE = True  # KL div + log-prob features  (fast, recommended)
USE_EMBEDDINGS = False  # distances between ref/alt final embeddings around the variant
MLM_WINDOW = (
    4  # positions each side for embedding window only (KL uses variant_span, not this).
)
WIN_LETTERS = 50  # WIN_ columns: annotation elements read this many letters each side of the variant
WIN_ELEMENTS = ("splice_donor", "splice_acceptor", "start_codon", "stop_codon")
LFC_PSEUDOCOUNT = 1.0  # LFC_ columns: log2(signal + 1) as in Borzoi's logSED; damps near-zero signal
SHARE_WINDOWS = False  # nearby variants share one reference pass (see run_inference)
COMPILE = False  # torch.compile the whole model (~1.5x after a ~1 min compile; README)
CHUNK_ROWS = 1024  # variants per checkpointed part (main writes OUTPUT_RESULTS_FILE.parts/)
SPECIES = "human"

INPUT_DATA_FILE = "clinvar_input.parquet"
OUTPUT_RESULTS_FILE = "clinvar_new_deltas.parquet"
GENOME_FILE = "hg38.fa"
METADATA_FILE = "functional_tracks_metadata_human.csv"
MODEL_NAME = "InstaDeepAI/NTv3_650M_post"  # InstaDeepAI/NTv3_100M_post


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def fetch_sequence_ucsc(chrom, start, end, genome_build="hg38"):
    """Fallback: fetch sequence from UCSC API when local genome is missing."""
    url = "https://api.genome.ucsc.edu/getData/sequence"
    params = {"genome": genome_build, "chrom": chrom, "start": start, "end": end}
    for attempt in range(3):
        try:
            resp = requests.get(url, params=params, timeout=5)
            if resp.status_code == 200:
                return resp.json()["dna"].upper()
            if resp.status_code == 429:
                time.sleep(1 + attempt)
        except Exception:
            time.sleep(0.35)
    return "N" * (end - start)


def _fetch(genome, chrom, start, end):
    """Reference letters [start, end) (short at a contig end), or None if unavailable."""
    if not genome:
        seq = fetch_sequence_ucsc(chrom, start, end)
        return None if set(seq) <= {"N"} else seq
    chrom = str(chrom)
    if chrom not in genome:
        chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    try:
        return genome[chrom][start:end].seq.upper()
    except Exception:
        return None


def _splice(seq, c, chrom, pos, ref, alt):
    """seq with `ref` at letter c replaced by `alt`; warns on a REF mismatch."""
    actual = seq[c : c + len(ref)]
    if actual.upper() != ref.upper():
        warnings.warn(f"REF mismatch at {chrom}:{pos}: expected '{ref}' but genome has '{actual}'", stacklevel=3)
    return seq[:c] + alt + seq[c + len(ref) :]


def get_genomic_sequence(genome, chrom, pos, ref, alt, context_len=1024):
    """Extract Ref/Alt sequences and return variant index in returned sequence."""
    variant_idx = pos - 1
    start = max(0, variant_idx - context_len // 2)
    ref_seq = _fetch(genome, chrom, start, variant_idx + context_len // 2)
    if not ref_seq:
        return None, None, None

    center = max(0, min(variant_idx - start, len(ref_seq) - 1))
    alt_seq = _splice(ref_seq, center, chrom, pos, ref, alt)
    target_len = min(len(ref_seq), len(alt_seq))
    center = min(center, target_len - 1)
    return ref_seq[:target_len], alt_seq[:target_len], center


KEY_MARKS = ("H3K4me3", "H3K27ac", "H3K36me3", "H3K27me3", "H3K9me3", "H3K4me1", "H3K9ac")


def get_track_indices(bigwig_names, metadata_file=None):
    """Filter BigWig tracks to a meaningful subset (<4k) via metadata."""
    if not metadata_file or not os.path.exists(metadata_file):
        print(f"⚠️  Metadata not found: {metadata_file} — using all tracks")
        return list(range(len(bigwig_names))), bigwig_names

    metadata = pd.read_csv(metadata_file)
    print(f"📂 Loaded metadata: {len(metadata)} tracks")
    metadata_lookup = metadata.drop_duplicates("file_id").set_index("file_id").to_dict("index")

    sel_idx, sel_names = [], []

    for i, tid in enumerate(bigwig_names):
        r = metadata_lookup.get(tid)
        if r is None:
            continue
        assay = str(r.get("assay", ""))
        target = str(r.get("experiment_target", ""))
        dataset = str(r.get("dataset", ""))

        keep = (
            (
                "ChIP" in assay
                and any(m in target for m in KEY_MARKS)
                and dataset in ("encode_v3", "geo")
            )
            or (("ATAC" in assay or "DNase" in assay) and dataset == "encode_v3")
            or dataset == "fantom5"
            or "CAGE" in assay
            or dataset == "gtex"
        )
        if keep:
            sel_idx.append(i)
            sel_names.append(tid)

    print(f"📊 Tracks: {len(sel_idx)}/{len(bigwig_names)} selected")
    return sel_idx, sel_names


def track_groups(track_ids, metadata_file):
    """Assay group -> positions in track_ids: each KEY_MARKS histone mark, DNase, ATAC, CAGE, RNA.

    Empty groups are left out; {} without metadata.
    """
    if not metadata_file or not os.path.exists(metadata_file):
        return {}
    meta = pd.read_csv(metadata_file).drop_duplicates("file_id").set_index("file_id")
    groups = {}
    for i, tid in enumerate(track_ids):
        if tid not in meta.index:
            continue
        assay, target = str(meta.at[tid, "assay"]), str(meta.at[tid, "experiment_target"])
        g = target if "ChIP" in assay and target in KEY_MARKS else next(
            (a for a in ("DNase", "ATAC", "CAGE", "RNA") if a in assay), None
        )
        if g:
            groups.setdefault(g, []).append(i)
    return groups


def log2_fold_change(ref, alt):
    """Per-track log2 fold change of raw BigWig signal (softplus head output, >= 0)."""
    return np.log2(alt + LFC_PSEUDOCOUNT) - np.log2(ref + LFC_PSEUDOCOUNT)


def _signed_max_abs(m):
    """Per row: the value with the largest |value| (sign kept) and its column; NaN entries skipped."""
    a = np.where(np.isnan(m), -1, np.abs(m))
    i = a.argmax(1)
    return np.where(a.max(1) >= 0, m[np.arange(len(m)), i], np.nan), i


def lfc_summary(lfc, track_ids, groups):
    """LFC_ columns from a [rows, tracks] log2_fold_change matrix (NaN tracks skipped).

    LFC_max_abs: the lfc with the largest |lfc|, signed; LFC_top_track: its track;
    LFC_mean_abs: mean |lfc|; LFC_max_abs_<group>: LFC_max_abs within each track_groups group.
    """
    top, i = _signed_max_abs(lfc)
    n = (~np.isnan(lfc)).sum(1)
    out = {
        "LFC_max_abs": top,
        "LFC_mean_abs": np.nansum(np.abs(lfc), 1) / np.where(n > 0, n, np.nan),
        "LFC_top_track": np.where(np.isnan(top), None, np.asarray(track_ids, object)[i]),
    }
    for g, idx in groups.items():
        out[f"LFC_max_abs_{g}"] = _signed_max_abs(lfc[:, idx])[0]
    return out


def clean_data(df, require_valid_pos=False):
    """Remove rows with NA ref/alt  (and pos <= 0 when require_valid_pos)."""
    n = len(df)
    bad = df["ref"].isna() | df["alt"].isna()
    if require_valid_pos:
        bad |= df["pos"] <= 0
    df_clean = df[~bad].copy().reset_index(drop=True)
    removed = n - len(df_clean)
    if removed:
        print(f"🧹 Cleaned: {n} -> {len(df_clean)} ({removed} removed)")
    return df_clean


def to_track_probabilities(track_values):
    """MAGI's track values: sigmoid of the head output (see Gotcha 4).

    BED heads output 2 classes per element ([..., 21, 2]): sigmoid of class 1.
    BigWig heads output one value per track: sigmoid of it.
    """
    # --- Binary-class output (e.g. BED: [..., num_elements, 2]) ---
    if track_values.shape[-1] == 2:
        return torch.sigmoid(track_values[..., 1])

    # --- Single-logit output (e.g. BigWig: [..., num_tracks]) ---
    return torch.sigmoid(track_values)


# ============================================================================
# MODEL LOADING
# ============================================================================
def load_model_and_tokenizer(model_name, device="cuda"):
    """Load NTv3 model + tokenizer; with COMPILE, torch.compile the whole model."""
    print(f" Loading {model_name} on {device}...")
    # TF32 matmuls on Ampere+ (fp32 accumulate, ~1e-3 relative input rounding): not byte-exact
    # with the paper's fp32 runs, but far below the benign-background noise of the deltas.
    torch.set_float32_matmul_precision("high")
    model = (
        AutoModel.from_pretrained(model_name, trust_remote_code=True).to(device).eval()
    )
    if COMPILE and str(device).startswith("cuda"):
        # One eager pass fills the lazy rotary cache; otherwise the second batch compiles again.
        with torch.inference_mode():
            model(input_ids=torch.full((1, CONTEXT_LEN), 10, device=device), species_ids=model.encode_species(["human"]).to(device))
        model = torch.compile(model, dynamic=False)  # static shapes: run_inference pads the last batch
        print("  torch.compile: the first batch compiles (~1 min)")
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    return model, tokenizer


def extract_feature_names(
    model, target_species="human", use_bed=True, use_bigwigs=True
):
    """Extract BED / BigWig feature names from model config."""
    bed_names = []
    if use_bed:
        for attr in ("bed_elements_names", "bed_tracks", "bed_track_labels"):
            if hasattr(model.config, attr):
                bed_names = getattr(model.config, attr)
                break
    bw_names = []
    if use_bigwigs and hasattr(model.config, "bigwigs_per_species"):
        bw_names = model.config.bigwigs_per_species.get(target_species, [])
    print(f"   Features: {len(bed_names)} BED, {len(bw_names)} BigWig")
    return bed_names, bw_names


# ============================================================================
# TOKENIZATION AND FORWARD PASS
# ============================================================================
def build_token_lut(tokenizer):
    """Byte -> token id, as the NTv3 tokenizer maps single letters (case-folded; anything else -> <unk>).

    ~300x faster than the tokenizer on 64 kb windows (the tokenizer took as long as the GPU forward).
    """
    vocab = tokenizer.get_vocab()
    lut = np.full(256, vocab.get("<unk>", 0), np.int64)
    for tok, tid in vocab.items():
        if len(tok) == 1:
            lut[ord(tok.upper())] = lut[ord(tok.lower())] = tid
    return lut


def tokenize(seqs, lut, length, pad_id):
    """Right-pad / truncate to `length`, as tokenizer(padding="max_length", truncation=True) does."""
    ids = np.full((len(seqs), length), pad_id, np.int64)
    for i, s in enumerate(seqs):
        b = np.frombuffer(s[:length].encode("ascii", "replace"), np.uint8)
        ids[i, : len(b)] = lut[b]
    return torch.from_numpy(ids)


def forward_trunk(model, input_ids, species_ids):
    """LM logits and the final per-letter embedding, without the track heads.

    The heads run later on the few letters we read: the stock forward runs every species'
    BigWig head over all 24k output letters (~1/3 of the forward time).
    """
    core = model.core
    heads = core.bigwig_head, core.bed_head
    core.bigwig_head = core.bed_head = None
    try:
        out = model(input_ids=input_ids, species_ids=species_ids)
    finally:
        core.bigwig_head, core.bed_head = heads
    return out.logits, out.embedding


def species_token_id(model, name):
    """NTv3 species token; raises on an unknown name (never falls back to human)."""
    return int(model.encode_species([name]).item())


def bigwig_head_for(model, species_name, n_tracks):
    """The species' BigWig head, checked by track count so a wrong head cannot pass silently."""
    if not n_tracks:
        return None
    sid = species_token_id(model, species_name) - model.config.num_species_special_tokens
    head = model.core.bigwig_head.species_heads[sid]
    n = getattr(getattr(head, "head", None), "out_features", 0)
    if n != n_tracks:
        raise ValueError(f"'{species_name}' BigWig head has {n} tracks, expected {n_tracks}")
    return head


# ============================================================================
# MLM & EMBEDDING FEATURES
# ============================================================================
def compute_mlm_features(
    ref, alt, ref_allele, alt_allele, nuc_token_map, use_kl=True, use_embeddings=False
):
    """
    Unified MLM features for SNPs and INDELs, from the letters read around the variant.

    ref / alt: per-allele reads from `run_inference`:
      logits  [span, V]  LM logits over the variant span (row 0 = variant letter)
      bases   str        the allele sequence's letters over the same span
      emb     [w, C]     final embedding over [vc - MLM_WINDOW, vc + MLM_WINDOW) (use_embeddings)
      kmer    str        5-mer centred on the variant

    Always:  LLR / MLM_Prior / MLM_Delta (NaN for non-SNPs), REF_5mer, ALT_5mer
    use_kl:  MLM_KL_mean/max = KL(alt || ref) over the span, MLM_logprob_ref/alt/delta
    use_emb: EMB_cosine_dist, EMB_l2_dist, EMB_max/mean_pos_dist
    """
    feat = {}
    ref_lp = F.log_softmax(ref["logits"].float(), dim=-1)
    alt_lp = F.log_softmax(alt["logits"].float(), dim=-1)

    # --- LLR (single-nucleotide substitutions only) ---
    is_snp = len(ref_allele) == 1 and len(alt_allele) == 1
    if is_snp and ref_allele in nuc_token_map and alt_allele in nuc_token_map:
        rp = float(ref_lp[0, nuc_token_map[ref_allele]].exp())
        ap = float(ref_lp[0, nuc_token_map[alt_allele]].exp())
        feat["LLR"] = np.log(ap / (rp + 1e-10) + 1e-10)
        feat["MLM_Prior"] = rp
        feat["MLM_Delta"] = ap - rp
    else:
        feat["LLR"] = feat["MLM_Prior"] = feat["MLM_Delta"] = np.nan

    feat["REF_5mer"] = ref["kmer"]
    feat["ALT_5mer"] = alt["kmer"]

    # --- KL divergence + log-prob over the variant span (no flanking dilution) ---
    if use_kl:
        kl = (alt_lp.exp() * (alt_lp - ref_lp)).sum(-1)
        feat["MLM_KL_mean"] = float(kl.mean())
        feat["MLM_KL_max"] = float(kl.max())
        rlp = [float(ref_lp[p, nuc_token_map[b]]) for p, b in enumerate(ref["bases"]) if b in nuc_token_map]
        alp = [float(alt_lp[p, nuc_token_map[b]]) for p, b in enumerate(alt["bases"]) if b in nuc_token_map]
        feat["MLM_logprob_ref"] = np.mean(rlp) if rlp else np.nan
        feat["MLM_logprob_alt"] = np.mean(alp) if alp else np.nan
        feat["MLM_logprob_delta"] = feat["MLM_logprob_alt"] - feat["MLM_logprob_ref"]

    # --- Embedding distances (optional) ---
    if use_embeddings:
        hr, ha = ref["emb"].float(), alt["emb"].float()
        hrm, ham = hr.mean(0), ha.mean(0)
        feat["EMB_cosine_dist"] = float(1.0 - F.cosine_similarity(hrm, ham, dim=0))
        feat["EMB_l2_dist"] = float(torch.norm(hrm - ham, p=2))
        per_pos = torch.norm(hr - ha, p=2, dim=-1)
        feat["EMB_max_pos_dist"] = float(per_pos.max())
        feat["EMB_mean_pos_dist"] = float(per_pos.mean())

    return feat


# ============================================================================
# INFERENCE PIPELINE
# ============================================================================
TOKEN_LETTERS = 128  # NTv3 pools 2**7 letters per transformer token


def _shared_windows(rows, context_len, max_span):
    """Group nearby variants under one reference window (`share_windows`).

    Only variants at the same place in a 128-letter token share: a 1-letter window shift changes
    NTv3's deltas far more than a whole-token shift (README). Rows come sorted by
    (species, chrom, pos mod 128, pos); a group spans at most `max_span` letters (half the track
    output) and its window moves from centred by whole tokens, so a lone variant gets the usual
    centred window. Yields (chrom, start, [row positions]).
    """
    def key(r):
        return r["chrom"], r.get("_sid"), (int(r["pos"]) - 1) % TOKEN_LETTERS

    i = 0
    while i < len(rows):
        first = last = int(rows[i]["pos"]) - 1
        j = i + 1
        while j < len(rows) and key(rows[j]) == key(rows[i]):
            p = int(rows[j]["pos"]) - 1
            if p + max(1, len(rows[j]["ref"])) - first > max_span:
                break
            last, j = p, j + 1
        shift = round((last - first) / 2 / TOKEN_LETTERS) * TOKEN_LETTERS
        yield rows[i]["chrom"], max(0, first + shift - context_len // 2), list(range(i, j))
        i = j


def _variant_pairs(rows, genome, context_len, share_windows, max_span):
    """(ref_seq, [(row position, alt_seq, vc)]) per reference window; vc = variant letter.

    Centred mode (default): one window per variant, as `get_genomic_sequence` builds it.
    Shared mode fills the alt after a deletion with real downstream letters (centred mode
    crops both alleles instead). Rows without a sequence are left out (features stay NaN).
    """
    if not share_windows:
        for k, row in enumerate(rows):
            r, a, vc = get_genomic_sequence(genome, row["chrom"], int(row["pos"]), row["ref"], row["alt"], context_len)
            if r:
                yield r, [(k, a, vc)]
        return
    for chrom, start, ks in _shared_windows(rows, context_len, max_span):
        g = _fetch(genome, chrom, start, start + context_len + max(len(rows[k]["ref"]) for k in ks))
        if not g:
            continue
        pairs = []
        for k in ks:
            c = int(rows[k]["pos"]) - 1 - start
            if 0 <= c < min(len(g), context_len):
                alt_seq = _splice(g, c, chrom, rows[k]["pos"], rows[k]["ref"].upper(), rows[k]["alt"].upper())
                pairs.append((k, alt_seq[:context_len], c))
        yield g[:context_len], pairs


def _win_offsets(ref, alt, allele):
    """WIN_ letters of one allele, as offsets from the variant letter.

    The WIN_LETTERS letters before the variant; the variant letters only when ref and alt have the
    same length (SNV, MNV); then WIN_LETTERS letters after the allele. Ref and alt reads line up
    entry by entry, and _win_offsets(ref, alt, ref) gives the offsets in reference letters.
    """
    m = len(allele) if len(ref) == len(alt) else 0
    j = np.arange(2 * WIN_LETTERS + m)
    return np.where(j < WIN_LETTERS + m, j - WIN_LETTERS, len(allele) + j - WIN_LETTERS - m)


def _score_rows(rows, model, ntm, lut, pad_id, cfg):
    """Forward passes and features for a list of row dicts (sorted by locus).

    Returns (track REF matrix, track delta matrix, MLM feature dicts, variant letter per row,
    LFC_ / WIN_ columns).
    """
    c, L = cfg, cfg["context_len"]
    track_len = int(L * model.config.keep_target_center_fraction)
    track_start = (L - track_len) // 2

    # One forward per unique sequence; each serves one or more (row, role) reads.
    seq_jobs = {}  # (sequence, species id) -> [(row position, role, vc)]
    for ref_seq, alts in _variant_pairs(rows, c["genome"], L, c["share_windows"], track_len // 2):
        for k, alt_seq, vc in alts:
            seq_jobs.setdefault((ref_seq, rows[k]["_sid"]), []).append((k, "ref", vc))
            seq_jobs.setdefault((alt_seq, rows[k]["_sid"]), []).append((k, "alt", vc))
    jobs = list(seq_jobs.items())

    reads = [{} for _ in rows]
    for b0 in range(0, len(jobs), c["batch_size"]):
        batch = jobs[b0 : b0 + c["batch_size"]]
        fill = batch + [batch[-1]] * (c["batch_size"] - len(batch))  # one batch shape (torch.compile)
        ids = tokenize([s for (s, _), _ in fill], lut, L, pad_id).to(c["device"])
        sp = torch.tensor([sid for (_, sid), _ in fill], device=c["device"])
        with torch.inference_mode():
            logits, emb = forward_trunk(model, ids, sp)
            track_rows, track_reads, win_rows, win_reads = [], [], [], []
            for b, ((seq, _), seq_reads) in enumerate(batch):
                for k, role, vc in seq_reads:
                    span = max(1, len(rows[k]["ref"]), len(rows[k]["alt"]))
                    read = {
                        "vc": vc,
                        "bases": seq[vc : vc + span],
                        "kmer": seq[max(0, vc - 2) : vc + 3] if len(seq) >= vc + 3 else "NNNNN",
                        "logits": logits[b, vc : vc + span].cpu(),
                    }
                    if c["use_embeddings"]:
                        read["emb"] = emb[b, max(0, vc - c["mlm_window"]) : vc + c["mlm_window"]].cpu()
                    if 0 <= vc - track_start < track_len:  # tracks: mean over the span, clipped to the output
                        track_rows.append(emb[b, vc : min(vc + span, track_start + track_len)])
                        track_reads.append(read)
                        if c["win_el"]:  # window letters inside the track output
                            p = vc + _win_offsets(rows[k]["ref"], rows[k]["alt"], rows[k][role])
                            read["win_ok"] = (p >= track_start) & (p < min(track_start + track_len, len(seq)))
                            win_rows.append(emb[b, torch.as_tensor(p[read["win_ok"]], device=emb.device)])
                            win_reads.append(read)
                    reads[k][role] = read
            heads = [h for h in (c["bed_head"], c["bw_head"]) if h is not None]
            if track_rows and heads:
                x, n = torch.cat(track_rows), [len(t) for t in track_rows]
                raw = [h(x) for h in heads]
                out = [to_track_probabilities(r) for r in raw]
                if c["bw_head"] is not None:
                    out[-1] = out[-1][:, c["bw_sel"]]
                    bw_raw = torch.stack([m.mean(0) for m in raw[-1][:, c["bw_sel"]].split(n)]).cpu()
                    for read, m in zip(track_reads, bw_raw):
                        read["bw_raw"] = m  # softplus signal, for LFC_
                means = torch.stack([m.mean(0) for m in torch.cat(out, -1).split(n)]).cpu()
                for read, m in zip(track_reads, means):
                    read["tracks"] = m
            if win_rows:  # its own head call, so the track values above stay as they were
                # standard probability (softmax over the 2 classes), not MAGI's sigmoid (AGENTS.md)
                w = c["bed_head"](torch.cat(win_rows))[:, c["win_el"]].softmax(-1)[..., 1].cpu().numpy()
                for read, part in zip(win_reads, np.split(w, np.cumsum([r["win_ok"].sum() for r in win_reads])[:-1])):
                    read["win"] = np.full((len(read["win_ok"]), w.shape[1]), np.nan)
                    read["win"][read["win_ok"]] = part

    tr_ref = np.full((len(rows), c["n_tracks"]), np.nan)
    tr_d = np.full((len(rows), c["n_tracks"]), np.nan)
    mlm, vcs = [], np.full(len(rows), np.nan)
    n_bw = len(c["bw_ids"])
    bw_ref, bw_alt = np.full((len(rows), n_bw), np.nan), np.full((len(rows), n_bw), np.nan)
    win = np.full((len(rows), len(c["win_el"]), 4), np.nan)
    for k, row in enumerate(rows):
        if "ref" not in reads[k]:
            mlm.append({})
            continue
        ref, alt = reads[k]["ref"], reads[k]["alt"]
        vcs[k] = ref["vc"]
        if c["n_tracks"] and "tracks" in ref and "tracks" in alt:
            tr_ref[k] = ref["tracks"].numpy()
            tr_d[k] = (alt["tracks"] - ref["tracks"]).numpy()
        if n_bw and "bw_raw" in ref and "bw_raw" in alt:
            bw_ref[k], bw_alt[k] = ref["bw_raw"].numpy(), alt["bw_raw"].numpy()
        if "win" in ref and "win" in alt:
            d = alt["win"] - ref["win"]  # NaN where a letter falls outside the track output
            if not np.isnan(d).all():
                g, lo, e = np.nanargmax(d, 0), np.nanargmin(d, 0), np.arange(d.shape[1])
                off = _win_offsets(row["ref"], row["alt"], row["ref"])
                gain, loss = np.maximum(d[g, e], 0), np.minimum(d[lo, e], 0)
                win[k] = np.stack([gain, loss, np.where(gain > 0, off[g], np.nan), np.where(loss < 0, off[lo], np.nan)], 1)
        mlm.append(
            compute_mlm_features(
                ref, alt, row["ref"], row["alt"], ntm, use_kl=c["use_kl"], use_embeddings=c["use_embeddings"]
            )
        )
    extra = lfc_summary(log2_fold_change(bw_ref, bw_alt), c["bw_ids"], c["lfc_groups"]) if n_bw else {}
    for e, el in enumerate(c["win_names"]):
        for s, kind in enumerate(("gain", "loss", "gain_pos", "loss_pos")):
            extra[f"WIN_{el}_{kind}"] = win[:, e, s]
    return tr_ref, tr_d, mlm, vcs, extra


def run_inference(
    df,
    model,
    tokenizer,
    genome,
    bed_names,
    bigwig_names,
    selected_bw_indices,
    device="cuda",
    *,
    extra_cols=None,
    use_kl=None,
    use_embeddings=None,
    context_len=None,
    batch_size=None,
    use_bed=None,
    use_bigwigs=None,
    mlm_window=None,
    species=None,
    share_windows=None,
    parts_dir=None,
    chunk_rows=None,
    lfc_groups=None,
):
    """
    Run inference on any variant DataFrame and compute deltas. Output rows follow df's order.

    Keyword-only params default to module-level globals.
    extra_cols: additional df columns to preserve in output.
    species: NTv3 species token string (e.g. 'canis_lupus_familiaris'); bigwig_names must belong
             to it. Defaults to 'human'. 'auto' reads df['species'] (no BigWig tracks then).
    share_windows: nearby variants share one reference pass (see _shared_windows and README:
             features move slightly, so never mix the two modes in one cohort). Adds `window_offset`.
    parts_dir: write each chunk of `chunk_rows` variants there and skip chunks already written,
             so a preempted run resumes; a run with other variants or settings refuses them.
    lfc_groups: track_groups() of the selected tracks, for the LFC_max_abs_<group> columns.
    """
    use_bed = USE_BED if use_bed is None else use_bed
    use_bigwigs = USE_BIGWIGS if use_bigwigs is None else use_bigwigs
    chunk_rows = CHUNK_ROWS if chunk_rows is None else chunk_rows
    extra_cols = extra_cols or []
    auto_species = species == "auto" and "species" in df.columns
    run_species = "human" if species in (None, "auto") else species
    bw_head = bigwig_head_for(model, run_species, len(bigwig_names) if use_bigwigs else 0)
    if auto_species and bw_head is not None:
        raise ValueError("species='auto' has no BigWig tracks (each species has its own head); pass use_bigwigs=False")
    bed_head = model.core.bed_head if use_bed and bed_names else None
    names = [f"BED_{nm}" for nm in (bed_names if bed_head is not None else [])]
    bw_ids = [bigwig_names[gi] for gi in (selected_bw_indices if bw_head is not None else [])]
    names += [f"BW_{t}" for t in bw_ids]
    win_names = [e for e in WIN_ELEMENTS if e in (bed_names if bed_head is not None else [])]
    cfg = dict(
        genome=genome,
        device=device,
        context_len=CONTEXT_LEN if context_len is None else context_len,
        batch_size=BATCH_SIZE if batch_size is None else batch_size,
        mlm_window=MLM_WINDOW if mlm_window is None else mlm_window,
        use_kl=USE_KL_DIVERGENCE if use_kl is None else use_kl,
        use_embeddings=USE_EMBEDDINGS if use_embeddings is None else use_embeddings,
        share_windows=SHARE_WINDOWS if share_windows is None else share_windows,
        bed_head=bed_head,
        bw_head=bw_head,
        bw_sel=torch.as_tensor(selected_bw_indices, dtype=torch.long, device=device),
        n_tracks=len(names),
        bw_ids=bw_ids,
        lfc_groups=lfc_groups or {},
        win_names=win_names,
        win_el=[bed_names.index(e) for e in win_names],
    )
    context_len, share_windows = cfg["context_len"], cfg["share_windows"]
    lut = build_token_lut(tokenizer)
    pad_id = tokenizer.pad_token_id if tokenizer.pad_token_id is not None else 1
    ntm = {n: int(lut[ord(n)]) for n in "ACGTN"}
    print(
        f"🧠 {len(df)} variants  (KL={cfg['use_kl']} EMB={cfg['use_embeddings']} ctx={context_len // 1024}k"
        f" shared windows={share_windows})"
    )

    # Sorted by locus (and by place in a token when windows are shared): neighbours share
    # windows, identical reference windows run once.
    df = df.reset_index(drop=True)
    keys = (["species"] if auto_species else []) + ["chrom"] + (["_phase"] if share_windows else []) + ["pos"]
    sdf = df.assign(_phase=(df["pos"] - 1) % TOKEN_LETTERS).sort_values(
        keys, kind="stable", key=lambda s: s.astype(str) if s.dtype == object else s
    ).drop(columns="_phase")
    sid = species_token_id(model, run_species)
    sdf = sdf.assign(
        _row=sdf.index,
        _sid=sdf["species"].astype(str).map(lambda s: species_token_id(model, s)) if auto_species else sid,
        ref=sdf["ref"].astype(str),
        alt=sdf["alt"].astype(str),
    )
    if parts_dir:
        variants = pd.util.hash_pandas_object(df[["chrom", "pos", "ref", "alt"]].astype(str), index=False)
        manifest = dict(
            variants=hashlib.sha1(variants.to_numpy().tobytes()).hexdigest(),
            columns=names + sorted(cfg["lfc_groups"]) + cfg["win_names"],
            model=str(getattr(model, "name_or_path", "")),
            compiled=hasattr(model, "_orig_mod"),
            genome=str(getattr(genome, "filename", genome)),
            species=str(species),
            chunk_rows=chunk_rows,
            **{k: v for k, v in cfg.items() if isinstance(v, (bool, int, str)) and k != "device"},
        )
        os.makedirs(parts_dir, exist_ok=True)
        path = os.path.join(parts_dir, "manifest.json")
        if os.path.exists(path):
            with open(path) as f:
                if json.load(f) != manifest:
                    raise ValueError(f"{parts_dir} holds parts of a different run; delete it to start over")
        else:
            with open(path, "w") as f:
                json.dump(manifest, f)

    out, failed = [], 0
    for ci, c0 in enumerate(tqdm(range(0, len(sdf), chunk_rows), desc="chunks")):
        part = os.path.join(parts_dir, f"part-{ci:05d}.parquet") if parts_dir else None
        if part and os.path.exists(part):
            out.append(pd.read_parquet(part))
            continue
        chunk = sdf.iloc[c0 : c0 + chunk_rows]
        tr_ref, tr_d, mlm, vcs, extra = _score_rows(chunk.to_dict("records"), model, ntm, lut, pad_id, cfg)
        failed += int(np.isnan(vcs).sum())
        cols = {c: chunk[c].to_numpy() for c in ("chrom", "pos", "ref", "alt", "label") if c in df.columns}
        for c in extra_cols:
            cols[c] = chunk[c].to_numpy() if c in df.columns else np.nan
        cols["indel_size"] = (chunk["alt"].str.len() - chunk["ref"].str.len()).to_numpy()
        if share_windows:
            cols["window_offset"] = vcs - context_len // 2
        for j, nm in enumerate(names):
            cols[f"REF_{nm}"] = tr_ref[:, j]
            cols[f"D_{nm}"] = tr_d[:, j]
        res = pd.concat([pd.DataFrame(cols), pd.DataFrame(mlm, index=range(len(chunk))), pd.DataFrame(extra)], axis=1)
        res["_row"] = chunk["_row"].to_numpy()
        if part:
            res.to_parquet(part + ".tmp", index=False)
            os.replace(part + ".tmp", part)  # a preempted write never leaves a truncated part
        out.append(res)
    if failed:
        print(f"⚠️  {failed} variants: no sequence (or outside a shared window); their features are NaN")
    res = pd.concat(out, ignore_index=True)
    del out
    return res.sort_values("_row").drop(columns="_row").reset_index(drop=True)


# ============================================================================
# VALIDATION
# ============================================================================
def validate_output(df):
    """Print summary of inference results (works for both SNPs and INDELs)."""
    print("\n" + "=" * 80)
    print("🔍 OUTPUT VALIDATION")
    print("=" * 80)

    # MLM
    print("\n  MLM Features:")
    for col in (
        "LLR",
        "MLM_Prior",
        "MLM_Delta",
        "MLM_KL_mean",
        "MLM_KL_max",
        "MLM_logprob_ref",
        "MLM_logprob_alt",
        "MLM_logprob_delta",
    ):
        if col not in df.columns:
            continue
        v = df[col].dropna()
        tag = f"{len(v)}/{len(df)}" if len(v) < len(df) else "all"
        if len(v):
            print(
                f"   ✅ {col:22} | {tag} | [{v.min():.4f}, {v.max():.4f}] mean={v.mean():.4f}"
            )
        else:
            print(f"   ⚠️  {col:22} | all NaN")

    # Embeddings
    ecols = [c for c in df.columns if c.startswith("EMB_")]
    if ecols:
        print("\n  Embedding Features:")
        for c in ecols:
            v = df[c].dropna()
            if len(v):
                print(f"   ✅ {c:22} | [{v.min():.4f}, {v.max():.4f}]")
            else:
                print(f"   ⚠️  {c:22} | all NaN")

    # BED / BigWig
    for pfx, label in [("REF_BED_", "BED"), ("REF_BW_", "BigWig")]:
        rcols = [c for c in df.columns if c.startswith(pfx)]
        dcols = [c for c in df.columns if c.startswith(pfx.replace("REF_", "D_"))]
        if rcols:
            vals = df[rcols].values.flatten()
            prob = " ✅ probs" if vals.min() >= 0 and vals.max() <= 1 else ""
            print(
                f"\n   {label}: {len(rcols)} REF + {len(dcols)} Delta"
                f"  range=[{vals.min():.4f},{vals.max():.4f}]{prob}"
            )

    # Indel size
    if "indel_size" in df.columns and (df["indel_size"] != 0).any():
        s = df["indel_size"]
        print(
            f"\n   Indel sizes: [{s.min()},{s.max()}]  del={(s < 0).sum()} ins={(s > 0).sum()} sub={(s == 0).sum()}"
        )

    # Metadata
    for c in ("variant_type", "gene", "variant_id", "quality_stars"):
        if c in df.columns:
            print(f"   {c}: {df[c].notna().sum()}/{len(df)}")

    if "label" in df.columns:
        print(f"\n   Labels: {dict(df['label'].value_counts())}")
    print("=" * 80)


# ============================================================================
# MAIN
# ============================================================================
def main():
    """Variant inference pipeline (default entry point; settings are the module constants)."""
    print("=" * 80 + "\n NTv3 Variant Inference\n" + "=" * 80)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        print(f"  GPU: {torch.cuda.get_device_name(0)}")

    print(f"  DEBUG={DEBUG_MODE}  KL={USE_KL_DIVERGENCE}  EMB={USE_EMBEDDINGS}  shared windows={SHARE_WINDOWS}")
    print(f"  Input: {INPUT_DATA_FILE}   Output: {OUTPUT_RESULTS_FILE}")

    if not os.path.exists(INPUT_DATA_FILE):
        print(f"❌ Not found: {INPUT_DATA_FILE}")
        return

    df = pd.read_parquet(INPUT_DATA_FILE)
    labels = dict(df["label"].value_counts()) if "label" in df.columns else "no label column"
    print(f"  {len(df)} variants | {labels}")

    base = {"chrom", "pos", "ref", "alt", "label"}
    extra = [c for c in df.columns if c not in base]
    if extra:
        print(f"  Extra cols: {extra}")

    df = clean_data(df, require_valid_pos=(df["pos"] <= 0).any())
    if DEBUG_MODE:
        df = df.head(33)
        print(f"🐞 DEBUG: {len(df)} variants")

    genome = Fasta(GENOME_FILE) if os.path.exists(GENOME_FILE) else None
    if genome is None and SPECIES != "human":  # the fallback fetches hg38 from UCSC
        print(f"❌ Not found: {GENOME_FILE} (needed for species {SPECIES})")
        return
    model, tokenizer = load_model_and_tokenizer(MODEL_NAME, device)
    bed_names, bw_names = extract_feature_names(
        model, target_species=SPECIES, use_bed=USE_BED, use_bigwigs=USE_BIGWIGS and SPECIES != "auto"
    )

    if USE_BIGWIGS and bw_names:
        sel_idx, sel_names = get_track_indices(bw_names, METADATA_FILE)
    else:
        sel_idx, sel_names = list(range(len(bw_names))), bw_names
    groups = track_groups(sel_names, METADATA_FILE)

    print(f"  {len(sel_names)} BW selected, {len(bed_names)} BED")

    parts_dir = OUTPUT_RESULTS_FILE + ".parts"  # resume point if the job is preempted
    rdf = run_inference(
        df,
        model,
        tokenizer,
        genome,
        bed_names,
        bw_names,
        sel_idx,
        device,
        extra_cols=extra,
        species=SPECIES,
        parts_dir=parts_dir,
        lfc_groups=groups,
    )
    validate_output(rdf)

    print(f"\n💾 Saving {rdf.shape} -> {OUTPUT_RESULTS_FILE}")
    rdf.to_parquet(OUTPUT_RESULTS_FILE, index=False)
    shutil.rmtree(parts_dir)
    print(f"✅ {len(rdf)} variants, {rdf.shape[1]} columns")
    return rdf


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="NTv3 variant inference. Defaults are the module constants.")
    p.add_argument("--input", default=INPUT_DATA_FILE)
    p.add_argument("--output", default=OUTPUT_RESULTS_FILE)
    p.add_argument("--model", default=MODEL_NAME)
    p.add_argument("--genome", default=GENOME_FILE)
    p.add_argument("--species", default=SPECIES, help="NTv3 species token (needs --genome unless human)")
    p.add_argument("--context-len", type=int, default=CONTEXT_LEN)
    p.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    p.add_argument("--embeddings", action=argparse.BooleanOptionalAction, default=USE_EMBEDDINGS)
    p.add_argument("--bigwigs", action=argparse.BooleanOptionalAction, default=USE_BIGWIGS)
    p.add_argument("--share-windows", action=argparse.BooleanOptionalAction, default=SHARE_WINDOWS)
    p.add_argument("--compile", action=argparse.BooleanOptionalAction, default=COMPILE)
    p.add_argument("--debug", action="store_true", default=DEBUG_MODE)
    a = p.parse_args()
    if a.species != "human" and a.genome == GENOME_FILE:
        p.error(f"--species {a.species} needs --genome (the default {GENOME_FILE} is human)")
    INPUT_DATA_FILE, OUTPUT_RESULTS_FILE, MODEL_NAME, GENOME_FILE = a.input, a.output, a.model, a.genome
    SPECIES, CONTEXT_LEN, BATCH_SIZE, USE_EMBEDDINGS = a.species, a.context_len, a.batch_size, a.embeddings
    USE_BIGWIGS, SHARE_WINDOWS, COMPILE, DEBUG_MODE = a.bigwigs, a.share_windows, a.compile, a.debug
    main()
