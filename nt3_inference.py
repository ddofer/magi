#!/usr/bin/env python3
"""
NTv3 Genomic Variant Inference Pipeline (shared module)

Purpose
  Predict functional track changes for genomic variants using NTv3 and compare
  ref vs alt signals. Supports SNPs and INDELs using a shared pipeline.

Inputs
  - Parquet with required columns: chrom, pos (1-based), ref, alt, label
  - Optional extra columns are preserved in the output
  - Reference genome FASTA (hg38.fa) or UCSC fallback
  - Track metadata CSV for BigWig filtering (functional_tracks_metadata_human.csv)

Outputs
  - Parquet with ref/alt deltas for BED and BigWig tracks
  - MLM features:
      LLR, MLM_Prior, MLM_Delta (SNPs only)
      MLM_KL_mean, MLM_KL_max (KL(alt || ref) across window)
      MLM_logprob_ref/alt/delta (sequence log-prob around variant)
  - Embedding features (optional): EMB_* distances from hidden states
  - Metadata: indel_size plus any extra input columns

Method notes
  - Ref/alt sequences are centered at pos and cropped to equal length
  - All features (BED, BigWig, MLM/KL) use variant-span extraction:
    single position for SNPs, mean over max(len(ref), len(alt)) for INDELs
  - Track heads output at ~0.375x input resolution (deconv); center = variant pos
  - KL features capture surprise/novelty in prediction distributions

Quick start
  - SNPs:   python nt3_inference.py
  - INDELs: python nt3_inference_indel.py

Examples
  - Override config in Python:
      import nt3_inference as base
      base.DEBUG_MODE = True
      base.OUTPUT_RESULTS_FILE = "tmp.parquet"
      base.main()
  - INDEL run with embeddings:
      import nt3_inference_indel  # sets base.USE_EMBEDDINGS = True
      import nt3_inference as base
      base.DEBUG_MODE = True
      base.main()

Config
  - CONTEXT_LEN, BATCH_SIZE, USE_BED, USE_BIGWIGS
  - USE_KL_DIVERGENCE, USE_EMBEDDINGS, MLM_WINDOW

Dependencies
  torch, transformers, pandas, numpy, pyfaidx, tqdm, requests

Gotchas
  1) Species ID: must use model.encode_species(["human"])
  2) BED tracks can be missing: use getattr(output, "bed_tracks_logits", None)
  3) Use AutoModel (not AutoModelForMaskedLM) for NTv3
  4) BigWig/BED outputs are converted to probabilities via sigmoid
     BED tracks have shape (B, L', num_elements, 2): extract positive-class
     logit [..., 1] then apply sigmoid.  BigWig tracks have shape
     (B, L', num_tracks): apply sigmoid directly.  Never use softmax —
     all track heads are multilabel / independent.

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

import os
import time
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
USE_EMBEDDINGS = False  # hidden-state distances  (needs output_hidden_states)
MLM_WINDOW = (
    4  # positions each side for embedding window only (KL uses variant_span, not this).
)

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


def get_genomic_sequence(genome, chrom, pos, ref, alt, context_len=1024):
    """Extract Ref/Alt sequences and return variant index in returned sequence."""
    variant_idx = pos - 1
    half = context_len // 2
    start = max(0, variant_idx - half)
    end = variant_idx + half

    if genome:
        try:
            if chrom not in genome:
                chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            ref_seq = genome[chrom][start:end].seq.upper()
        except Exception:
            return None, None, None
    else:
        ref_seq = fetch_sequence_ucsc(chrom, start, end)

    if not ref_seq:
        fallback = "N" * context_len
        return fallback, fallback, context_len // 2

    center = variant_idx - start
    center = max(0, min(center, len(ref_seq) - 1))
    # Validate REF allele against genome
    ref_end = min(center + len(ref), len(ref_seq))
    actual_ref = ref_seq[center:ref_end]
    if actual_ref.upper() != ref.upper():
        import warnings

        warnings.warn(
            f"REF mismatch at {chrom}:{pos}: expected '{ref}' but genome has '{actual_ref}'",
            stacklevel=2,
        )
    alt_seq = ref_seq[:center] + alt + ref_seq[ref_end:]

    target_len = min(len(ref_seq), len(alt_seq))
    center = min(center, target_len - 1)
    return ref_seq[:target_len], alt_seq[:target_len], center


def get_track_indices(bigwig_names, metadata_file=None):
    """Filter BigWig tracks to a meaningful subset (<4k) via metadata."""
    if not metadata_file or not os.path.exists(metadata_file):
        print(f"⚠️  Metadata not found: {metadata_file} — using all tracks")
        return list(range(len(bigwig_names))), bigwig_names

    metadata = pd.read_csv(metadata_file)
    print(f"📂 Loaded metadata: {len(metadata)} tracks")

    key_marks = {
        "H3K4me3",
        "H3K27ac",
        "H3K36me3",
        "H3K27me3",
        "H3K9me3",
        "H3K4me1",
        "H3K9ac",
    }
    sel_idx, sel_names = [], []

    for i, tid in enumerate(bigwig_names):
        rows = metadata[metadata["file_id"] == tid]
        if rows.empty:
            continue
        r = rows.iloc[0]
        assay = str(r.get("assay", ""))
        target = str(r.get("experiment_target", ""))
        dataset = str(r.get("dataset", ""))

        keep = (
            (
                "ChIP" in assay
                and any(m in target for m in key_marks)
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


def build_nuc_token_map(tokenizer):
    """Pre-compute nucleotide -> token ID mapping."""
    return {
        nuc: tokenizer(nuc, add_special_tokens=False)["input_ids"][0]
        for nuc in "ACGTN"
        if tokenizer(nuc, add_special_tokens=False)["input_ids"]
    }


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
    """Convert NTv3 track logits to probabilities.

    NTv3 post-trained track heads output logits; convert with sigmoid.
    Never use softmax (tracks are multilabel/independent).

    Rules:
    - single logit per track, e.g. BigWig: apply sigmoid directly.
    - If last dim is 2: binary classification per element
      (confirmed: BED shape is (B, L', 21, 2); index 1 = positive class).
      Extract the positive-class logit ([..., 1]) then apply sigmoid.

    """
    if track_values is None:
        return None

    # --- Binary-class output (e.g. BED: [..., num_elements, 2]) ---
    if track_values.shape[-1] == 2:
        return torch.sigmoid(track_values[..., 1])

    # --- Single-logit output (e.g. BigWig: [..., num_tracks]) ---
    return torch.sigmoid(track_values)


# ============================================================================
# MODEL LOADING
# ============================================================================
def load_model_and_tokenizer(model_name, device="cuda"):
    """Load NTv3 model + tokenizer."""
    print(f" Loading {model_name} on {device}...")
    model = (
        AutoModel.from_pretrained(model_name, trust_remote_code=True).to(device).eval()
    )
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
# MLM & EMBEDDING FEATURE EXTRACTION
# ============================================================================
def compute_mlm_features(
    out_ref,
    out_alt,
    ref_seq,
    alt_seq,
    idx,
    variant_center,
    ref_allele,
    alt_allele,
    nuc_token_map,
    use_kl=True,
    use_embeddings=False,
    window=50,
):
    """
    Unified MLM features for SNPs and INDELs.

    Always:  LLR / MLM_Prior / MLM_Delta (NaN for non-SNPs), REF_5mer, ALT_5mer
    use_kl:  MLM_KL_mean/max, MLM_logprob_ref/alt/delta
           KL is computed as KL(alt || ref) across the window.
    use_emb: EMB_cosine_dist, EMB_l2_dist, EMB_max/mean_pos_dist
           Requires output_hidden_states=True (or hidden_states present).
    """
    feat = {}
    ref_logits = out_ref.logits[idx]
    alt_logits = out_alt.logits[idx]
    seq_len = min(ref_logits.shape[0], alt_logits.shape[0])

    # KL/logprob window: cover exactly the variant span (1 pos for SNPs, full INDEL for INDELs)
    variant_span = max(1, len(ref_allele), len(alt_allele))
    variant_center = int(max(0, min(variant_center, seq_len - 1)))
    kl_ws = variant_center
    kl_we = min(seq_len, variant_center + variant_span)
    # Embedding window: symmetric flanking context (unchanged)
    ws = max(0, variant_center - window)
    we = min(seq_len, variant_center + window)

    # --- LLR (single-nucleotide substitutions only) ---
    is_snp = len(ref_allele) == 1 and len(alt_allele) == 1
    if is_snp and ref_allele in nuc_token_map and alt_allele in nuc_token_map:
        probs = F.softmax(ref_logits[variant_center], dim=-1)
        rp = float(probs[nuc_token_map[ref_allele]].cpu())
        ap = float(probs[nuc_token_map[alt_allele]].cpu())
        feat["LLR"] = np.log(ap / (rp + 1e-10) + 1e-10)
        feat["MLM_Prior"] = rp
        feat["MLM_Delta"] = ap - rp
    else:
        feat["LLR"] = feat["MLM_Prior"] = feat["MLM_Delta"] = np.nan

    # --- Context k-mers ---
    for pf, seq in [("REF", ref_seq), ("ALT", alt_seq)]:
        feat[f"{pf}_5mer"] = (
            seq[max(0, variant_center - 2) : variant_center + 3]
            if len(seq) >= variant_center + 3
            else "NNNNN"
        )

    # --- KL divergence + log-prob (optional) ---
    if use_kl:
        if kl_we <= kl_ws:
            feat["MLM_KL_mean"] = np.nan
            feat["MLM_KL_max"] = np.nan
            feat["MLM_logprob_ref"] = np.nan
            feat["MLM_logprob_alt"] = np.nan
            feat["MLM_logprob_delta"] = np.nan
            return feat

        # Compute KL over the variant span (no flanking dilution)
        rp_w = F.softmax(ref_logits[kl_ws:kl_we], dim=-1)
        ap_w = F.softmax(alt_logits[kl_ws:kl_we], dim=-1)
        kl = F.kl_div(rp_w.log(), ap_w, reduction="none", log_target=False).sum(-1)
        feat["MLM_KL_mean"] = float(kl.mean().cpu())
        feat["MLM_KL_max"] = float(kl.max().cpu())

        # Reuse softmax tensors for log-prob (avoid redundant softmax calls)
        rlp = [
            float(torch.log(rp_w[p - kl_ws, nuc_token_map[ref_seq[p]]] + 1e-10).cpu())
            for p in range(kl_ws, kl_we)
            if p < len(ref_seq) and ref_seq[p] in nuc_token_map
        ]
        alp = [
            float(torch.log(ap_w[p - kl_ws, nuc_token_map[alt_seq[p]]] + 1e-10).cpu())
            for p in range(kl_ws, kl_we)
            if p < len(alt_seq) and alt_seq[p] in nuc_token_map
        ]
        feat["MLM_logprob_ref"] = np.mean(rlp) if rlp else np.nan
        feat["MLM_logprob_alt"] = np.mean(alp) if alp else np.nan
        feat["MLM_logprob_delta"] = feat["MLM_logprob_alt"] - feat["MLM_logprob_ref"]

    # --- Embedding distances (optional) ---
    if use_embeddings:
        hr = getattr(out_ref, "last_hidden_state", None)
        ha = getattr(out_alt, "last_hidden_state", None)
        if hr is not None and ha is not None:
            hr, ha = hr[idx, ws:we, :], ha[idx, ws:we, :]
            hrm, ham = hr.mean(0), ha.mean(0)
            feat["EMB_cosine_dist"] = float(
                1.0 - F.cosine_similarity(hrm.unsqueeze(0), ham.unsqueeze(0)).cpu()
            )
            feat["EMB_l2_dist"] = float(torch.norm(hrm - ham, p=2).cpu())
            per_pos = torch.norm(hr - ha, p=2, dim=-1)
            feat["EMB_max_pos_dist"] = float(per_pos.max().cpu())
            feat["EMB_mean_pos_dist"] = float(per_pos.mean().cpu())
        else:
            for k in (
                "EMB_cosine_dist",
                "EMB_l2_dist",
                "EMB_max_pos_dist",
                "EMB_mean_pos_dist",
            ):
                feat[k] = np.nan

    return feat


# ============================================================================
# INFERENCE PIPELINE
# ============================================================================
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
):
    """
    Run inference on any variant DataFrame and compute deltas.

    Keyword-only params default to module-level globals.
    extra_cols: additional df columns to preserve in output.
    species: NTv3 species token string (e.g. 'canis_lupus_familiaris').
             Defaults to 'human'. Pass 'auto' to read from df['species'] column
             (requires a 'species' column with consistent NTv3 token names).
    """
    # defaults from module globals
    if use_kl is None:
        use_kl = USE_KL_DIVERGENCE
    if use_embeddings is None:
        use_embeddings = USE_EMBEDDINGS
    if context_len is None:
        context_len = CONTEXT_LEN
    if batch_size is None:
        batch_size = BATCH_SIZE
    if use_bed is None:
        use_bed = USE_BED
    if use_bigwigs is None:
        use_bigwigs = USE_BIGWIGS
    if mlm_window is None:
        mlm_window = MLM_WINDOW
    extra_cols = extra_cols or []

    # Species ID (default: human for backward-compatibility)
    _species_name = species if species and species != "auto" else "human"
    try:
        _species_id_default = model.encode_species([_species_name]).item()
    except (AttributeError, TypeError):
        _species_id_default = model.config.species_to_token_id.get(_species_name, 27)
    _auto_species = species == "auto" and "species" in df.columns

    ntm = build_nuc_token_map(tokenizer)
    need_hs = use_embeddings
    results = []
    print(
        f"🧠 {len(df)} variants  (KL={use_kl} EMB={use_embeddings} ctx={context_len // 1024}k)"
    )

    for i in tqdm(range(0, len(df), batch_size)):
        batch = df.iloc[i : i + batch_size]
        rseqs, aseqs, vcenters = [], [], []
        for _, row in batch.iterrows():
            r, a, vc = get_genomic_sequence(
                genome, row["chrom"], row["pos"], row["ref"], row["alt"], context_len
            )
            rseqs.append(r if r else "N" * context_len)
            aseqs.append(a if a else "N" * context_len)
            vcenters.append(vc if vc is not None else context_len // 2)

        tok_kw = dict(
            return_tensors="pt",
            padding="max_length",
            max_length=context_len,
            truncation=True,
            add_special_tokens=False,
            pad_to_multiple_of=128,
        )
        inp_r = tokenizer(rseqs, **tok_kw).to(device)
        inp_a = tokenizer(aseqs, **tok_kw).to(device)

        with torch.no_grad():
            if _auto_species:
                batch_species_ids = []
                for _, _row in batch.iterrows():
                    _sname = str(_row.get("species", _species_name))
                    try:
                        _sid = model.encode_species([_sname]).item()
                    except Exception:
                        _sid = model.config.species_to_token_id.get(
                            _sname, _species_id_default
                        )
                    batch_species_ids.append(_sid)
                sp = torch.tensor(batch_species_ids, device=device)
            else:
                sp = torch.tensor([_species_id_default] * len(rseqs), device=device)
            if need_hs:
                try:
                    out_r = model(**inp_r, species_ids=sp, output_hidden_states=True)
                    out_a = model(**inp_a, species_ids=sp, output_hidden_states=True)
                    for o in (out_r, out_a):
                        if (
                            getattr(o, "last_hidden_state", None) is None
                            and hasattr(o, "hidden_states")
                            and o.hidden_states
                        ):
                            o.last_hidden_state = o.hidden_states[-1]
                except TypeError:
                    out_r = model(**inp_r, species_ids=sp)
                    out_a = model(**inp_a, species_ids=sp)
            else:
                out_r = model(**inp_r, species_ids=sp)
                out_a = model(**inp_a, species_ids=sp)

        for idx, (_, row) in enumerate(batch.iterrows()):
            res = {c: row[c] for c in ("chrom", "pos", "ref", "alt", "label")}
            for c in extra_cols:
                res[c] = row.get(c, np.nan)
            ref_allele = str(row["ref"])
            alt_allele = str(row["alt"])
            res["indel_size"] = len(alt_allele) - len(ref_allele)
            # Variant span in bases (1 for SNPs, max allele len for INDELs)
            variant_span = max(1, max(len(ref_allele), len(alt_allele)))
            vcenter = int(vcenters[idx])
            in_len = int(inp_r["input_ids"].shape[1])

            # === BED ===
            bed_r = getattr(out_r, "bed_tracks_logits", None)
            bed_a = getattr(out_a, "bed_tracks_logits", None)
            if use_bed and bed_r is not None and bed_a is not None:
                bed_r_probs = to_track_probabilities(bed_r[idx])
                bed_a_probs = to_track_probabilities(bed_a[idx])
                track_len = int(bed_r_probs.shape[0])
                track_start = max(0, (in_len - track_len) // 2)
                bed_pos = vcenter - track_start
                if 0 <= bed_pos < track_len:
                    be = min(bed_pos + variant_span, track_len)
                    br = bed_r_probs[bed_pos:be].mean(0).cpu().numpy()
                    ba = bed_a_probs[bed_pos:be].mean(0).cpu().numpy()
                    for j, nm in enumerate(bed_names):
                        res[f"REF_BED_{nm}"] = float(br[j])
                        res[f"D_BED_{nm}"] = float(ba[j] - br[j])
                else:
                    for nm in bed_names:
                        res[f"REF_BED_{nm}"] = np.nan
                        res[f"D_BED_{nm}"] = np.nan

            # === BigWig ===
            bw_r_logits = getattr(out_r, "bigwig_tracks_logits", None)
            bw_a_logits = getattr(out_a, "bigwig_tracks_logits", None)
            if use_bigwigs and bw_r_logits is not None and bw_a_logits is not None:
                bw_r_probs = to_track_probabilities(bw_r_logits[idx])
                bw_a_probs = to_track_probabilities(bw_a_logits[idx])
                track_len = int(bw_r_probs.shape[0])
                track_start = max(0, (in_len - track_len) // 2)
                bw_pos = vcenter - track_start
                if 0 <= bw_pos < track_len:
                    bwe = min(bw_pos + variant_span, track_len)
                    bwr = bw_r_probs[bw_pos:bwe].mean(0).cpu().numpy()
                    bwa = bw_a_probs[bw_pos:bwe].mean(0).cpu().numpy()
                    for gi in selected_bw_indices:
                        res[f"REF_BW_{bigwig_names[gi]}"] = float(bwr[gi])
                        res[f"D_BW_{bigwig_names[gi]}"] = float(bwa[gi] - bwr[gi])
                else:
                    for gi in selected_bw_indices:
                        res[f"REF_BW_{bigwig_names[gi]}"] = np.nan
                        res[f"D_BW_{bigwig_names[gi]}"] = np.nan

            # === MLM ===
            res.update(
                compute_mlm_features(
                    out_r,
                    out_a,
                    rseqs[idx],
                    aseqs[idx],
                    idx,
                    vcenter,
                    ref_allele,
                    alt_allele,
                    ntm,
                    use_kl=use_kl,
                    use_embeddings=use_embeddings,
                    window=mlm_window,
                )
            )
            results.append(res)

    return pd.DataFrame(results)


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

    print(f"\n   Labels: {dict(df['label'].value_counts())}")
    print("=" * 80)


# ============================================================================
# MAIN
# ============================================================================
def main():
    """SNP inference pipeline (default entry point)."""
    print("=" * 80 + "\n NTv3 Variant Inference\n" + "=" * 80)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        print(f"  GPU: {torch.cuda.get_device_name(0)}")

    print(f"  DEBUG={DEBUG_MODE}  KL={USE_KL_DIVERGENCE}  EMB={USE_EMBEDDINGS}")
    print(f"  Input: {INPUT_DATA_FILE}   Output: {OUTPUT_RESULTS_FILE}")

    if not os.path.exists(INPUT_DATA_FILE):
        print(f"❌ Not found: {INPUT_DATA_FILE}")
        return

    df = pd.read_parquet(INPUT_DATA_FILE)
    print(f"  {len(df)} variants | {dict(df['label'].value_counts())}")

    base = {"chrom", "pos", "ref", "alt", "label"}
    extra = [c for c in df.columns if c not in base]
    if extra:
        print(f"  Extra cols: {extra}")

    df = clean_data(df, require_valid_pos=(df["pos"] <= 0).any())
    if DEBUG_MODE:
        df = df.head(33)
        print(f"🐞 DEBUG: {len(df)} variants")

    genome = Fasta(GENOME_FILE) if os.path.exists(GENOME_FILE) else None
    model, tokenizer = load_model_and_tokenizer(MODEL_NAME, device)
    bed_names, bw_names = extract_feature_names(
        model, use_bed=USE_BED, use_bigwigs=USE_BIGWIGS
    )

    if USE_BIGWIGS and bw_names:
        sel_idx, sel_names = get_track_indices(bw_names, METADATA_FILE)
    else:
        sel_idx, sel_names = list(range(len(bw_names))), bw_names

    print(f"  {len(sel_names)} BW selected, {len(bed_names)} BED")

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
    )
    validate_output(rdf)

    print(f"\n💾 Saving {rdf.shape} -> {OUTPUT_RESULTS_FILE}")
    rdf.to_parquet(OUTPUT_RESULTS_FILE, index=False)
    print(f"✅ {len(rdf)} variants, {rdf.shape[1]} columns")
    return rdf


if __name__ == "__main__":
    main()
