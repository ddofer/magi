#!/usr/bin/env python3
"""
NTv3 inference module for the MAGI Gradio app.

Adapted from the top-level inference.py pipeline for web deployment.

Key features:
- Uses local hg38.fa via pyfaidx when available
- Loads the configured NTv3 model once and caches it across requests
- Returns BED/BigWig deltas together with sequence-model metrics

Usage:
    from inference import predict_variants
    results_df = predict_variants(variants_df)
"""

import os
import time
import warnings
from pathlib import Path
from typing import Optional, Tuple, List, Any, Dict, cast
from functools import lru_cache

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import requests
from pyfaidx import Fasta
from transformers import AutoModel, AutoTokenizer

# ============================================================================
# CONFIGURATION
# ============================================================================
MODEL_NAME = "InstaDeepAI/NTv3_650M_post"
CONTEXT_LEN = 32 * 1024  # 32 kb sequence window
USE_BED = True
USE_BIGWIGS = True
USE_KL_DIVERGENCE = True
USE_EMBEDDINGS = True  # Useful for indels
MLM_WINDOW = 3  # ±3 positions for embedding window

# Paths (relative to this file)
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
GENOME_FILE = DATA_DIR / "hg38.fa"
GENOME_GZ_FILE = DATA_DIR / "hg38.fa.gz"
METADATA_FILE = DATA_DIR / "functional_tracks_metadata_human.csv"
UCSC_SEQUENCE_URL = "https://api.genome.ucsc.edu/getData/sequence"
ENSEMBL_SEQUENCE_URL = "https://rest.ensembl.org/sequence/region"
ENSEMBL_PLANTS_SEQUENCE_URL = "https://rest.plants.ensembl.org/sequence/region"
FORCE_UCSC = os.environ.get("NTV3_FORCE_UCSC", "0") == "1"

SUPPORTED_UI_SPECIES: Dict[str, str] = {
    "human": "Human",
    "mouse": "Mouse",
    "rattus_norvegicus": "Rat",
    "canis_lupus_familiaris": "Dog",
    "felis_catus": "Cat",
    "gallus_gallus": "Chicken",
    "danio_rerio": "Zebrafish",
    "gorilla_gorilla": "Gorilla",
    "macaca_nemestrina": "Pig-tailed macaque",
    "bison_bison_bison": "Bison",
    "chinchilla_lanigera": "Chinchilla",
    "serinus_canaria": "Canary",
    "salmo_trutta": "Brown trout",
    "tetraodon_nigroviridis": "Green spotted puffer",
    "amphiprion_ocellaris": "Clownfish",
    "arabidopsis_thaliana": "Arabidopsis (Thale cress)",
    "oryza_sativa": "Rice",
    "glycine_max": "Soybean",
    "gossypium_hirsutum": "Cotton",
    "triticum_aestivum": "Wheat",
    "zea_mays": "Maize",
}

SPECIES_TO_ENSEMBL: Dict[str, str] = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "rattus_norvegicus": "rattus_norvegicus",
    "canis_lupus_familiaris": "canis_lupus_familiaris",
    "felis_catus": "felis_catus",
    "gallus_gallus": "gallus_gallus",
    "danio_rerio": "danio_rerio",
    "gorilla_gorilla": "gorilla_gorilla",
    "macaca_nemestrina": "macaca_nemestrina",
    "bison_bison_bison": "bison_bison_bison",
    "chinchilla_lanigera": "chinchilla_lanigera",
    "serinus_canaria": "serinus_canaria",
    "salmo_trutta": "salmo_trutta",
    "tetraodon_nigroviridis": "tetraodon_nigroviridis",
    "amphiprion_ocellaris": "amphiprion_ocellaris",
    "arabidopsis_thaliana": "arabidopsis_thaliana",
    "oryza_sativa": "oryza_sativa",
    "glycine_max": "glycine_max",
    "gossypium_hirsutum": "gossypium_hirsutum",
    "triticum_aestivum": "triticum_aestivum",
    "zea_mays": "zea_mays",
}

PLANT_SPECIES = {
    "arabidopsis_thaliana",
    "oryza_sativa",
    "glycine_max",
    "gossypium_hirsutum",
    "triticum_aestivum",
    "zea_mays",
}

# Global cache for model and genome
_MODEL_CACHE: Dict[str, Any] = {
    "model": None,
    "tokenizer": None,
    "genome": None,
    "bed_names": None,
    "bigwig_names": None,
    "selected_bw_indices": None,
    "nuc_token_map": None,
    "human_id": None,
    "species_to_token_id": None,
    "sequence_source": "ucsc",
}

# Cache for the most recent track profiles (used by tracks.py for plotting)
_LAST_TRACK_PROFILES: Dict[str, Any] = {}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def _normalize_ucsc_chrom(chrom: str) -> str:
    """Normalize chromosome string for UCSC API (expects chr-prefixed names)."""
    chrom = str(chrom).strip()
    if not chrom:
        return chrom
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _normalize_ensembl_chrom(chrom: str) -> str:
    """Normalize chromosome string for Ensembl REST (expects no chr prefix)."""
    clean = str(chrom).strip()
    if clean.lower().startswith("chr"):
        clean = clean[3:]
    if clean.upper() == "M":
        return "MT"
    return clean


@lru_cache(maxsize=2048)
def _fetch_ucsc_window(chrom: str, start: int, end: int) -> Optional[str]:
    """Fetch sequence window from UCSC API with in-session caching."""
    try:
        response = requests.get(
            UCSC_SEQUENCE_URL,
            params={
                "genome": "hg38",
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
            },
            timeout=10,
        )
        response.raise_for_status()
        payload = response.json()
        dna = payload.get("dna", "")
        if not dna:
            return None
        return str(dna).upper()
    except Exception:
        return None


@lru_cache(maxsize=2048)
def _fetch_ensembl_window(
    species: str, chrom: str, start: int, end: int
) -> Optional[str]:
    """Fetch sequence window from Ensembl REST for non-human species.

    For plant species, first tries Ensembl REST, then falls back to Ensembl Plants REST.
    """
    ensembl_species = SPECIES_TO_ENSEMBL.get(species)
    if not ensembl_species:
        return None

    clean_chrom = _normalize_ensembl_chrom(chrom)
    if not clean_chrom:
        return None

    start_1based = max(1, int(start) + 1)
    end_1based = max(start_1based, int(end))

    is_plant = species in PLANT_SPECIES
    urls = []

    if is_plant:
        urls.append(
            f"{ENSEMBL_PLANTS_SEQUENCE_URL}/{ensembl_species}/"
            f"{clean_chrom}:{start_1based}..{end_1based}"
        )
        urls.append(
            f"{ENSEMBL_SEQUENCE_URL}/{ensembl_species}/"
            f"{clean_chrom}:{start_1based}..{end_1based}"
        )
    else:
        urls.append(
            f"{ENSEMBL_SEQUENCE_URL}/{ensembl_species}/"
            f"{clean_chrom}:{start_1based}..{end_1based}"
        )

    for url in urls:
        try:
            response = requests.get(
                url,
                headers={"Content-Type": "text/plain", "Accept": "text/plain"},
                timeout=15,
            )
            response.raise_for_status()
            dna = response.text.strip()
            if dna:
                return dna.upper()
        except Exception:
            continue

    return None


def _fetch_sequence_from_local(
    genome: Fasta, chrom: str, start: int, end: int
) -> Optional[str]:
    """Fetch sequence window from local pyfaidx genome."""
    query_chrom = chrom
    if query_chrom not in genome:
        query_chrom = (
            query_chrom if query_chrom.startswith("chr") else f"chr{query_chrom}"
        )
    if query_chrom not in genome:
        query_chrom = (
            query_chrom.replace("chr", "")
            if "chr" in query_chrom
            else f"chr{query_chrom}"
        )
    if query_chrom not in genome:
        return None

    try:
        return str(genome[query_chrom][start:end]).upper()
    except Exception:
        return None


def get_genomic_sequence(
    genome: Optional[Fasta],
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    context_len: int = 4096,
    species: str = "human",
) -> Tuple[Optional[str], Optional[str], Optional[int]]:
    """
    Extract Ref/Alt sequences centered at variant position.

    Returns:
        (ref_seq, alt_seq, variant_center_idx) or (None, None, None) on error
    """
    variant_idx = pos - 1  # Convert to 0-based
    half = context_len // 2
    start = max(0, variant_idx - half)
    end = variant_idx + half

    ref_seq = None
    if species == "human":
        if genome is not None:
            ref_seq = _fetch_sequence_from_local(genome, str(chrom), start, end)

        if ref_seq is None:
            ucsc_chrom = _normalize_ucsc_chrom(str(chrom))
            ref_seq = _fetch_ucsc_window(ucsc_chrom, start, end)
            if ref_seq is None:
                warnings.warn(
                    f"Failed to fetch sequence for {chrom}:{pos} from local genome and UCSC API"
                )
                return None, None, None
    else:
        ref_seq = _fetch_ensembl_window(species, str(chrom), start, end)
        if ref_seq is None:
            warnings.warn(
                f"Failed to fetch sequence for {species} {chrom}:{pos} from Ensembl REST API"
            )
            return None, None, None

    if not ref_seq:
        return None, None, None

    center = variant_idx - start
    center = max(0, min(center, len(ref_seq) - 1))

    # Validate REF allele
    ref_end = min(center + len(ref), len(ref_seq))
    actual_ref = ref_seq[center:ref_end]
    if actual_ref.upper() != ref.upper():
        warnings.warn(
            f"REF mismatch at {chrom}:{pos}: expected '{ref}' but genome has '{actual_ref}'"
        )

    # Build ALT sequence
    alt_seq = ref_seq[:center] + alt + ref_seq[ref_end:]

    # Crop to equal length
    target_len = min(len(ref_seq), len(alt_seq))
    center = min(center, target_len - 1)

    return ref_seq[:target_len], alt_seq[:target_len], center


def get_track_indices(
    bigwig_names: List[str], metadata_file: Path
) -> Tuple[List[int], List[str]]:
    """Filter BigWig tracks to functional subset using metadata."""
    if not metadata_file.exists():
        print(f"⚠️  Metadata not found, using all {len(bigwig_names)} tracks")
        return list(range(len(bigwig_names))), bigwig_names

    metadata = pd.read_csv(metadata_file)

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

    print(f"📊 Selected {len(sel_idx)}/{len(bigwig_names)} BigWig tracks")
    return sel_idx, sel_names


def build_nuc_token_map(tokenizer):
    """Pre-compute nucleotide -> token ID mapping for LLR calculation."""
    return {
        nuc: tokenizer(nuc, add_special_tokens=False)["input_ids"][0]
        for nuc in "ACGTN"
        if tokenizer(nuc, add_special_tokens=False)["input_ids"]
    }


def to_track_probabilities(track_values):
    """
    Convert NTv3 track logits to probabilities via sigmoid.

    BED tracks: shape (B, L', 21, 2) → extract positive class [..., 1]
    BigWig tracks: shape (B, L', N) → apply sigmoid directly
    """
    if track_values is None:
        return None
    if track_values.shape[-1] == 2:  # Binary classification (BED)
        return torch.sigmoid(track_values[..., 1])
    return torch.sigmoid(track_values)


def compute_mlm_features(
    out_ref,
    out_alt,
    ref_seq: str,
    alt_seq: str,
    idx: int,
    variant_center: int,
    ref_allele: str,
    alt_allele: str,
    nuc_token_map: dict,
    use_kl: bool = True,
    use_embeddings: bool = False,
    window: int = 50,
) -> dict:
    """
    Compute MLM language model features (LLR, KL divergence, log-probs, embeddings).

    Returns dict with keys:
        - LLR, MLM_Prior, MLM_Delta (SNPs only)
        - MLM_KL_mean, MLM_KL_max
        - MLM_logprob_ref, MLM_logprob_alt, MLM_logprob_delta
        - REF_5mer, ALT_5mer
        - EMB_* (if use_embeddings=True)
    """
    feat = {}
    ref_logits = out_ref.logits[idx]
    alt_logits = out_alt.logits[idx]
    seq_len = min(ref_logits.shape[0], alt_logits.shape[0])

    # Determine variant span (1 for SNP, max allele length for indel)
    variant_span = max(1, len(ref_allele), len(alt_allele))
    variant_center = int(max(0, min(variant_center, seq_len - 1)))

    # KL window: cover exactly the variant span
    kl_ws = variant_center
    kl_we = min(seq_len, variant_center + variant_span)

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
        if len(seq) >= variant_center + 3:
            feat[f"{pf}_5mer"] = seq[max(0, variant_center - 2) : variant_center + 3]
        else:
            feat[f"{pf}_5mer"] = "NNNNN"

    # --- KL divergence + log-prob ---
    if use_kl:
        if kl_we <= kl_ws:
            feat["MLM_KL_mean"] = np.nan
            feat["MLM_KL_max"] = np.nan
            feat["MLM_logprob_ref"] = np.nan
            feat["MLM_logprob_alt"] = np.nan
            feat["MLM_logprob_delta"] = np.nan
            return feat

        rp_w = F.softmax(ref_logits[kl_ws:kl_we], dim=-1)
        ap_w = F.softmax(alt_logits[kl_ws:kl_we], dim=-1)
        kl = F.kl_div(rp_w.log(), ap_w, reduction="none", log_target=False).sum(-1)
        feat["MLM_KL_mean"] = float(kl.mean().cpu())
        feat["MLM_KL_max"] = float(kl.max().cpu())

        # Log-probs
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

    # --- Embedding distances ---
    if use_embeddings:
        hr = getattr(out_ref, "last_hidden_state", None)
        ha = getattr(out_alt, "last_hidden_state", None)
        if hr is not None and ha is not None:
            ws = max(0, variant_center - window)
            we = min(seq_len, variant_center + window)
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
# MODEL LOADING AND CACHING
# ============================================================================
def load_model_and_resources(device: str = "cuda"):
    """
    Load model, tokenizer, genome, and metadata once at startup.
    Caches everything in global _MODEL_CACHE.

    Handles:
    - HF_TOKEN authentication for gated models
    - bf16 mixed precision when supported by GPU
    """
    if _MODEL_CACHE["model"] is not None:
        return  # Already loaded

    # Get HF token from environment (needed for gated models like NTv3)
    hf_token = os.environ.get("HF_TOKEN")
    if not hf_token:
        print("⚠️  HF_TOKEN not set; model loading may fail if the model is gated")

    print(f"🧠 Loading NTv3 model '{MODEL_NAME}' on {device}...")

    # Prepare bf16 kwargs if GPU supports bfloat16
    bf16_kwargs = {}
    if device == "cuda" and torch.cuda.is_bf16_supported():
        print("💾 Enabling bf16 mixed precision to reduce memory usage")
        bf16_kwargs = {
            "stem_compute_dtype": "bfloat16",
            "down_convolution_compute_dtype": "bfloat16",
            "transformer_qkvo_compute_dtype": "bfloat16",
            "transformer_ffn_compute_dtype": "bfloat16",
            "up_convolution_compute_dtype": "bfloat16",
            "modulation_compute_dtype": "bfloat16",
        }

    try:
        model = (
            AutoModel.from_pretrained(
                MODEL_NAME, trust_remote_code=True, token=hf_token, **bf16_kwargs
            )
            .to(device)
            .eval()
        )
        tokenizer = AutoTokenizer.from_pretrained(
            MODEL_NAME, trust_remote_code=True, token=hf_token
        )
    except Exception as e:
        error_msg = f"Failed to load model: {str(e)}"
        if "gated" in str(e).lower() or "401" in str(e):
            error_msg += "\n\nThe NTv3_650M_post model is gated. You need to:\n1. Accept the model terms at https://huggingface.co/InstaDeepAI/NTv3_650M_post\n2. Set HF_TOKEN as an environment variable with your HuggingFace API token"
        raise RuntimeError(error_msg) from e

    species_to_token_id = dict(getattr(model.config, "species_to_token_id", {}) or {})

    # Extract human species ID
    try:
        human_id = model.encode_species(["human"]).item()
    except (AttributeError, TypeError):
        human_id = species_to_token_id.get("human", 6)

    # Extract BED/BigWig names
    bed_names = []
    if USE_BED:
        for attr in ("bed_elements_names", "bed_tracks", "bed_track_labels"):
            if hasattr(model.config, attr):
                bed_names = getattr(model.config, attr)
                break
    bw_names = []
    if USE_BIGWIGS and hasattr(model.config, "bigwigs_per_species"):
        bw_names = model.config.bigwigs_per_species.get("human", [])

    # Filter BigWig tracks
    selected_bw_indices, selected_bw_names = get_track_indices(bw_names, METADATA_FILE)

    # Load genome optionally (for fast local lookups); otherwise UCSC fallback
    genome = None
    sequence_source = "ucsc"
    if FORCE_UCSC:
        print("⚠️  NTV3_FORCE_UCSC=1 set; using UCSC API for sequence retrieval")
    elif GENOME_FILE.exists():
        print(f"🧬 Loading local reference genome from {GENOME_FILE}...")
        genome = Fasta(str(GENOME_FILE))
        sequence_source = "local+ucsc-fallback"
    elif GENOME_GZ_FILE.exists():
        print(
            f"⚠️  Found compressed genome at {GENOME_GZ_FILE}; using UCSC API (decompress to enable local fast path)"
        )
    else:
        print("⚠️  No local hg38.fa found; using UCSC API for sequence retrieval")

    # Build token map
    nuc_token_map = build_nuc_token_map(tokenizer)

    # Cache everything
    _MODEL_CACHE.update(
        {
            "model": model,
            "tokenizer": tokenizer,
            "genome": genome,
            "bed_names": bed_names,
            "bigwig_names": bw_names,
            "selected_bw_indices": selected_bw_indices,
            "nuc_token_map": nuc_token_map,
            "human_id": human_id,
            "species_to_token_id": species_to_token_id,
            "sequence_source": sequence_source,
        }
    )

    print(
        f"✅ Model loaded: {len(bed_names)} BED, {len(selected_bw_indices)} BigWig tracks | sequence source: {sequence_source}"
    )


# ============================================================================
# INFERENCE
# ============================================================================
def predict_variants(
    df: pd.DataFrame,
    device: str = "cuda",
    species: str = "human",
    cache_profiles: bool = True,
) -> pd.DataFrame:
    """
    Run NTv3 inference on variants DataFrame.

    Args:
        df: DataFrame with columns ['chrom', 'pos', 'ref', 'alt']
        device: 'cuda' or 'cpu'

    Returns:
        DataFrame with original columns plus:
        - D_BED_* (21 BED element deltas)
        - REF_BED_* (21 BED element ref probabilities)
        - D_BW_* (filtered BigWig deltas)
        - REF_BW_* (filtered BigWig ref probabilities)
        - LLR, MLM_Prior, MLM_Delta, MLM_KL_mean, MLM_KL_max
        - MLM_logprob_ref, MLM_logprob_alt, MLM_logprob_delta
        - REF_5mer, ALT_5mer
        - EMB_* (if embeddings enabled)
        - indel_size
    """
    # Load model if not already loaded
    if _MODEL_CACHE["model"] is None:
        load_model_and_resources(device)

    model = _MODEL_CACHE.get("model")
    tokenizer = _MODEL_CACHE.get("tokenizer")
    genome = cast(Optional[Fasta], _MODEL_CACHE.get("genome"))
    bed_names = cast(List[str], _MODEL_CACHE.get("bed_names") or [])
    bigwig_names = cast(List[str], _MODEL_CACHE.get("bigwig_names") or [])
    selected_bw_indices = cast(List[int], _MODEL_CACHE.get("selected_bw_indices") or [])
    nuc_token_map = cast(dict, _MODEL_CACHE.get("nuc_token_map") or {})
    species_to_token_id = cast(
        Dict[str, int], _MODEL_CACHE.get("species_to_token_id") or {}
    )

    species = str(species).strip()
    species_id = species_to_token_id.get(species)
    active_bigwig_names = bigwig_names if species == "human" else []
    active_bw_indices = selected_bw_indices if species == "human" else []

    if model is None or tokenizer is None or species_id is None:
        raise RuntimeError("Model resources are not initialized correctly")

    # Validate input
    required_cols = {"chrom", "pos", "ref", "alt"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Clean data
    df = df.copy()
    df = df[df["ref"].notna() & df["alt"].notna()].reset_index(drop=True)

    results = []

    for idx, row in df.iterrows():
        # Get sequences
        ref_seq, alt_seq, vcenter = get_genomic_sequence(
            genome,
            row["chrom"],
            row["pos"],
            row["ref"],
            row["alt"],
            CONTEXT_LEN,
            species=species,
        )

        if ref_seq is None or alt_seq is None or vcenter is None:
            # Failed to fetch sequence - return NaN results
            res = {c: row[c] for c in df.columns}
            res["indel_size"] = len(str(row["alt"])) - len(str(row["ref"]))
            for nm in bed_names:
                res[f"REF_BED_{nm}"] = np.nan
                res[f"D_BED_{nm}"] = np.nan
            for gi in active_bw_indices:
                res[f"REF_BW_{active_bigwig_names[gi]}"] = np.nan
                res[f"D_BW_{active_bigwig_names[gi]}"] = np.nan
            for k in (
                "LLR",
                "MLM_Prior",
                "MLM_Delta",
                "MLM_KL_mean",
                "MLM_KL_max",
                "MLM_logprob_ref",
                "MLM_logprob_alt",
                "MLM_logprob_delta",
                "REF_5mer",
                "ALT_5mer",
            ):
                res[k] = np.nan if k not in ("REF_5mer", "ALT_5mer") else "NNNNN"
            results.append(res)
            continue

        # Tokenize
        tok_kw = dict(
            return_tensors="pt",
            padding="max_length",
            max_length=CONTEXT_LEN,
            truncation=True,
            add_special_tokens=False,
            pad_to_multiple_of=128,
        )
        inp_r = tokenizer([ref_seq], **tok_kw).to(device)
        inp_a = tokenizer([alt_seq], **tok_kw).to(device)

        # Forward pass
        with torch.no_grad():
            sp = torch.tensor([species_id], device=device)
            if USE_EMBEDDINGS:
                try:
                    out_r = model(**inp_r, species_ids=sp, output_hidden_states=True)
                    out_a = model(**inp_a, species_ids=sp, output_hidden_states=True)
                    for o in (out_r, out_a):
                        if getattr(o, "last_hidden_state", None) is None and hasattr(
                            o, "hidden_states"
                        ):
                            o.last_hidden_state = o.hidden_states[-1]
                except TypeError:
                    out_r = model(**inp_r, species_ids=sp)
                    out_a = model(**inp_a, species_ids=sp)
            else:
                out_r = model(**inp_r, species_ids=sp)
                out_a = model(**inp_a, species_ids=sp)

        # Build result
        res = {c: row[c] for c in df.columns}
        ref_allele = str(row["ref"])
        alt_allele = str(row["alt"])
        res["indel_size"] = len(alt_allele) - len(ref_allele)
        variant_span = max(1, len(ref_allele), len(alt_allele))
        in_len = int(inp_r["input_ids"].shape[1])

        # === BED tracks ===
        bed_r = getattr(out_r, "bed_tracks_logits", None)
        bed_a = getattr(out_a, "bed_tracks_logits", None)
        if USE_BED and bed_r is not None and bed_a is not None:
            bed_r_probs = to_track_probabilities(bed_r[0])
            bed_a_probs = to_track_probabilities(bed_a[0])
            track_len = int(bed_r_probs.shape[0])
            track_start = max(0, (in_len - track_len) // 2)
            bed_pos = vcenter - track_start
            if 0 <= bed_pos < track_len:
                be = min(bed_pos + variant_span, track_len)
                br = bed_r_probs[bed_pos:be].mean(0).float().cpu().numpy()
                ba = bed_a_probs[bed_pos:be].mean(0).float().cpu().numpy()
                for j, nm in enumerate(bed_names):
                    res[f"REF_BED_{nm}"] = float(br[j])
                    res[f"D_BED_{nm}"] = float(ba[j] - br[j])
            else:
                for nm in bed_names:
                    res[f"REF_BED_{nm}"] = np.nan
                    res[f"D_BED_{nm}"] = np.nan
        else:
            for nm in bed_names:
                res[f"REF_BED_{nm}"] = np.nan
                res[f"D_BED_{nm}"] = np.nan

        # === BigWig tracks ===
        bw_r = getattr(out_r, "bigwig_tracks_logits", None)
        bw_a = getattr(out_a, "bigwig_tracks_logits", None)
        if species == "human" and USE_BIGWIGS and bw_r is not None and bw_a is not None:
            bw_r_probs = to_track_probabilities(bw_r[0])
            bw_a_probs = to_track_probabilities(bw_a[0])
            track_len = int(bw_r_probs.shape[0])
            track_start = max(0, (in_len - track_len) // 2)
            bw_pos = vcenter - track_start
            if 0 <= bw_pos < track_len:
                bwe = min(bw_pos + variant_span, track_len)
                bwr = bw_r_probs[bw_pos:bwe].mean(0).float().cpu().numpy()
                bwa = bw_a_probs[bw_pos:bwe].mean(0).float().cpu().numpy()
                for gi in active_bw_indices:
                    res[f"REF_BW_{active_bigwig_names[gi]}"] = float(bwr[gi])
                    res[f"D_BW_{active_bigwig_names[gi]}"] = float(bwa[gi] - bwr[gi])
            else:
                for gi in active_bw_indices:
                    res[f"REF_BW_{active_bigwig_names[gi]}"] = np.nan
                    res[f"D_BW_{active_bigwig_names[gi]}"] = np.nan
        else:
            for gi in active_bw_indices:
                res[f"REF_BW_{active_bigwig_names[gi]}"] = np.nan
                res[f"D_BW_{active_bigwig_names[gi]}"] = np.nan

        # === MLM features ===
        res.update(
            compute_mlm_features(
                out_r,
                out_a,
                ref_seq,
                alt_seq,
                0,
                vcenter,
                ref_allele,
                alt_allele,
                nuc_token_map,
                use_kl=USE_KL_DIVERGENCE,
                use_embeddings=USE_EMBEDDINGS,
                window=MLM_WINDOW,
            )
        )

        # === Cache full track profiles for plotting (single-variant mode only) ===
        if cache_profiles:
            _cache_track_profiles(
                row,
                vcenter,
                in_len,
                bed_names,
                active_bigwig_names,
                active_bw_indices,
                bed_r,
                bed_a,
                bw_r,
                bw_a,
        )

        results.append(res)

    return pd.DataFrame(results)


def _cache_track_profiles(
    row,
    vcenter,
    in_len,
    bed_names,
    bigwig_names,
    selected_bw_indices,
    bed_r_logits,
    bed_a_logits,
    bw_r_logits,
    bw_a_logits,
):
    """Cache the most recent variant's full track-length logit profiles for plotting."""
    global _LAST_TRACK_PROFILES

    profiles: Dict[str, Any] = {
        "chrom": str(row["chrom"]),
        "pos": int(row["pos"]),
        "ref": str(row["ref"]),
        "alt": str(row["alt"]),
        "variant_center": vcenter,
        "input_len": in_len,
        "bed_names": list(bed_names),
        "bigwig_names": list(bigwig_names),
        "selected_bw_indices": list(selected_bw_indices),
    }

    # BED profiles: convert to probabilities and store as numpy (L, 21)
    if bed_r_logits is not None and bed_a_logits is not None:
        bed_ref_probs = to_track_probabilities(bed_r_logits[0]).float().cpu().numpy()
        bed_alt_probs = to_track_probabilities(bed_a_logits[0]).float().cpu().numpy()
        profiles["bed_ref"] = bed_ref_probs
        profiles["bed_alt"] = bed_alt_probs
        profiles["bed_track_len"] = bed_ref_probs.shape[0]
        profiles["bed_track_start"] = max(0, (in_len - bed_ref_probs.shape[0]) // 2)
    else:
        profiles["bed_ref"] = profiles["bed_alt"] = None

    # BigWig profiles: convert to probabilities and store as numpy (L, T)
    if bw_r_logits is not None and bw_a_logits is not None:
        bw_ref_probs = to_track_probabilities(bw_r_logits[0]).float().cpu().numpy()
        bw_alt_probs = to_track_probabilities(bw_a_logits[0]).float().cpu().numpy()
        profiles["bw_ref"] = bw_ref_probs
        profiles["bw_alt"] = bw_alt_probs
        profiles["bw_track_len"] = bw_ref_probs.shape[0]
        profiles["bw_track_start"] = max(0, (in_len - bw_ref_probs.shape[0]) // 2)
    else:
        profiles["bw_ref"] = profiles["bw_alt"] = None

    _LAST_TRACK_PROFILES = profiles
