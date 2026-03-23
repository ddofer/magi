#!/usr/bin/env python3
"""
OMIA Multi-Species Variant Inference Pipeline
==============================================

Purpose
  Fetch single-gene Mendelian variants from OMIA (Online Mendelian Inheritance
  in Animals), wrangle them into the same format as ClinVar SNP/INDEL inputs,
  and run NTv3 inference for supported non-human species.

    Acts as a counterpart to inference.py (human ClinVar) and
    inference_indel.py (human INDELs), extending the pipeline to animals.

Supported Species (NTv3 model config species_to_token_id):
  dog       canis_lupus_familiaris  (token 17)
  cat       felis_catus             (token  9)
  chicken   gallus_gallus           (token 23)
  zebrafish danio_rerio             (token 10)
  rat       rattus_norvegicus       (token 28)
  mouse     mouse                   (token 29)

  Cattle, horses, pigs and sheep are NOT in the NTv3 model config and are
  skipped with a warning.

Usage
    python inference_omia.py --validate-only
        # validation-only: fetch + parse + coordinate checks, no model inference

    python inference_omia.py
        # sample inference: default species, 30 variants/species
    python inference_omia.py --species dog cat --sample-size 50
        # sample inference: selected species with custom sample size

    CUDA_VISIBLE_DEVICES=2,3 python inference_omia.py --full
        # full inference (all variants for run species) on GPUs 2 and 3

  python inference_omia.py --species dog cat     # specific species
  python inference_omia.py --fetch-only          # download data + genomes
  python inference_omia.py --model-size 100M     # smaller model

Data Sources
  OMIA results table (cached locally for reproducibility):
    https://www.omia.org/results/?search_type=advanced&result_type=variant&singlelocus=yes&limit=10000
  Reference genomes matched to OMIA assemblies (Ensembl + Ensembl Rapid Release)

Output Format
        Parquet files in data/ and data/results/.
        Metadata columns include:
        chrom, pos, ref, alt, label, rationale, species, gene, omia_id,
        label_source, variant_phenotype, variant_effect, variant_type
        Feature columns include REF_BED_*, D_BED_* (42 total), LLR, and MLM_* columns.
        Full output is typically 66 columns:
            13 metadata + 42 BED + 1 LLR + 7 MLM + 3 convenience fields.
        Note: BigWig columns are NOT produced for non-human species.
            The NTv3 MultiSpeciesHead uses ZeroHead for species without
            BigWig training data, yielding constant sigmoid(0)=0.5 outputs.

Gotchas
  1. NTv3 BED tracks (21 genomic elements) work for all species (shared head).
     BigWig tracks (7k) are **human-only by architecture**: the MultiSpeciesHead
     uses a separate LinearHead per species, and non-human species get ZeroHead
     (returns 0-dim tensors, padded to max tracks → logit=0 → sigmoid=0.5).
     Ref and alt both get 0.5, so deltas are always 0.0.  USE_BIGWIGS=False.
  2. Ensembl chroms use NO 'chr' prefix (1, 2, X, A1, etc.).
     OMIA data sometimes includes 'chr' — we strip it.
  3. Mouse NTv3 token is 'mouse' (not 'mus_musculus').
  4. Chromosome validation uses samtools faidx for spot-checks.
  5. OMIA 'Allele' column is often blank — extract from HGVS g. string.
  6. Only SNPs are supported for LLR feature; INDELs get KL/logprob features.
  7. Many REF mismatch warnings with expected N or '-' are normal for INDEL
      placeholders; mismatches where expected REF is A/C/G/T should be reviewed.

Authors: Dan Ofer / Michal Linial Lab (based on inference.py)
"""

import argparse
import io
import os
import re
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests

import inference as base

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths (all OMIA artifacts are kept under the project-level data/ folder)
DATA_DIR = Path("data")
GENOME_DIR = DATA_DIR / "genomes"
RESULTS_DIR = DATA_DIR / "results"
OMIA_CACHE_FILE = DATA_DIR / "omia_raw_variants.csv"
OMIA_URL = (
    "https://www.omia.org/results/?search_type=advanced"
    "&result_type=variant&singlelocus=yes&limit=10000"
)
DEFAULT_RUN_SPECIES = ["dog", "cat", "chicken"]

# Inference defaults — override via CLI
SAMPLE_SIZE = 30  # variants per species for smoke test (None = all)
MODEL_SIZE = "650M"  # '100M' or '650M'
CONTEXT_LEN = 16 * 1024  # smaller than human (64K) for speed in animal runs
BATCH_SIZE = 2
USE_BED = True
USE_BIGWIGS = False  # BigWig head uses species-specific LinearHead/ZeroHead;
# non-human species get ZeroHead → sigmoid(0)=0.5 constant
# outputs, making deltas always 0.  Human-only by architecture.
USE_KL_DIVERGENCE = True
USE_EMBEDDINGS = False  # set True for INDEL-heavy species
METADATA_FILE = (
    "functional_tracks_metadata_human.csv"  # for BigWig subset filter (human only)
)

# ============================================================================
# SPECIES CONFIG
# ============================================================================
# All NTv3-supported animal species that plausibly appear in OMIA.
# Keys are short names used as CLI args and filename stems.
# 'token' must exactly match model config species_to_token_id keys.
# 'omia_names': substrings to match against OMIA "Species Name" column
#   (case-insensitive); first match wins.
# 'unsupported': if True, skip with warning instead of running.

SPECIES_CONFIG = {
    "dog": {
        "token": "canis_lupus_familiaris",
        "omia_names": ["canis lupus familiaris", "dog", "canis familiaris"],
        "common_name": "Dog",
        "target_assemblies": ["canfam3.1"],
        "genome_url": (
            "https://ftp.ensembl.org/pub/release-104/fasta/"
            "canis_lupus_familiaris/dna/"
            "Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa.gz"
        ),
        "genome_fa": "Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa",
        "ucsc_build": "canFam3",
    },
    "cat": {
        "token": "felis_catus",
        "omia_names": ["felis catus", "cat", "felis silvestris catus"],
        "common_name": "Cat",
        "target_assemblies": ["f.catus_fca126_mat1.0", "f.catus_fcat126_mat1.0"],
        "genome_url": (
            "https://ftp.ensembl.org/pub/rapid-release/species/"
            "Felis_catus/GCA_018350175.1/ensembl/genome/"
            "Felis_catus-GCA_018350175.1-unmasked.fa.gz"
        ),
        "genome_fa": "Felis_catus-GCA_018350175.1-unmasked.fa",
        "ucsc_build": "felCat9",
    },
    "chicken": {
        "token": "gallus_gallus",
        "omia_names": ["gallus gallus", "chicken", "red junglefowl"],
        "common_name": "Chicken",
        "target_assemblies": ["grcg6a"],
        "genome_url": (
            "https://ftp.ensembl.org/pub/release-104/fasta/"
            "gallus_gallus/dna/"
            "Gallus_gallus.GRCg6a.dna.toplevel.fa.gz"
        ),
        "genome_fa": "Gallus_gallus.GRCg6a.dna.toplevel.fa",
        "ucsc_build": "galGal6",
    },
    "zebrafish": {
        "token": "danio_rerio",
        "omia_names": ["danio rerio", "zebrafish"],
        "common_name": "Zebrafish",
        "genome_url": (
            "https://ftp.ensembl.org/pub/release-113/fasta/"
            "danio_rerio/dna/"
            "Danio_rerio.GRCz11.dna.toplevel.fa.gz"
        ),
        "genome_fa": "Danio_rerio.GRCz11.dna.toplevel.fa",
        "ucsc_build": "danRer11",
    },
    "rat": {
        "token": "rattus_norvegicus",
        "omia_names": ["rattus norvegicus", "norway rat", "rat"],
        "common_name": "Rat",
        "genome_url": (
            "https://ftp.ensembl.org/pub/release-113/fasta/"
            "rattus_norvegicus/dna/"
            "Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
        ),
        "genome_fa": "Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa",
        "ucsc_build": "rn7",
    },
    "mouse": {
        "token": "mouse",  # NTv3 token is literally 'mouse', not 'mus_musculus'
        "omia_names": ["mus musculus", "mouse", "house mouse"],
        "common_name": "Mouse",
        "genome_url": (
            "https://ftp.ensembl.org/pub/release-113/fasta/"
            "mus_musculus/dna/"
            "Mus_musculus.GRCm39.dna.toplevel.fa.gz"
        ),
        "genome_fa": "Mus_musculus.GRCm39.dna.toplevel.fa",
        "ucsc_build": "mm39",
    },
    # ---------- species NOT in NTv3 model config → skipped with warning ----------
    "cattle": {
        "token": None,
        "omia_names": ["bos taurus", "cattle", "cow", "bovine"],
        "common_name": "Cattle",
        "unsupported": True,
    },
    "horse": {
        "token": None,
        "omia_names": ["equus caballus", "horse"],
        "common_name": "Horse",
        "unsupported": True,
    },
    "pig": {
        "token": None,
        "omia_names": ["sus scrofa", "pig", "swine"],
        "common_name": "Pig",
        "unsupported": True,
    },
    "sheep": {
        "token": None,
        "omia_names": ["ovis aries", "sheep"],
        "common_name": "Sheep",
        "unsupported": True,
    },
}

# Pre-compiled HGVS → chrom/pos/ref/alt patterns
# For RefSeq HGVS: NC_058375.1:g.147442389C>T
_HGVS_G_FULL = re.compile(r"g\.(\d+)([ACGT])>([ACGT])", re.IGNORECASE)
# HGVS insertion: g.147442389insACGT
_HGVS_G_INS = re.compile(r"g\.(\d+)ins([ACGT]+)", re.IGNORECASE)
# HGVS deletion: g.147442389del or g.147442389_147442391del
_HGVS_G_DEL = re.compile(r"g\.(\d+)(?:_\d+)?del([ACGT]*)", re.IGNORECASE)
# HGVS duplication: g.147442389dup or g.147442389dupC
_HGVS_G_DUP = re.compile(r"g\.(\d+)(?:_\d+)?dup([ACGT]*)", re.IGNORECASE)
# HGVS delins: g.147442389_147442390delinsAC
_HGVS_G_DELINS = re.compile(r"g\.(\d+)(?:_\d+)?delins([ACGT]+)", re.IGNORECASE)
# Generic SNP fallback: G>A or G/A
_GENERIC_SNP = re.compile(r"\b([ACGT])>([ACGT])\b", re.IGNORECASE)
# Numeric position standalone (last resort)
_POS_NUM = re.compile(r"\b(\d{4,})\b")


# ============================================================================
# DIRECTORY SETUP
# ============================================================================


def ensure_dirs():
    """Create required output directories."""
    DATA_DIR.mkdir(exist_ok=True)
    GENOME_DIR.mkdir(exist_ok=True)
    RESULTS_DIR.mkdir(exist_ok=True)


# ============================================================================
# GENOME MANAGEMENT
# ============================================================================


def _genome_path(species_key: str) -> Path:
    """Return path to the uncompressed FASTA for a species."""
    cfg = SPECIES_CONFIG[species_key]
    return GENOME_DIR / species_key / cfg["genome_fa"]


def download_genome(species_key: str) -> Optional[Path]:
    """
    Download, decompress, and index reference genome from Ensembl FTP.

    Uses wget -nc (no-clobber) so re-runs are safe.
    Requires samtools and gunzip on PATH.

    Returns path to the .fa file, or None if download fails.
    """
    cfg = SPECIES_CONFIG[species_key]
    if cfg.get("unsupported"):
        return None

    species_genome_dir = GENOME_DIR / species_key
    species_genome_dir.mkdir(parents=True, exist_ok=True)

    fa_path = species_genome_dir / cfg["genome_fa"]
    gz_path = species_genome_dir / (cfg["genome_fa"] + ".gz")
    fai_path = Path(str(fa_path) + ".fai")

    if fa_path.exists() and fai_path.exists():
        print(f"   ✅ Genome already indexed: {fa_path}")
        return fa_path

    # Download
    if not gz_path.exists() and not fa_path.exists():
        print(f"   📥 Downloading {cfg['common_name']} genome (~1-3 GB)...")
        print(f"      URL: {cfg['genome_url']}")
        result = subprocess.run(
            [
                "wget",
                "-nc",
                "-q",
                "--show-progress",
                "-P",
                str(species_genome_dir),
                cfg["genome_url"],
            ],
            check=False,
        )
        if result.returncode != 0:
            print(
                f"   ❌ wget failed for {species_key}. "
                f"Try manually: wget -c {cfg['genome_url']}"
            )
            return None

    # Decompress
    if not fa_path.exists() and gz_path.exists():
        print(f"   📦 Decompressing {gz_path.name}...")
        result = subprocess.run(["gunzip", "-k", "-f", str(gz_path)], check=False)
        if result.returncode != 0:
            print(f"   ❌ gunzip failed for {gz_path}")
            return None

    if not fa_path.exists():
        print(f"   ❌ FASTA not found after download: {fa_path}")
        return None

    # Index with samtools
    if not fai_path.exists():
        print("   🔍 Indexing with samtools faidx...")
        result = subprocess.run(["samtools", "faidx", str(fa_path)], check=False)
        if result.returncode != 0:
            print("   ⚠️  samtools faidx failed — coordinate validation will be skipped")

    print(f"   ✅ Genome ready: {fa_path}")
    return fa_path


# ============================================================================
# OMIA DATA FETCHING
# ============================================================================


def fetch_raw_omia(force_refresh: bool = False) -> pd.DataFrame:
    """
    Fetch variant data from OMIA HTML table with local caching.

    OMIA CSV endpoints are currently unavailable (404); HTML table is used as
    canonical source and cached locally for reproducibility.
    """
    if OMIA_CACHE_FILE.exists() and not force_refresh:
        print(f"   📂 Loading cached OMIA data from {OMIA_CACHE_FILE}")
        cached = pd.read_csv(OMIA_CACHE_FILE, low_memory=False)
        marker_cols = [
            c
            for c in cached.columns
            if str(c).lower() in ("verbal description", "omia id")
        ]
        is_synthetic = False
        for col in marker_cols:
            vals = cached[col].astype(str).str.lower()
            if vals.str.contains("synthetic omia-style cached row").any():
                is_synthetic = True
                break
            if str(col).lower() == "omia id" and vals.str.startswith("syn").any():
                is_synthetic = True
                break
        if not is_synthetic:
            return cached
        print(
            "   ⚠️  Cached OMIA file appears synthetic/placeholder; "
            "refetching live data."
        )

    print(f"   🌐 Fetching OMIA variant data from {OMIA_URL} ...")
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (compatible; NTv3-OMIA-script/1.0; "
            "+https://github.com/instadeep)"
        )
    }

    try:
        r = requests.get(OMIA_URL, headers=headers, timeout=120)
        r.raise_for_status()
        tables = pd.read_html(io.StringIO(r.text))
        for tbl in tables:
            cols = [str(c).lower() for c in tbl.columns]
            if any("phenotype" in c or "variant" in c or "gene" in c for c in cols):
                print(
                    f"   ✅ Extracted HTML table with {len(tbl)} rows, "
                    f"{len(tbl.columns)} cols"
                )
                tbl.to_csv(OMIA_CACHE_FILE, index=False)
                return tbl
        raise ValueError("No suitable table found in OMIA HTML response")
    except Exception as e:
        raise ConnectionError(
            f"Failed to fetch OMIA data from HTML endpoint: {e}"
        ) from e


# ============================================================================
# VARIANT PARSING
# ============================================================================


def _normalise_columns(df: pd.DataFrame) -> dict:
    """Return a lower-stripped col-name → original-col-name mapping."""
    return {str(c).lower().strip(): c for c in df.columns}


def _extract_snp_from_hgvs(g_str: str):
    """
    Try to extract (pos, ref, alt) from an HGVS genomic string.

    Priority:
      1. HGVS substitution:  g.123456C>T
      2. Returns (None, None, None) if not a simple SNP substitution.
    """
    m = _HGVS_G_FULL.search(g_str)
    if m:
        return int(m.group(1)), m.group(2).upper(), m.group(3).upper()
    return None, None, None


def _extract_indel_from_hgvs(g_str: str):
    """
    Try to extract (pos, ref, alt) for simple INDELs from HGVS.

    Returns (pos, ref, alt) or (None, None, None).
      - insertion:  pos, "-", inserted_seq
      - deletion:   pos, deleted_seq_or_"N", "-"
      - duplication: pos, "-", duplicated_seq_or_"N"
      - delins: pos, "N", inserted_seq
    """
    m_ins = _HGVS_G_INS.search(g_str)
    if m_ins:
        return int(m_ins.group(1)), "-", m_ins.group(2).upper()
    m_dup = _HGVS_G_DUP.search(g_str)
    if m_dup:
        dup_seq = m_dup.group(2).upper() if m_dup.group(2) else "N"
        return int(m_dup.group(1)), "-", dup_seq
    m_del = _HGVS_G_DEL.search(g_str)
    if m_del:
        deleted = m_del.group(2).upper() if m_del.group(2) else "N"
        return int(m_del.group(1)), deleted, "-"
    m_delins = _HGVS_G_DELINS.search(g_str)
    if m_delins:
        return int(m_delins.group(1)), "N", m_delins.group(2).upper()
    return None, None, None


def _is_valid_allele(allele: str) -> bool:
    """Accept only nucleotide strings (plus '-' for INDEL placeholder)."""
    if not allele or allele in ("nan", "None", "."):
        return False
    if allele == "-":
        return True
    return bool(re.match(r"^[ACGTN\-]+$", allele.upper()))


def _norm_text(value: str) -> str:
    return re.sub(r"\s+", " ", str(value).strip().lower())


def _species_match(species_value: str, patterns: list[str]) -> bool:
    text = _norm_text(species_value)
    for pat in patterns:
        p = _norm_text(pat)
        if " " in p:
            if p in text:
                return True
        else:
            if re.search(rf"\b{re.escape(p)}\b", text):
                return True
    return False


def _assembly_match(assembly_value: str, target_patterns: list[str]) -> bool:
    text = _norm_text(assembly_value)
    return any(_norm_text(pat) in text for pat in target_patterns)


def parse_omia_variants(
    df: pd.DataFrame,
    species_key: str,
    include_indels: bool = True,
) -> pd.DataFrame:
    """
    Filter OMIA dataframe to one species and parse into NTv3-compatible format.

    Columns produced:
      chrom, pos, ref, alt, label, rationale, species, gene,
      omia_id, variant_phenotype, variant_effect, variant_type

    HGVS extraction priority:
      1. 'g. or m.' column (genomic HGVS) — most reliable
      2. 'c. or n.' column + fallback generic SNP pattern + position column
    """
    cfg = SPECIES_CONFIG[species_key]
    cmap = _normalise_columns(df)

    # Filter to this species using strict matching on Species Name column
    sp_col = cmap.get("species name", cmap.get("species", None))
    if sp_col is None:
        for k, v in cmap.items():
            if "species" in k:
                sp_col = v
                break
    if sp_col is None:
        print("   ⚠️  No species column found in OMIA data; trying all rows")
        species_df = df.copy()
    else:
        omia_name_patterns = cfg["omia_names"]
        mask = (
            df[sp_col]
            .astype(str)
            .apply(lambda s: _species_match(s, omia_name_patterns))
        )
        species_df = df[mask].copy()

    if species_df.empty:
        print(f"   ⚠️  No OMIA rows found for {cfg['common_name']}")
        return pd.DataFrame()

    print(f"   📊 {cfg['common_name']}: {len(species_df)} OMIA variant rows")

    # Column aliases
    cmap2 = _normalise_columns(species_df)

    def gcol(key, *aliases):
        for k in [key, *aliases]:
            if k in cmap2:
                return cmap2[k]
        return None

    chrom_col = gcol("chromosome", "chr.", "chr", "chrom")
    gm_col = gcol("g. or m.", "g.or m.", "hgvs genomic", "genomic hgvs")
    cn_col = gcol("c. or n.", "c.or n.", "hgvs coding", "hgvs")
    pos_col = gcol("position", "pos", "start")
    gene_col = gcol("gene", "gene symbol", "gene_symbol")
    phene_col = gcol("variant phenotype", "phenotype", "phene", "trait")
    effect_col = gcol("variant effect", "variant_effect", "effect")
    desc_col = gcol(
        "verbal description", "description", "verbal_description", "summary"
    )
    omia_id_col = gcol("omia id", "omia_id", "id", "omia-id")
    allele_col = gcol("allele", "alleles", "allele1", "variant allele")
    refseq_col = gcol("reference sequence", "reference", "assembly")

    target_assemblies = cfg.get("target_assemblies", [])
    if refseq_col and target_assemblies:
        before = len(species_df)
        species_df = species_df[
            species_df[refseq_col]
            .astype(str)
            .apply(lambda x: _assembly_match(x, target_assemblies))
        ].copy()
        print(
            f"   🧭 Assembly filter ({target_assemblies}): "
            f"{before} → {len(species_df)} rows"
        )
        if species_df.empty:
            print(
                f"   ⚠️  No rows for {cfg['common_name']} matched target assembly "
                f"{target_assemblies}"
            )
            return pd.DataFrame()

    records = []

    for _, row in species_df.iterrows():
        chrom = str(row[chrom_col]).strip() if chrom_col else ""
        g_str = str(row[gm_col]).strip() if gm_col else ""
        c_str = str(row[cn_col]).strip() if cn_col else ""
        allele_str = str(row[allele_col]).strip() if allele_col else ""

        ref, alt, pos = None, None, None

        # --- Priority 1: HGVS genomic substitution ---
        pos, ref, alt = _extract_snp_from_hgvs(g_str)

        # --- Priority 2: HGVS INDEL (if SNP extraction failed) ---
        if pos is None and include_indels:
            pos, ref, alt = _extract_indel_from_hgvs(g_str)

        # --- Priority 3: Try c./n. string + generic SNP pattern ---
        if pos is None:
            search_text = " | ".join([g_str, c_str, allele_str])
            gm = _GENERIC_SNP.search(search_text)
            if gm:
                ref_cand, alt_cand = gm.group(1).upper(), gm.group(2).upper()
                # Get position from 'position' col or g_str number
                if pos_col:
                    pm = _POS_NUM.search(str(row.get(pos_col, "")))
                    if pm:
                        pos = int(pm.group(1))
                if pos is None:
                    pm = _POS_NUM.search(g_str)
                    if pm:
                        pos = int(pm.group(1))
                if pos:
                    ref, alt = ref_cand, alt_cand

        # --- Skip if we still couldn't parse a variant ---
        if pos is None or pos <= 0:
            continue
        if ref == "-" and alt == "-":
            continue
        if not ref or not alt:
            continue
        # Skip clearly bad alleles (numeric garbage, etc.)
        if not _is_valid_allele(ref) or not _is_valid_allele(alt):
            continue
        if chrom.lower() in ("", "nan", "none", "?"):
            continue

        # --- Clean chromosome: strip 'chr' prefix for Ensembl format ---
        clean_chrom = re.sub(r"^(?:chr|Chr|CHR)", "", chrom).strip()
        if not clean_chrom or clean_chrom.lower() in ("nan", "none"):
            continue

        # --- Build rationale ---
        gene = str(row.get(gene_col, "Unknown")).strip() if gene_col else "Unknown"
        phene = str(row.get(phene_col, "Unknown")).strip() if phene_col else "Unknown"
        effect = (
            str(row.get(effect_col, "Unknown")).strip() if effect_col else "Unknown"
        )
        desc = str(row.get(desc_col, "")).strip() if desc_col else ""

        if desc and desc.lower() not in ("nan", "none", "") and len(desc) > 15:
            rationale = desc
        else:
            rationale = f"[{gene}] {phene} — {effect}"

        # --- Label: OMIA Mendelian single-locus variants are curated disease/causal ---
        is_path = True

        # --- OMIA ID ---
        omia_id = str(row.get(omia_id_col, np.nan)) if omia_id_col else np.nan

        # --- Variant type ---
        if ref == "-" or alt == "-" or len(ref) != len(alt):
            v_type = "indel"
        else:
            v_type = "snp"

        records.append(
            {
                "chrom": clean_chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "label": is_path,
                "label_source": "omia_curated_mendelian",
                "rationale": rationale,
                "species": cfg["token"],
                "gene": gene,
                "omia_id": omia_id,
                "variant_phenotype": phene,
                "variant_effect": effect,
                "variant_type": v_type,
            }
        )

    result_df = pd.DataFrame(records)
    if result_df.empty:
        print(f"   ⚠️  No parseable variants found for {cfg['common_name']}")
        return result_df

    # Deduplicate on (chrom, pos, ref, alt)
    n_before = len(result_df)
    result_df = result_df.drop_duplicates(
        subset=["chrom", "pos", "ref", "alt"]
    ).reset_index(drop=True)
    if n_before > len(result_df):
        print(f"   🔄 Deduplicated: {n_before} → {len(result_df)} variants")

    print(
        f"   ✅ {cfg['common_name']}: {len(result_df)} variants parsed "
        f"({result_df['variant_type'].value_counts().to_dict()})"
    )
    return result_df


# ============================================================================
# COORDINATE VALIDATION
# ============================================================================


def validate_reference_bases(
    df: pd.DataFrame,
    fasta_path: Path,
    max_checks: int = 50,
) -> float:
    """
    Spot-check ref alleles against reference genome using samtools faidx.

    Returns the fraction of matches (1.0 = perfect).
    Skips if faidx index is missing.
    """
    fai_path = Path(str(fasta_path) + ".fai")
    if not fasta_path.exists() or not fai_path.exists():
        print(
            f"   ⚠️  Skipping coordinate validation (genome not indexed): {fasta_path}"
        )
        return float("nan")

    # Only validate SNPs with real nucleotides (skip N placeholders from indels)
    snp_df = df[
        (df["ref"].str.len() == 1)
        & (df["ref"].str.upper().isin(["A", "C", "G", "T"]))
        & (df["alt"] != "-")
    ].head(max_checks)
    if snp_df.empty:
        print("   ℹ️  No SNPs to validate (all indels?)")
        return float("nan")

    checks = len(snp_df)
    matches = 0
    mismatches = []

    for _, row in snp_df.iterrows():
        region = f"{row['chrom']}:{row['pos']}-{row['pos']}"
        result = subprocess.run(
            ["samtools", "faidx", str(fasta_path), region],
            capture_output=True,
            text=True,
        )
        lines = result.stdout.strip().split("\n")
        if len(lines) >= 2:
            fasta_ref = lines[1].strip().upper()
            if fasta_ref == row["ref"].upper():
                matches += 1
            else:
                mismatches.append(
                    f"{row['chrom']}:{row['pos']} expected {row['ref']} got {fasta_ref}"
                )

    frac = matches / checks if checks > 0 else 0.0
    if frac >= 0.9:
        print(f"   ✅ Coordinate validation: {matches}/{checks} matched ({frac:.0%})")
    else:
        print(f"   ⚠️  Coordinate validation: {matches}/{checks} matched ({frac:.0%})")
        print("       Possible assembly mismatch! Mismatches (first 5):")
        for m in mismatches[:5]:
            print(f"         {m}")
    return frac


# ============================================================================
# EDA / SUMMARY
# ============================================================================


def eda_summary(df: pd.DataFrame, species_key: str) -> None:
    """
    Print exploratory data analysis summary for a wrangled species DataFrame.
    """
    cfg = SPECIES_CONFIG[species_key]
    print(f"\n{'=' * 70}")
    print(f"  EDA: {cfg['common_name']} ({cfg['token']})")
    print(f"{'=' * 70}")
    print(f"  Total variants : {len(df)}")
    if "label" in df.columns:
        print(f"  Labels         : {df['label'].value_counts(dropna=False).to_dict()}")
    if "variant_type" in df.columns:
        print(f"  Variant types  : {df['variant_type'].value_counts().to_dict()}")

    # Chromosomes
    chrom_counts = df["chrom"].value_counts()
    print(f"  Chromosomes    : {len(chrom_counts)} unique")
    print(f"  Top chroms     : {chrom_counts.head(5).to_dict()}")

    # Genes
    if "gene" in df.columns:
        gene_counts = df["gene"].value_counts()
        print(f"  Top genes      : {gene_counts.head(8).to_dict()}")

    # Phenotypes
    if "variant_phenotype" in df.columns:
        phene_counts = df["variant_phenotype"].value_counts()
        print(f"  Top phenotypes : {phene_counts.head(5).to_dict()}")

    # Rationale quality
    if "rationale" in df.columns:
        long_rat = (df["rationale"].str.len() > 30).sum()
        print(f"  Rich rationale : {long_rat}/{len(df)} (>30 chars)")
        print("\n  Sample rationales:")
        for _, row in df.dropna(subset=["rationale"]).head(3).iterrows():
            rat = row["rationale"][:120]
            print(f"    [{row.get('gene', '?')}] {rat}")

    # Position ranges
    print(f"\n  Pos range      : {df['pos'].min():,} – {df['pos'].max():,}")
    print(f"  Ref alleles    : {df['ref'].value_counts().head(6).to_dict()}")
    print(f"  Alt alleles    : {df['alt'].value_counts().head(6).to_dict()}")
    print(f"{'=' * 70}")


def validate_inference_output(rdf: pd.DataFrame, species_key: str) -> None:
    """
    Validate NTv3 inference output for an animal species.

    Reuses base validate_output and adds animal-specific checks.
    """
    print(f"\n  Inference validation — {SPECIES_CONFIG[species_key]['common_name']}")
    base.validate_output(rdf)

    # Check: no BigWig columns should be present for non-human
    bw_cols = [c for c in rdf.columns if c.startswith(("REF_BW_", "D_BW_"))]
    if bw_cols:
        print(
            f"  ⚠️  Unexpected BigWig columns present ({len(bw_cols)}): "
            f"{bw_cols[:3]} ..."
        )
    else:
        print("  ✅ No BigWig columns (expected for non-human species)")

    # Check: species column preserved
    if "species" in rdf.columns:
        sp_vals = rdf["species"].unique()
        print(f"  ✅ Species column present: {sp_vals}")
    else:
        print("  ⚠️  Species column missing from output")

    # Check: BED columns present (21 elements)
    bed_cols = [c for c in rdf.columns if c.startswith("REF_BED_")]
    print(f"  BED REF columns: {len(bed_cols)}")
    if bed_cols:
        vals = rdf[bed_cols].values.flatten()
        vals = vals[~np.isnan(vals)]
        if len(vals):
            in_bounds = ((vals >= 0) & (vals <= 1)).mean()
            print(f"  BED prob bounds: {in_bounds:.0%} in [0,1]")


# ============================================================================
# INFERENCE RUNNER
# ============================================================================


def run_species_inference(
    species_key: str,
    df: pd.DataFrame,
    model,
    tokenizer,
    device: str,
    sample_size: Optional[int] = SAMPLE_SIZE,
) -> Optional[pd.DataFrame]:
    """
    Run NTv3 inference for one animal species.

    Loads the species-specific genome FASTA, builds feature name lists and
    calls base.run_inference() with the correct species token.

    BED names come from the shared model config (21 elements, all species).
    BigWig is disabled for non-human species: the NTv3 MultiSpeciesHead
    uses a ZeroHead per non-human species, producing logit=0 → sigmoid=0.5
    constant for both ref and alt, so all deltas would be 0.0.

    Returns:
      DataFrame of inference results, or None on error.
    """
    cfg = SPECIES_CONFIG[species_key]
    token = cfg["token"]

    if not token:
        print(f"   ⏭️  Skipping {cfg['common_name']} — not in NTv3 model config")
        return None

    if df.empty:
        print(f"   ⏭️  Skipping {cfg['common_name']} — no variants")
        return None

    # Sample if requested
    run_df = df.head(sample_size).copy() if sample_size else df.copy()
    print(
        f"\n🧬 Running inference: {cfg['common_name']} "
        f"({len(run_df)} variants, species='{token}')"
    )

    # Load genome FASTA
    fa_path = _genome_path(species_key)
    if fa_path.exists():
        from pyfaidx import Fasta

        try:
            genome = Fasta(str(fa_path))
            print(
                f"   📂 Genome loaded: {fa_path.name} ({len(genome.keys())} sequences)"
            )
        except Exception as e:
            warnings.warn(f"Failed to load genome {fa_path}: {e}")
            genome = None
    else:
        print(f"   ❌ Genome not found at {fa_path} — refusing to run inference")
        return None

    # BED names: shared across all species (21 elements in model config).
    # BigWig names: only defined for 'human' in model config; the head itself
    # runs for any species, so we pull human names and apply the same subset
    # filter (ChIP-seq key marks + ATAC/DNase + CAGE/FANTOM5/GTEx) as a
    # cross-species regulatory proxy.
    _bw_lookup_species = "human" if USE_BIGWIGS else token
    bed_names, bw_names = base.extract_feature_names(
        model,
        target_species=_bw_lookup_species,
        use_bed=USE_BED,
        use_bigwigs=USE_BIGWIGS,
    )
    if USE_BIGWIGS and bw_names:
        sel_idx, sel_names = base.get_track_indices(bw_names, METADATA_FILE)
        if not sel_idx:  # metadata missing / no matches — fall back to all
            sel_idx = list(range(len(bw_names)))
        else:
            print(
                f"   📊 BigWig subset: {len(sel_idx)}/{len(bw_names)} tracks selected"
            )
    else:
        sel_idx = []
        bw_names = []

    # Extra columns to preserve
    base_cols = {"chrom", "pos", "ref", "alt", "label"}
    extra_cols = [c for c in run_df.columns if c not in base_cols]

    try:
        rdf = base.run_inference(
            run_df,
            model,
            tokenizer,
            genome,
            bed_names,
            bw_names,
            sel_idx,
            device,
            extra_cols=extra_cols,
            use_kl=USE_KL_DIVERGENCE,
            use_embeddings=USE_EMBEDDINGS,
            context_len=CONTEXT_LEN,
            batch_size=BATCH_SIZE,
            use_bed=USE_BED,
            use_bigwigs=USE_BIGWIGS,
            species=token,
        )
    except Exception as e:
        print(f"   ❌ Inference failed for {cfg['common_name']}: {e}")
        import traceback

        traceback.print_exc()
        return None

    return rdf


# ============================================================================
# MAIN PIPELINE
# ============================================================================


def main(
    species_keys=None,
    fetch_only: bool = False,
    validate_only: bool = False,
    full_run: bool = False,
    sample_size: Optional[int] = None,
    model_size: str = MODEL_SIZE,
    force_refresh: bool = False,
    include_indels: bool = True,
):
    """
    Full OMIA pipeline: fetch → parse → validate → infer (sample or full).

    Parameters
    ----------
    species_keys : list[str] or None
        Species to process. None = all supported species.
    fetch_only : bool
        If True, stop after downloading data and genomes (no inference).
    validate_only : bool
        If True, run data parse + validation and stop before model inference.
    full_run : bool
        If True, run inference on all variants (overrides sample_size).
    sample_size : int or None
        Variants per species for inference sample. None = use SAMPLE_SIZE global.
    model_size : str
        '100M' or '650M'.
    force_refresh : bool
        Re-download OMIA data even if cache exists.
    """
    ensure_dirs()

    # Determine which species to run
    supported = [k for k, v in SPECIES_CONFIG.items() if not v.get("unsupported")]
    requested = species_keys if species_keys else DEFAULT_RUN_SPECIES
    # Warn about unsupported
    for k in species_keys or []:
        if k in SPECIES_CONFIG and SPECIES_CONFIG[k].get("unsupported"):
            print(
                f"⚠️  {k} ({SPECIES_CONFIG[k]['common_name']}) is not in NTv3 "
                f"model config — skipping."
            )
    run_keys = [
        k
        for k in requested
        if k in SPECIES_CONFIG and not SPECIES_CONFIG[k].get("unsupported")
    ]

    print("=" * 70)
    print("  OMIA Multi-Species NTv3 Inference Pipeline")
    print("=" * 70)
    print(f"  Species : {run_keys}")
    print(f"  Model   : InstaDeepAI/NTv3_{model_size}_post")
    print(
        f"  Full run: {full_run}  |  Sample size: "
        f"{'all' if full_run else (sample_size or SAMPLE_SIZE)}"
    )
    print(f"  Validate only: {validate_only}")

    # ── 1. Fetch OMIA raw data ────────────────────────────────────────────────
    print("\n── Fetching OMIA data ──")
    raw_df = fetch_raw_omia(force_refresh=force_refresh)
    print(f"   OMIA raw shape: {raw_df.shape}")
    print(f"   Columns: {list(raw_df.columns[:12])} ...")

    # ── 2. Download genomes (if needed) ──────────────────────────────────────
    print("\n── Downloading / checking reference genomes ──")
    genome_paths = {}
    for sk in run_keys:
        print(f"  [{SPECIES_CONFIG[sk]['common_name']}]")
        gp = download_genome(sk)
        genome_paths[sk] = gp

    if fetch_only:
        print("\n✅ --fetch-only mode complete. No inference run.")
        return

    # ── 3. Parse + wrangle per species ───────────────────────────────────────
    print("\n── Parsing OMIA variants per species ──")
    species_dfs = {}
    validation_scores = {}
    for sk in run_keys:
        print(f"\n  [{SPECIES_CONFIG[sk]['common_name']}]")
        df_sp = parse_omia_variants(raw_df, sk, include_indels=include_indels)
        if df_sp.empty:
            continue

        # Save wrangled parquet
        out_path = DATA_DIR / f"{sk}_variants.parquet"
        df_sp.to_parquet(out_path, index=False)
        print(f"   💾 Saved wrangled variants: {out_path} ({len(df_sp)} rows)")

        # EDA
        eda_summary(df_sp, sk)

        # Coordinate validation
        gp = genome_paths.get(sk)
        match_frac = float("nan")
        if gp:
            match_frac = validate_reference_bases(df_sp, gp, max_checks=50)
        validation_scores[sk] = match_frac
        if not np.isnan(match_frac) and match_frac < 0.8:
            print(
                f"   ❌ Validation failed for {SPECIES_CONFIG[sk]['common_name']}: "
                f"{match_frac:.0%} < 80% reference match"
            )
            continue

        species_dfs[sk] = df_sp

    if not species_dfs:
        print("\n❌ No species had parseable variants. Check OMIA data columns.")
        return

    if validate_only:
        print("\n✅ Validation-only run complete. No inference executed.")
        return {
            sk: {
                "n_variants": len(df_sp),
                "label_dist": df_sp["label"].value_counts(dropna=False).to_dict(),
                "validation_match": validation_scores.get(sk, np.nan),
            }
            for sk, df_sp in species_dfs.items()
        }

    # ── 4. Load model once (shared across all species) ───────────────────────
    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        print(f"\n  GPU: {torch.cuda.get_device_name(0)}")
    else:
        print("\n  Running on CPU — will be slow for large datasets")

    model_name = f"InstaDeepAI/NTv3_{model_size}_post"
    print(f"\n── Loading model: {model_name} ──")
    model, tokenizer = base.load_model_and_tokenizer(model_name, device)

    # ── 5. Run inference per species ─────────────────────────────────────────
    print("\n── Running NTv3 inference ──")
    n_inference = None if full_run else (sample_size or SAMPLE_SIZE)
    results_summary = {}

    for sk, df_sp in species_dfs.items():
        rdf = run_species_inference(
            sk, df_sp, model, tokenizer, device, sample_size=n_inference
        )
        if rdf is None:
            continue

        # Validate output
        validate_inference_output(rdf, sk)

        # Save results
        suffix = "" if full_run else f"_sample{n_inference}"
        out_path = RESULTS_DIR / f"{sk}_deltas{suffix}.parquet"
        rdf.to_parquet(out_path, index=False)
        print(f"   💾 Saved: {out_path}  shape={rdf.shape}")

        results_summary[sk] = {
            "n_input": len(df_sp),
            "n_inferred": len(rdf),
            "n_cols": rdf.shape[1],
            "label_dist": rdf["label"].value_counts(dropna=False).to_dict(),
            "bed_cols": len([c for c in rdf.columns if c.startswith("REF_BED_")]),
        }

    # ── 6. Final summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  OMIA Pipeline Complete — Summary")
    print("=" * 70)
    for sk, stats in results_summary.items():
        cfg = SPECIES_CONFIG[sk]
        print(f"\n  {cfg['common_name']} ({cfg['token']}):")
        print(f"    Input variants : {stats['n_input']}")
        print(f"    Inferred       : {stats['n_inferred']}")
        print(
            f"    Output columns : {stats['n_cols']}"
            f"  (BED: {stats['bed_cols']}, BigWig: 0)"
        )
        print(f"    Labels         : {stats['label_dist']}")

    print("\n  Output files:")
    for sk in results_summary:
        suffix = "" if full_run else f"_sample{n_inference}"
        print(f"    Variants : {DATA_DIR}/{sk}_variants.parquet")
        print(f"    Results  : {RESULTS_DIR}/{sk}_deltas{suffix}.parquet")

    print("=" * 70)
    return results_summary


# ============================================================================
# CLI
# ============================================================================


def _parse_args():
    all_species = [k for k, v in SPECIES_CONFIG.items() if not v.get("unsupported")]
    parser = argparse.ArgumentParser(
        description="OMIA multi-species NTv3 variant inference pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--species",
        nargs="+",
        choices=list(SPECIES_CONFIG.keys()),
        default=None,
        metavar="SPECIES",
        help=(
            f"Species to process. Choices: {all_species}. "
            f"Default: {DEFAULT_RUN_SPECIES}."
        ),
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Run full inference on all variants (default: sample only)",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=None,
        metavar="N",
        help=f"Variants per species for sample run (default: {SAMPLE_SIZE})",
    )
    parser.add_argument(
        "--model-size",
        choices=["100M", "650M"],
        default=MODEL_SIZE,
        help=f"NTv3 model size (default: {MODEL_SIZE})",
    )
    parser.add_argument(
        "--fetch-only",
        action="store_true",
        help="Download OMIA data and genomes only — no inference",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Fetch + parse + coordinate-validate only; skip model inference",
    )
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="Re-download OMIA data even if local cache exists",
    )
    parser.add_argument(
        "--no-indels",
        action="store_true",
        help="Skip INDEL variants (keep SNPs only)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    main(
        species_keys=args.species,
        fetch_only=args.fetch_only,
        validate_only=args.validate_only,
        full_run=args.full,
        sample_size=args.sample_size,
        model_size=args.model_size,
        force_refresh=args.force_refresh,
        include_indels=not args.no_indels,
    )
