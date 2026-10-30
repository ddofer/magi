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
        # validation-only: fetch + parse + REF-letter check on 50 variants per
        # species (only where the genome is in data/genomes); writes nothing

    Any run that writes goes through SLURM and should pass --out-dir <folder>:
    the default output folder is the repo's data/, which holds tracked files.

    python inference_omia.py
        # sample inference: default species, 30 variants/species
    python inference_omia.py --species dog cat --sample-size 50
        # sample inference: selected species with custom sample size

    srun -n 1 -p low -w sblxa123c --gres=gpu:1 -c 4 --mem=32G -t 4:00:00 \
        uv run --no-sync python inference_omia.py --full --out-dir scratch/omia
        # full inference (all variants for run species), through SLURM

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
        Note: BigWig columns are off by default (USE_BIGWIGS). Each species has its
            own BigWig head; among these species only mouse has one (2450 tracks).

Gotchas
  1. NTv3 BED tracks (21 genomic elements) work for all species (shared head).
     BigWig tracks come from a per-species head. Dog, cat, chicken, zebrafish and
     rat have none (an empty head), so they never get BigWig columns.
  2. Ensembl chroms use NO 'chr' prefix (1, 2, X, A1, etc.).
     OMIA data sometimes includes 'chr' — we strip it.
  3. Mouse NTv3 token is 'mouse' (not 'mus_musculus').
  4. Chromosome validation uses pyfaidx for spot-checks.
  5. OMIA 'Allele' column is often blank — extract from HGVS g. string.
  6. Only SNPs are supported for LLR feature; INDELs get KL/logprob features.
  7. INDELs are written VCF-style (anchor letter from the species genome), the format
      inference.py expects. HGVS events with unknown letters (e.g. insN[226]), letters that
      disagree with the genome, or spans longer than half the window are skipped and counted.

Authors: Dan Ofer / Michal Linial Lab (based on inference.py)
"""

import argparse
import gzip
import io
import os
import re
import shutil
import urllib.request
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
from pyfaidx import Fasta

import inference as base

# ============================================================================
# CONFIGURATION
# ============================================================================

# Inputs (OMIA cache, genomes) live in the repo's data/, wherever the script runs from.
# Outputs (wrangled variants, results) go to OUT_DIR; --out-dir moves them, so a rerun
# does not overwrite the tracked data/*_variants.parquet and data/results/.
DATA_DIR = Path(__file__).resolve().parent / "data"
GENOME_DIR = DATA_DIR / "genomes"
OMIA_CACHE_FILE = DATA_DIR / "omia_raw_variants.csv"
OUT_DIR = DATA_DIR
RESULTS_DIR = OUT_DIR / "results"
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
USE_BIGWIGS = False  # per-species BigWig heads: of these species only mouse has one (2450 tracks)
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
# HGVS INDEL: g.<start>[_<end>]<del|dup|delins|ins>[ACGT letters]. The lookahead rejects
# anything after the letters (N[226], legacy delAinsGG) rather than matching a shorter event.
_HGVS_G_INDEL = re.compile(r"g\.(\d+)(?:_(\d+))?(delins|del|dup|ins)([ACGT]*)(?![A-Z\[\d])", re.IGNORECASE)
# Any HGVS genomic indel, parseable or not: such rows are never read as SNVs.
_HGVS_G_ANY_INDEL = re.compile(r"g\.\d+(?:_\d+)?(?:del|dup|ins)", re.IGNORECASE)
# Generic SNP fallback: G>A or G/A
_GENERIC_SNP = re.compile(r"\b([ACGT])>([ACGT])\b", re.IGNORECASE)
# Numeric position standalone (last resort)
_POS_NUM = re.compile(r"\b(\d{4,})\b")


# ============================================================================
# DIRECTORY SETUP
# ============================================================================


def ensure_dirs():
    """Create required output directories."""
    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)


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

    Skips files that already exist, so re-runs are safe.
    pyfaidx builds the .fai index.

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
        part = gz_path.with_name(gz_path.name + ".part")  # renamed when complete: a killed run leaves no short file
        try:
            urllib.request.urlretrieve(cfg["genome_url"], part)
            os.replace(part, gz_path)
        except Exception as e:
            part.unlink(missing_ok=True)
            print(
                f"   ❌ Download failed for {species_key}: {e}. "
                f"Try manually: wget -c {cfg['genome_url']}"
            )
            return None

    # Decompress (keeps the .gz)
    if not fa_path.exists() and gz_path.exists():
        print(f"   📦 Decompressing {gz_path.name}...")
        part = fa_path.with_name(fa_path.name + ".part")
        try:
            with gzip.open(gz_path, "rb") as f_in, open(part, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.replace(part, fa_path)
        except Exception as e:
            part.unlink(missing_ok=True)
            if isinstance(e, (EOFError, gzip.BadGzipFile)):
                gz_path.unlink(missing_ok=True)  # a corrupt download would otherwise block every later run
            print(f"   ❌ Decompression failed for {gz_path}: {e}")
            return None

    if not fa_path.exists():
        print(f"   ❌ FASTA not found after download: {fa_path}")
        return None

    # Index with pyfaidx
    if not fai_path.exists():
        print("   🔍 Indexing with pyfaidx...")
        try:
            Fasta(str(fa_path)).close()
        except Exception as e:
            print(f"   ⚠️  Indexing failed ({e}) — coordinate validation will be skipped")

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


def _extract_indel_from_hgvs(g_str: str, genome, chrom: str):
    """
    VCF-style (pos, ref, alt) for an HGVS genomic INDEL, with letters from `genome`.

    Deletions and duplications are anchored on the preceding / last letter, as in VCF.
    Returns (None, None, None) when unparseable, the sequence is unknown or disagrees with the
    genome, the span leaves the contig or is longer than half the window, or letters are not ACGT.
    """
    none = None, None, None
    m = _HGVS_G_INDEL.search(g_str)
    if not m or genome is None:
        return none
    start, end = int(m.group(1)), int(m.group(2) or m.group(1))
    kind, seq = m.group(3).lower(), m.group(4).upper()
    span = end - start + 1
    if start < 2 or span < 1 or max(span, len(seq)) > CONTEXT_LEN // 2:
        return none
    fetch = lambda a, b: base._fetch(genome, chrom, a - 1, b) or ""  # 1-based inclusive
    region = fetch(start, end)
    if len(region) != span:
        return none
    if kind == "delins":
        pos, ref, alt = start, region, seq
    elif kind == "del":
        pos, ref = start - 1, fetch(start - 1, start - 1) + region
        alt = ref[:1]
    elif kind == "dup":
        pos, ref, alt = end, region[-1], region[-1] + region
    else:  # ins between start and start+1
        if m.group(2) and end != start + 1:
            return none
        pos, ref = start, region[:1]
        alt = ref + seq
    if kind in ("del", "dup") and seq and seq != region:
        return none
    if not alt or len(ref) < 1 or not set(ref + alt) <= set("ACGT"):
        return none
    return pos, ref, alt


def _is_valid_allele(allele: str) -> bool:
    """Accept only nucleotide strings."""
    if not allele or allele in ("nan", "None", "."):
        return False
    return bool(re.match(r"^[ACGT]+$", allele.upper()))


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
    genome=None,
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
    n_indel_skipped = 0
    if include_indels and genome is None:
        print("   ⚠️  No genome: HGVS indels are skipped (their letters come from the genome)")

    for _, row in species_df.iterrows():
        chrom = str(row[chrom_col]).strip() if chrom_col else ""
        clean_chrom = re.sub(r"^(?:chr|Chr|CHR)", "", chrom).strip()
        g_str = str(row[gm_col]).strip() if gm_col else ""
        c_str = str(row[cn_col]).strip() if cn_col else ""
        allele_str = str(row[allele_col]).strip() if allele_col else ""

        ref, alt, pos = None, None, None

        # --- Priority 1: HGVS genomic substitution ---
        pos, ref, alt = _extract_snp_from_hgvs(g_str)

        # --- Priority 2: HGVS INDEL (if SNP extraction failed) ---
        if pos is None and _HGVS_G_ANY_INDEL.search(g_str):
            if include_indels:
                pos, ref, alt = _extract_indel_from_hgvs(g_str, genome, clean_chrom)
            if pos is None:  # never fall back to reading an indel row as an SNV
                n_indel_skipped += include_indels
                continue

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
        if not ref or not alt:
            continue
        # Skip clearly bad alleles (numeric garbage, etc.)
        if not _is_valid_allele(ref) or not _is_valid_allele(alt):
            continue
        if chrom.lower() in ("", "nan", "none", "?"):
            continue

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
        v_type = "snp" if len(ref) == len(alt) == 1 else "indel"

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

    if n_indel_skipped:
        print(f"   ⏭️  {n_indel_skipped} HGVS indels skipped (unknown or mismatched letters, too long, off-contig)")
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
    Spot-check ref alleles against reference genome using pyfaidx.

    Returns the fraction of matches (1.0 = perfect).
    Skips if faidx index is missing.
    """
    fai_path = Path(str(fasta_path) + ".fai")
    if not fasta_path.exists() or not fai_path.exists():
        print(
            f"   ⚠️  Skipping coordinate validation (genome not indexed): {fasta_path}"
        )
        return float("nan")

    # Only validate SNPs: INDEL ref letters were read from this genome, so they always match
    snp_df = df[df["variant_type"] == "snp"].head(max_checks)
    if snp_df.empty:
        print("   ℹ️  No SNPs to validate (all indels?)")
        return float("nan")

    checks = len(snp_df)
    matches = 0
    mismatches = []

    with Fasta(str(fasta_path)) as genome:
        for _, row in snp_df.iterrows():
            pos = int(row["pos"])
            fasta_ref = base._fetch(genome, row["chrom"], pos - 1, pos)
            if not fasta_ref:
                continue
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
    BigWig names come from the species' own head (USE_BIGWIGS; only mouse has one).

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
        try:
            genome = Fasta(str(fa_path))
            print(
                f"   📂 Genome loaded: {fa_path.name} ({len(genome.keys())} sequences)"
            )
        except Exception as e:
            # No fallback: without a genome, inference.py would fetch hg38 at these coordinates.
            print(f"   ❌ Failed to load genome {fa_path}: {e} — refusing to run inference")
            return None
    else:
        print(f"   ❌ Genome not found at {fa_path} — refusing to run inference")
        return None

    # BED names: shared across all species (21 elements in model config).
    # BigWig names: the species' own head. Of the OMIA species only mouse has one
    # (2450 tracks); the others have none, so they get no BigWig columns.
    bed_names, bw_names = base.extract_feature_names(
        model,
        target_species=token,
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
    out_dir: Optional[str] = None,
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
        Re-download OMIA data even if cache exists (overwrites data/omia_raw_variants.csv).
    out_dir : str or None
        Folder for wrangled variants and results. None = the repo's data/.
    """
    global OUT_DIR, RESULTS_DIR
    if out_dir:
        OUT_DIR = Path(out_dir).resolve()
        RESULTS_DIR = OUT_DIR / "results"
    if validate_only:
        # Read-only: existing OMIA cache and genomes only, no downloads, no writes.
        if force_refresh or not OMIA_CACHE_FILE.exists():
            raise SystemExit("--validate-only needs data/omia_raw_variants.csv and no --force-refresh")
    else:
        ensure_dirs()

    # Determine which species to run
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
        if validate_only:
            gp = _genome_path(sk) if _genome_path(sk).exists() else None
        else:
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
        gp = genome_paths.get(sk)
        genome = Fasta(str(gp)) if gp and gp.exists() else None
        df_sp = parse_omia_variants(raw_df, sk, include_indels=include_indels, genome=genome)
        if genome is not None:
            genome.close()
        if df_sp.empty:
            continue

        if not validate_only:
            out_path = OUT_DIR / f"{sk}_variants.parquet"
            df_sp.to_parquet(out_path, index=False)
            print(f"   💾 Saved wrangled variants: {out_path} ({len(df_sp)} rows)")

        # EDA
        eda_summary(df_sp, sk)

        # Coordinate validation
        gp = genome_paths.get(sk)
        match_frac = float("nan")
        if gp:
            match_frac = validate_reference_bases(df_sp, gp, max_checks=50)
        else:
            print("   ⚠️  No genome: REF letters NOT checked for this species.")
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
        print(f"    Variants : {OUT_DIR}/{sk}_variants.parquet")
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
        help="Fetch + parse + check REF letters (50 per species, needs the genome); "
        "writes nothing, skips model inference",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Folder for wrangled variants and results (default: the repo's data/, "
        "which holds tracked files)",
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
        out_dir=args.out_dir,
    )
