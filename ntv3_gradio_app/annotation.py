#!/usr/bin/env python3
"""
MANE Genomic Annotation Module
=======================================================
Adapted from stav_analysis_v2.py for the Gradio app.

Provides gene-level and structural annotation (exon, intron, UTR, splice, promoter)
using the MANE Select transcript dataset (RefSeq).

Usage:
    from annotation import annotate_dataframe
    df_annotated = annotate_dataframe(df)  # adds 30+ annotation columns
"""
from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd
import numpy as np

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
MANE_FILE = DATA_DIR / "MANE_processed.csv"
MANE_PARQUET = DATA_DIR / "MANE_processed.parquet"
PROMOTER_FILE = DATA_DIR / "Promoter_processed.csv"
PROMOTER_PARQUET = DATA_DIR / "Promoter_processed.parquet"

# Annotation columns (27 region flags + transcript sets)
ANNOTATION_COLUMNS = [
    'gene', 'mRNA', 'mRNA_promoter', 'mRNA_exon', 'coding_sequence',
    'start_codon', 'stop_codon', 'five_prime_UTR', 'three_prime_UTR',
    'mRNA_intron', 'mRNA_splice', 'lncRNA', 'lncRNA_promoter', 'lncRNA_exon',
    'snRNA', 'snRNA_promoter', 'snRNA_exon', 'antisenseRNA',
    'antisenseRNA_promoter', 'antisenseRNA_exon', 'telomeraseRNA',
    'telomeraseRNA_promoter', 'telomeraseRNA_exon', 'RNaseMRPRNA',
    'RNaseMRPRNA_promoter', 'RNaseMRPRNA_exon', 'snoRNA', 'snoRNA_promoter',
    'snoRNA_exon', 'other'
]

RNA_TYPES = ['lncRNA', 'snRNA', 'antisenseRNA', 'telomeraseRNA',
             'RNaseMRPRNA', 'snoRNA']

# Global cache for MANE data
_MANE_CACHE = {
    "mane_by_chrom": None,
    "promoter_by_chrom": None,
    "mane_parent_idx": None,
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def collapse_region_class(region: str) -> str:
    """
    Collapse a comma-separated region annotation into a high-level class.
    
    Priority order:
      CODING > SPLICE > UTR_5 > UTR_3 > PROMOTER > INTRONIC > GENIC_OTHER > OTHER
    
    Args:
        region: Comma-separated string of annotation flags
    
    Returns:
        High-level region class string
    """
    if not isinstance(region, str) or not region.strip():
        return "OTHER"

    parts = {r.strip() for r in region.split(",")}

    if {"coding_sequence", "start_codon", "stop_codon"} & parts:
        if "mRNA_splice" in parts:
            return "CODING, SPLICE"
        return "CODING"

    if "mRNA_splice" in parts:
        return "SPLICE"

    if "five_prime_UTR" in parts:
        return "UTR_5"

    if "three_prime_UTR" in parts:
        return "UTR_3"

    if "mRNA_promoter" in parts:
        return "PROMOTER"

    if "mRNA_intron" in parts:
        return "INTRONIC"

    if parts & {"gene", "mRNA", "mRNA_exon"}:
        return "GENIC_OTHER"

    return "OTHER"


def preprocess(mane_raw: pd.DataFrame, promoter_raw: pd.DataFrame) -> Tuple[Dict, Dict, Dict]:
    """
    Pre-cast types and build per-chromosome lookup structures.
    
    Returns:
        (mane_by_chrom, promoter_by_chrom, mane_parent_idx)
    """
    mane = mane_raw.copy()
    promoter = promoter_raw.copy()

    # Cast integer columns
    mane['Start'] = mane['Start'].astype(np.int64)
    mane['End'] = mane['End'].astype(np.int64)
    promoter['Promoter_Start'] = promoter['Promoter_Start'].astype(np.int64)
    promoter['Promoter_End'] = promoter['Promoter_End'].astype(np.int64)

    # Normalize chromosome naming
    mane['chrom_key'] = mane['Chromosome']
    promoter['chrom_key'] = promoter['Chromosome'].astype(str)

    # Group by chromosome
    mane_by_chrom = {k: g for k, g in mane.groupby('chrom_key')}
    promoter_by_chrom = {k: g for k, g in promoter.groupby('chrom_key')}

    # Build parent-indexed lookups
    mane_parent_idx = {}
    for chrom, df in mane_by_chrom.items():
        parent_groups = {}
        for parent_val, grp in df.groupby('Parent', sort=False):
            parent_groups[parent_val] = grp
        mane_parent_idx[chrom] = parent_groups

    return mane_by_chrom, promoter_by_chrom, mane_parent_idx


def annotate_variant(
    chrom: str,
    pos: int,
    mane_by_chrom: Dict,
    promoter_by_chrom: Dict,
    mane_parent_idx: Dict
) -> Tuple[dict, Set[str], Set[str]]:
    """
    Annotate a single variant at the given position.
    
    Args:
        chrom: Chromosome (e.g., 'chr17')
        pos: 1-based genomic position
        mane_by_chrom: MANE data grouped by chromosome
        promoter_by_chrom: Promoter data grouped by chromosome
        mane_parent_idx: Parent-indexed MANE data
    
    Returns:
        (annotation_dict, transcript_set, promoter_transcript_set)
    """
    result = {col: 0 for col in ANNOTATION_COLUMNS}
    transcript_set = set()
    promoter_transcript_set = set()

    chrom_str = str(chrom)

    # --- Overlap with MANE ---
    mane_df = mane_by_chrom.get(chrom_str)
    if mane_df is not None and len(mane_df) > 0:
        mask = (mane_df['Start'].values <= pos) & (mane_df['End'].values >= pos)
        annotation = mane_df[mask]
    else:
        annotation = pd.DataFrame()

    # --- Overlap with Promoter ---
    prom_df = promoter_by_chrom.get(chrom_str)
    if prom_df is not None and len(prom_df) > 0:
        mask = (prom_df['Promoter_Start'].values <= pos) & (prom_df['Promoter_End'].values >= pos)
        annotation_promoter = prom_df[mask]
    else:
        annotation_promoter = pd.DataFrame()

    if annotation.empty and annotation_promoter.empty:
        result['other'] = 1
        return result, transcript_set, promoter_transcript_set

    types = set(annotation['Feature'].unique()) if not annotation.empty else set()
    types_promoter = set(annotation_promoter['Feature'].unique()) if not annotation_promoter.empty else set()

    # --- gene ---
    if 'gene' in types:
        result['gene'] = 1

    # --- mRNA ---
    if 'mRNA' in types:
        result['mRNA'] = 1
        tids = set(annotation.loc[annotation['Feature'] == 'mRNA', 'transcript_id'].dropna())
        transcript_set.update(tids)

        parent_idx = mane_parent_idx.get(chrom_str, {})

        for tid in tids:
            rna_key = f'rna-{tid}'

            # Strand
            id_match = annotation[annotation['ID'] == rna_key]
            if id_match.empty:
                continue
            strand = id_match['Strand'].iloc[0]

            # Exon/CDS overlapping this position
            ann_exon = annotation[(annotation['Parent'] == rna_key) & (annotation['Feature'] == 'exon')]
            ann_cds = annotation[(annotation['Parent'] == rna_key) & (annotation['Feature'] == 'CDS')]

            # Full transcript exons/CDS
            full_exon = parent_idx.get(rna_key)
            if full_exon is not None:
                tr_exon = full_exon[full_exon['Feature'] == 'exon']
                tr_cds = full_exon[full_exon['Feature'] == 'CDS']
            else:
                tr_exon = pd.DataFrame()
                tr_cds = pd.DataFrame()

            if not ann_cds.empty and not ann_exon.empty:
                # Exon + CDS
                result['mRNA_exon'] = 1
                result['coding_sequence'] = 1

                if not tr_cds.empty:
                    cds_starts = tr_cds['Start'].values
                    cds_ends = tr_cds['End'].values
                    if strand == '+':
                        start_1 = cds_starts.min()
                        if start_1 <= pos <= start_1 + 2:
                            result['start_codon'] = 1
                        stop_3 = cds_ends.max()
                        if stop_3 - 2 <= pos <= stop_3:
                            result['stop_codon'] = 1
                    else:
                        start_1 = cds_ends.max()
                        if start_1 - 2 <= pos <= start_1:
                            result['start_codon'] = 1
                        stop_3 = cds_starts.min()
                        if stop_3 <= pos <= stop_3 + 2:
                            result['stop_codon'] = 1

            elif ann_cds.empty and not ann_exon.empty:
                # UTR
                result['mRNA_exon'] = 1
                if not tr_exon.empty and not tr_cds.empty:
                    exon_starts = tr_exon['Start'].values
                    exon_ends = tr_exon['End'].values
                    cds_starts = tr_cds['Start'].values
                    cds_ends = tr_cds['End'].values

                    if strand == '+':
                        five_start = exon_starts.min()
                        five_end = cds_starts.min() - 1
                        if five_start <= pos <= five_end:
                            result['five_prime_UTR'] = 1
                        three_start = cds_ends.max() + 1
                        three_end = exon_ends.max()
                        if three_start <= pos <= three_end:
                            result['three_prime_UTR'] = 1
                    else:
                        five_start = exon_ends.max()
                        five_end = cds_ends.max() + 1
                        if five_end <= pos <= five_start:
                            result['five_prime_UTR'] = 1
                        three_start = cds_starts.min() - 1
                        three_end = exon_starts.min()
                        if three_end <= pos <= three_start:
                            result['three_prime_UTR'] = 1

            elif ann_cds.empty and ann_exon.empty:
                # Intron
                result['mRNA_intron'] = 1
                if not tr_exon.empty:
                    ex_starts = tr_exon['Start'].values
                    ex_ends = tr_exon['End'].values
                    splice_positions = np.concatenate([
                        ex_starts - 1, ex_starts - 2,
                        ex_ends + 1, ex_ends + 2
                    ])
                    if pos in splice_positions:
                        result['mRNA_splice'] = 1

    # --- mRNA promoter ---
    if 'mRNA' in types_promoter:
        result['mRNA_promoter'] = 1
        tids = set(annotation_promoter.loc[
            annotation_promoter['Feature'] == 'mRNA', 'transcript_id'
        ].dropna())
        promoter_transcript_set.update(tids)

    # --- Other RNA types ---
    for rna in RNA_TYPES:
        if rna in types:
            result[rna] = 1
            tids = set(annotation.loc[annotation['Feature'] == rna, 'transcript_id'].dropna())
            transcript_set.update(tids)
            for tid in tids:
                rna_key = f'rna-{tid}'
                ann_exon = annotation[(annotation['Parent'] == rna_key) & (annotation['Feature'] == 'exon')]
                if not ann_exon.empty:
                    result[f'{rna}_exon'] = 1

        if rna in types_promoter:
            result[f'{rna}_promoter'] = 1
            tids = set(annotation_promoter.loc[
                annotation_promoter['Feature'] == rna, 'transcript_id'
            ].dropna())
            promoter_transcript_set.update(tids)

    return result, transcript_set, promoter_transcript_set


# ============================================================================
# PUBLIC API
# ============================================================================
def _load_or_convert(csv_path: Path, parquet_path: Path) -> pd.DataFrame:
    """Load from parquet if available, otherwise read CSV and cache as parquet."""
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    df = pd.read_csv(csv_path)
    try:
        df.to_parquet(parquet_path, index=False)
        print(f"   Cached {parquet_path.name} for faster future loads")
    except Exception as exc:
        print(f"   ⚠️  Parquet cache write failed: {exc}")
    return df


def load_mane_data():
    """Load and preprocess MANE and Promoter data. Caches globally."""
    if _MANE_CACHE["mane_by_chrom"] is not None:
        return  # Already loaded

    print(f"📚 Loading MANE annotation data from {DATA_DIR}...")
    mane_raw = _load_or_convert(MANE_FILE, MANE_PARQUET)
    promoter_raw = _load_or_convert(PROMOTER_FILE, PROMOTER_PARQUET)

    print(f"   MANE: {len(mane_raw):,} features")
    print(f"   Promoter: {len(promoter_raw):,} features")

    mane_by_chrom, promoter_by_chrom, mane_parent_idx = preprocess(mane_raw, promoter_raw)

    _MANE_CACHE.update({
        "mane_by_chrom": mane_by_chrom,
        "promoter_by_chrom": promoter_by_chrom,
        "mane_parent_idx": mane_parent_idx,
    })

    print(f"✅ MANE data loaded: {len(mane_by_chrom)} chromosomes indexed")


def annotate_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add MANE annotations to a variants DataFrame.
    
    Args:
        df: DataFrame with 'chrom' and 'pos' columns
    
    Returns:
        DataFrame with added annotation columns:
        - 30 binary flags (gene, mRNA, coding_sequence, start_codon, etc.)
        - 'transcript_set', 'promoter_transcript_set' (sets of transcript IDs)
        - 'region' (comma-separated list of active flags)
        - 'region_class' (high-level category)
        - 'gene_name' (from MANE, if available)
    """
    # Load MANE data if not already loaded
    if _MANE_CACHE["mane_by_chrom"] is None:
        load_mane_data()

    mane_by_chrom = _MANE_CACHE["mane_by_chrom"]
    promoter_by_chrom = _MANE_CACHE["promoter_by_chrom"]
    mane_parent_idx = _MANE_CACHE["mane_parent_idx"]

    # Validate input
    if "chrom" not in df.columns or "pos" not in df.columns:
        raise ValueError("DataFrame must have 'chrom' and 'pos' columns")

    df = df.copy()
    chroms = df['chrom'].values
    positions = df['pos'].astype(np.int64).values

    all_results = []
    all_tsets = []
    all_ptsets = []

    for i in range(len(df)):
        res, tset, ptset = annotate_variant(
            chroms[i], positions[i],
            mane_by_chrom, promoter_by_chrom, mane_parent_idx
        )
        all_results.append(res)
        all_tsets.append(tset)
        all_ptsets.append(ptset)

    # Add annotation columns
    ann_df = pd.DataFrame(all_results, index=df.index)
    for col in ANNOTATION_COLUMNS:
        df[col] = ann_df[col].values
    df['transcript_set'] = all_tsets
    df['promoter_transcript_set'] = all_ptsets

    # Combine into region string
    df['region'] = (
        df[ANNOTATION_COLUMNS]
        .apply(lambda r: ','.join(r.index[r == 1]), axis=1)
    )

    # Collapse to high-level class
    df["region_class"] = df["region"].apply(collapse_region_class)

    # Extract gene name from MANE (if available)
    gene_names = []
    for i in range(len(df)):
        chrom_str = str(chroms[i])
        pos = positions[i]
        mane_df = mane_by_chrom.get(chrom_str)
        gene_name = ""
        if mane_df is not None:
            mask = (mane_df['Start'].values <= pos) & (mane_df['End'].values >= pos)
            overlaps = mane_df[mask]
            if not overlaps.empty and 'gene' in overlaps.columns:
                genes = overlaps['gene'].dropna().unique()
                if len(genes) > 0:
                    gene_name = genes[0]
        gene_names.append(gene_name)
    
    df['gene_name'] = gene_names

    return df
