"""MANE / Promoter region annotation for filtered variants."""

from __future__ import annotations

import numpy as np
import pandas as pd

from clinvar.constants import ANNOTATION_COLUMNS, RNA_TYPES

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def collapse_region_class(region: str) -> str:
    """
    Collapse a comma-separated region annotation string into a single
    high-level genomic region class for prompting.

    Priority order (highest → lowest):
      CODING > SPLICE > UTR_5 > UTR_3 > PROMOTER > INTRONIC > GENIC_OTHER > OTHER
    """

    if not isinstance(region, str) or not region.strip():
        return "OTHER"

    parts = {r.strip() for r in region.split(",")}

    # --- Priority-based classification ---
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

# ---------------------------------------------------------------------------
# 1. One-time preprocessing — call once before annotating
# ---------------------------------------------------------------------------
def preprocess(MANE_raw, Promoter_raw):
    """
    Pre-cast types and build per-chromosome lookup structures.
    Returns (mane_by_chrom, promoter_by_chrom, mane_full_by_chrom).
    """
    MANE = MANE_raw.copy()
    Promoter = Promoter_raw.copy()

    # Cast once
    MANE['Start'] = MANE['Start'].astype(np.int64)
    MANE['End'] = MANE['End'].astype(np.int64)
    Promoter['Promoter_Start'] = Promoter['Promoter_Start'].astype(np.int64)
    Promoter['Promoter_End'] = Promoter['Promoter_End'].astype(np.int64)

    # Normalize chromosome naming — strip 'chr' prefix in MANE to match VCF style
    # Adjust this depending on your actual data
    MANE['chrom_key'] = MANE['Chromosome']
    Promoter['chrom_key'] = Promoter['Chromosome'].astype(str)

    # Group by chromosome for O(1) chrom lookup
    mane_by_chrom = {k: g for k, g in MANE.groupby('chrom_key')}
    promoter_by_chrom = {k: g for k, g in Promoter.groupby('chrom_key')}

    # For the "full MANE" lookups (transcript exons/CDS across whole chrom)
    # we also build parent-indexed lookups
    mane_parent_idx = {}
    for chrom, df in mane_by_chrom.items():
        parent_groups = {}
        for parent_val, grp in df.groupby('Parent', sort=False):
            parent_groups[parent_val] = grp
        mane_parent_idx[chrom] = parent_groups

    return mane_by_chrom, promoter_by_chrom, mane_parent_idx


# ---------------------------------------------------------------------------
# 2. Fast overlap query using sorted arrays + searchsorted
# ---------------------------------------------------------------------------
def build_interval_index(df, start_col='Start', end_col='End'):
    """
    Returns (starts, ends, indices) sorted by start for searchsorted.
    """
    starts = df[start_col].values
    ends = df[end_col].values
    order = np.argsort(starts)
    return starts[order], ends[order], df.index.values[order]


def query_overlaps(starts_sorted, ends_sorted, idx_sorted, pos):
    """
    Find all intervals containing pos.
    intervals where start <= pos AND end >= pos.
    """
    # All intervals that start <= pos
    right = np.searchsorted(starts_sorted, pos, side='right')  # first index where start > pos
    # Among those (0..right-1), filter end >= pos
    if right == 0:
        return np.array([], dtype=np.int64)
    candidate_ends = ends_sorted[:right]
    mask = candidate_ends >= pos
    return idx_sorted[:right][mask]


# ---------------------------------------------------------------------------
# 3. Per-variant annotation (operates on dicts, not pd.Series)
# ---------------------------------------------------------------------------
def annotate_variant(chrom, pos, mane_by_chrom, promoter_by_chrom, mane_parent_idx):
    """
    Annotate a single variant. Returns dict of annotation flags + transcript sets.
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

        # Get parent-index for this chrom
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

            # Full transcript exons/CDS (from full MANE, not just overlapping)
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
                    # Splice: pos within 2bp of any exon boundary
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


# ---------------------------------------------------------------------------
# 4. Main vectorized driver — no multiprocessing overhead for small-medium data
# ---------------------------------------------------------------------------
def annotate_clinvar(ClinVar, MANE_raw, Promoter_raw, use_parallel=False, n_workers=None):
    """
    Annotate all ClinVar variants.

    For datasets < 100k variants, single-process with pre-indexed lookups
    is usually faster than multiprocessing (avoids serialization overhead).

    For > 100k, set use_parallel=True.
    """
    from tqdm import tqdm

    # Preprocess once
    mane_by_chrom, promoter_by_chrom, mane_parent_idx = preprocess(MANE_raw, Promoter_raw)

    chroms = ClinVar['chrom'].values
    positions = ClinVar['pos'].astype(np.int64).values

    all_results = []
    all_tsets = []
    all_ptsets = []

    if use_parallel:
        from concurrent.futures import ProcessPoolExecutor
        import functools

        # For true parallelism, partition by chromosome to avoid sharing state
        # But for simplicity, use ThreadPoolExecutor (GIL-bound but avoids pickle)
        from concurrent.futures import ThreadPoolExecutor

        if n_workers is None:
            import os
            n_workers = os.cpu_count()

        def _worker(i):
            return annotate_variant(
                chroms[i], positions[i],
                mane_by_chrom, promoter_by_chrom, mane_parent_idx
            )

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = list(tqdm(
                executor.map(_worker, range(len(chroms))),
                total=len(chroms),
                desc="Annotating variants"
            ))
            for res, tset, ptset in futures:
                all_results.append(res)
                all_tsets.append(tset)
                all_ptsets.append(ptset)
    else:
        for i in tqdm(range(len(chroms)), desc="Annotating variants"):
            res, tset, ptset = annotate_variant(
                chroms[i], positions[i],
                mane_by_chrom, promoter_by_chrom, mane_parent_idx
            )
            all_results.append(res)
            all_tsets.append(tset)
            all_ptsets.append(ptset)

    # Build result DataFrame
    ann_df = pd.DataFrame(all_results, index=ClinVar.index)
    for col in ANNOTATION_COLUMNS:
        ClinVar[col] = ann_df[col].values
    ClinVar['transcript_set'] = all_tsets
    ClinVar['promoter_transcript_set'] = all_ptsets

    ClinVar['region'] = (
    ClinVar[ANNOTATION_COLUMNS]
    .apply(lambda r: ','.join(r.index[r == 1]), axis=1)
    )


    ClinVar["region_class"] = ClinVar["region"].apply(collapse_region_class)

    return ClinVar



