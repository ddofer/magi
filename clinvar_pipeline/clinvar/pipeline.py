"""High-level pipeline stages."""

from __future__ import annotations

import gc
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from clinvar.clinvar_enrichment import load_submission_summary
from clinvar.region_annotation import annotate_clinvar
from clinvar.signal_extraction_threshold import extract_signals_with_thresholds
from clinvar.variant_prep import process_indels, process_snps


def prepare_variants(paths: dict, *, run_snp: bool = True, run_indel: bool = True) -> None:
    """Stage 1: strict filter + MANE region annotation → annotated parquets."""
    Path(paths["snp_strict"]).parent.mkdir(parents=True, exist_ok=True)

    snp_strict = indel_strict = None
    sub_df = None
    funnel: dict[str, dict[str, int]] = {}

    if run_snp:
        sub_df = load_submission_summary(paths["submission_summary"])
        snp_strict, funnel["snp"] = process_snps(paths, sub_df=sub_df)
        del sub_df
        sub_df = None
        gc.collect()
    if run_indel:
        sub_df = load_submission_summary(paths["submission_summary"])
        indel_strict, funnel["indel"] = process_indels(paths, sub_df=sub_df)
        del sub_df
        gc.collect()

    if not run_snp and not run_indel:
        raise ValueError("At least one of run_snp / run_indel must be True")

    MANE = pd.read_csv(paths["MANE_processed"])
    Promoter = pd.read_csv(paths["Promoter_processed"])
    print("\n══ PART 3: Region Annotation ══")

    if run_snp and snp_strict is not None:
        snp_annotated = annotate_clinvar(snp_strict, MANE, Promoter)
        snp_annotated.to_parquet(paths["snp_annotated"], index=False)
        print(f"  Saved → {paths['snp_annotated']}  ({len(snp_annotated):,} variants)")

    if run_indel and indel_strict is not None:
        indel_annotated = annotate_clinvar(indel_strict, MANE, Promoter)
        indel_annotated.to_parquet(paths["indel_annotated"], index=False)
        print(f"  Saved → {paths['indel_annotated']}  ({len(indel_annotated):,} variants)")

    if funnel:
        from clinvar.variant_prep import write_cohort_funnel

        write_cohort_funnel(paths, funnel)

    print("\n✓ Variant prep complete.")


def extract_signals(
    paths: dict,
    variant_type: str,
    k_signals: int = 5,
) -> pd.DataFrame:
    """Stage 2: top-k BED/BW/MLM signal columns → signaled parquet (legacy threshold path)."""
    if variant_type == "indel":
        input_file = paths["indel_annotated"]
        output_file = paths["indel_annotated_signaled"]
        thresholds_file = paths["indel_z_thresholds"]
    elif variant_type == "snp":
        input_file = paths["snp_annotated"]
        output_file = paths["snp_annotated_signaled"]
        thresholds_file = paths["snp_z_thresholds"]
    else:
        raise ValueError(f"Unknown variant type: {variant_type}")

    print("Loading data...")
    variant_df = pq.read_table(input_file).to_pandas()
    print(f"Loaded {len(variant_df)} variants")

    print("Loading track metadata...")
    try:
        metadata_df = pd.read_csv(paths["tracks_metadata"])
        print(f"Loaded metadata for {len(metadata_df)} tracks")
    except OSError as e:
        print(f"Warning: Could not load metadata ({e}). BW signals will use track IDs.")
        metadata_df = None

    bed_n = sum(c.startswith("D_BED_") for c in variant_df.columns)
    bw_n = sum(c.startswith("D_BW_") for c in variant_df.columns)
    print(f"Found {bed_n} BED delta columns, {bw_n} BW delta columns")

    variant_df = extract_signals_with_thresholds(
        variant_df,
        metadata_df,
        thresholds_path=thresholds_file,
        k=k_signals,
        variant_type=variant_type,
    )

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    variant_df.to_parquet(output_file)
    print(f"\nResults saved to: {output_file}")
    return variant_df
