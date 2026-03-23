#!/usr/bin/env python3
"""
NTv3 Genomic Variant Analysis Pipeline

Downstream analysis of NTv3 inference results for semantic validation.

📊 FEATURES:
   - Normalized impact scores (Z-score, MAD, or separate BED/BW)
   - Statistical comparisons (pathogenic vs benign)
   - Per-element significance testing (Bonferroni-corrected)
   - Variant fingerprint visualizations
   - Effect size calculations (Cohen's d, AUC)

🔬 ANALYSIS METHODS:
   1. Impact Score Computation:
      - BED: Top-3 mean absolute delta (probabilities)
      - BW: Top-10 mean absolute delta with normalization
      - Z-score normalization using benign distribution

   2. Statistical Tests:
      - T-tests for pathogenic vs benign differences
      - Cohen's d for effect sizes
      - Bonferroni correction for multiple testing
      - Mann-Whitney U for AUC proxy

   3. Visualizations:
      - Variant fingerprints (bar plots)
      - Combined BED + BigWig views
      - Per-element comparison plots

💻 USAGE:
   # Basic analysis:
   python nt3_analysis.py --input clinvar_new_deltas.parquet

   # Debug mode (first 70 variants):
   python nt3_analysis.py --input clinvar_new_deltas.parquet --debug

   # In Python:
   from nt3_analysis import compute_normalized_impact_scores, plot_variant_fingerprint
   df = pd.read_parquet('results.parquet')
   df = compute_normalized_impact_scores(df, normalization='z_score')
   plot_variant_fingerprint(df, variant_index=0)

📈 EXPECTED OUTPUTS:
   - Normalized impact scores (pathogenic > benign)
   - Cohen's d effect sizes (typically 0.5-0.8 for good discrimination)
   - P-values < 0.05 for significant elements
   - Variant fingerprints showing top disrupted features

Authors: Dan Ofer / Michal Linial Lab
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

pd.set_option("display.max_columns", 40)


# ============================================================================
# NORMALIZED SCORING FUNCTIONS
# ============================================================================


def compute_normalized_impact_scores(df, normalization="z_score"):
    """
    Compute normalized impact scores that properly handle BED and BW scale differences.

    Parameters:
    -----------
    df : DataFrame
        Results dataframe with D_BED_* and D_BW_* columns
    normalization : str
        - 'z_score': Z-score normalize BW tracks using benign distribution
        - 'mad': Use MAD (Median Absolute Deviation) for robust normalization
        - 'separate': Keep BED and BW scores separate (no mixing)

    Returns:
    --------
    df : DataFrame (modified in place)
        Adds columns: Impact_Score_BED, Impact_Score_BW, Impact_Score_Normalized
    """
    delta_bed_cols = [c for c in df.columns if c.startswith("D_BED_")]
    delta_bw_cols = [c for c in df.columns if c.startswith("D_BW_")]

    # BED impact (top-k mean absolute delta)
    if delta_bed_cols:
        bed_matrix = df[delta_bed_cols].values
        k_bed = min(3, bed_matrix.shape[1])
        top_k_bed = np.sort(np.abs(bed_matrix), axis=1)[:, -k_bed:]
        df["Impact_Score_BED"] = np.mean(top_k_bed, axis=1)

    # BigWig impact with normalization
    if delta_bw_cols:
        bw_matrix = df[delta_bw_cols].values
        benign_mask = df["label"] == False  # Benign variants

        if normalization == "z_score":
            # Z-score: normalize by benign std deviation
            benign_matrix = bw_matrix[benign_mask]
            sigma = np.std(benign_matrix, axis=0).copy()
            sigma[sigma < 1e-6] = 1.0
            normalized_bw = np.abs(bw_matrix) / sigma[np.newaxis, :]
            df["Impact_Score_BW"] = np.mean(normalized_bw, axis=1)

        elif normalization == "mad":
            # MAD: robust to outliers
            benign_matrix = bw_matrix[benign_mask]
            median = np.median(np.abs(benign_matrix), axis=0)
            mad = np.median(
                np.abs(benign_matrix - median[np.newaxis, :]), axis=0
            ).copy()
            mad[mad < 1e-6] = 1.0
            normalized_bw = np.abs(bw_matrix) / mad[np.newaxis, :]
            df["Impact_Score_BW"] = np.mean(normalized_bw, axis=1)

        else:  # 'separate' or raw
            k_bw = min(10, bw_matrix.shape[1])
            top_k_bw = np.sort(np.abs(bw_matrix), axis=1)[:, -k_bw:]
            df["Impact_Score_BW"] = np.mean(top_k_bw, axis=1)

    # Combined score
    if delta_bed_cols and delta_bw_cols:
        df["Impact_Score_Normalized"] = (
            df["Impact_Score_BED"] + df["Impact_Score_BW"]
        ) / 2

    return df


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================


def plot_variant_fingerprint(df, variant_index=0, delta_prefix="D_BED_", top_k=15):
    """
    Plot top disrupted features for a single variant.

    Parameters:
    -----------
    df : DataFrame
        Results dataframe
    variant_index : int
        Row index of variant to plot
    delta_prefix : str
        Column prefix for deltas ("D_BED_" or "D_BW_")
    top_k : int
        Number of top features to show
    """
    delta_cols = [c for c in df.columns if c.startswith(delta_prefix)]

    if not delta_cols:
        print(f"⚠️ No columns found with prefix '{delta_prefix}'")
        return

    if variant_index >= len(df):
        print(f"⚠️ Index {variant_index} out of range")
        return

    row = df.iloc[variant_index]

    # Extract deltas
    deltas = {}
    for col in delta_cols:
        name = col.replace(delta_prefix, "")
        val = row[col]
        if isinstance(val, pd.Series):
            val = float(val.iloc[0])
        deltas[name] = float(val)

    # Sort by absolute magnitude
    sorted_deltas = sorted(deltas.items(), key=lambda x: abs(x[1]), reverse=True)[
        :top_k
    ]

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    names = [x[0] for x in sorted_deltas]
    vals = [x[1] for x in sorted_deltas]
    colors = ["#d73027" if v > 0 else "#4575b4" for v in vals]

    ax.barh(names, vals, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Delta (Alt - Ref)", fontsize=12)
    ax.set_title(
        f"Variant: {row['chrom']}:{row['pos']} {row['ref']}→{row['alt']} ({row['label']})",
        fontsize=14,
        fontweight="bold",
    )
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

    # Print top 3
    print(f"🔬 Variant: {row['chrom']}:{row['pos']} {row['ref']}→{row['alt']}")
    print(f"   Label: {row['label']}")
    print("\n   Top 3 Disruptions:")
    for i, (name, val) in enumerate(sorted_deltas[:3], 1):
        direction = "Gain" if val > 0 else "Loss"
        print(f"   {i}. {name}: {val:+.4f} ({direction})")


def plot_variant_fingerprint_normalized(
    df,
    variant_index=0,
    label_filter=None,
    feature_set="bed",
    top_k=15,
    normalize_bw=True,
    figsize=(12, 6),
):
    """
    Plot variant fingerprint with proper BW normalization.

    Parameters:
    -----------
    df : DataFrame
        Results dataframe
    variant_index : int
        Row index within the filtered subset
    label_filter : str or None
        Filter by label ("Pathogenic" or "Benign")
    feature_set : str
        'bed', 'bw', or 'all'
    top_k : int
        Number of features to display
    normalize_bw : bool
        If True, z-score normalize BW tracks for display
    figsize : tuple
        Figure size
    """
    delta_bed_cols = [c for c in df.columns if c.startswith("D_BED_")]
    delta_bw_cols = [c for c in df.columns if c.startswith("D_BW_")]

    # Select columns based on feature_set
    if feature_set == "bed":
        delta_cols = delta_bed_cols
        title_suffix = "BED Elements"
        xlabel = "ΔProbability"
    elif feature_set == "bw":
        delta_cols = delta_bw_cols
        title_suffix = "BigWig Tracks" + (
            " (Z-Normalized)" if normalize_bw else " (Raw Logits)"
        )
        xlabel = "Z-Score" if normalize_bw else "ΔLogit"
    else:
        delta_cols = delta_bed_cols + delta_bw_cols
        title_suffix = "All Features"
        xlabel = "Normalized Delta"

    if not delta_cols:
        print(f"⚠️ No columns found for feature_set='{feature_set}'")
        return

    # Filter by label if specified
    subset = df[df["label"] == label_filter] if label_filter else df
    if variant_index >= len(subset):
        print(f"⚠️ Index {variant_index} out of range")
        return

    target = subset.iloc[variant_index]

    # Extract and normalize deltas
    deltas = {}
    for col in delta_cols:
        val = target[col]
        # Ensure val is a scalar
        if isinstance(val, pd.Series):
            val = float(val.iloc[0])

        # Z-score normalize BW if requested
        if normalize_bw and col in delta_bw_cols:
            benign_mask = df["label"] == False
            mu = df.loc[benign_mask, col].mean()
            sigma = df.loc[benign_mask, col].std()
            val = (val - mu) / sigma if sigma > 0 else 0

        name = col.replace("D_BED_", "").replace("D_BW_", "")
        if col in delta_bw_cols:
            name = "BW:" + name[:30]
        else:
            name = "BED:" + name
        deltas[name] = float(val)

    # Sort by absolute magnitude
    sorted_deltas = sorted(deltas.items(), key=lambda x: abs(x[1]), reverse=True)[
        :top_k
    ]

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    names = [x[0] for x in sorted_deltas]
    vals = [x[1] for x in sorted_deltas]
    colors = ["#d73027" if v > 0 else "#4575b4" for v in vals]

    ax.barh(names, vals, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel(xlabel, fontsize=12)

    label_str = "Pathogenic" if target["label"] else "Benign"
    ax.set_title(
        f"{target['chrom']}:{target['pos']} {target['ref']}→{target['alt']} ({label_str}) - {title_suffix}",
        fontsize=14,
        fontweight="bold",
    )
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()


def plot_combined_fingerprint(df, variant_index=0, top_k_each=8, figsize=(16, 6)):
    """
    Plot BED and BigWig side-by-side for complete mechanistic view.

    Parameters:
    -----------
    df : DataFrame
        Results dataframe
    variant_index : int
        Row index of variant
    top_k_each : int
        Number of features to show in each panel
    figsize : tuple
        Figure size
    """
    row = df.iloc[variant_index]

    delta_bed_cols = [c for c in df.columns if c.startswith("D_BED_")]
    delta_bw_cols = [c for c in df.columns if c.startswith("D_BW_")]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # LEFT: BED Elements
    if delta_bed_cols:
        bed_deltas = {
            col.replace("D_BED_", ""): float(row[col]) for col in delta_bed_cols
        }
        sorted_bed = sorted(bed_deltas.items(), key=lambda x: abs(x[1]), reverse=True)[
            :top_k_each
        ]

        names = [x[0] for x in sorted_bed]
        vals = [x[1] for x in sorted_bed]
        colors = ["#d73027" if v > 0 else "#4575b4" for v in vals]

        ax1.barh(names, vals, color=colors)
        ax1.axvline(0, color="black", linewidth=0.8)
        ax1.set_xlabel("ΔProbability", fontsize=11)
        ax1.set_title(
            "BED Elements (Genomic Structure)", fontsize=12, fontweight="bold"
        )
        ax1.invert_yaxis()

    # RIGHT: BigWig Tracks (Z-normalized)
    if delta_bw_cols:
        benign_mask = df["label"] == False
        bw_deltas = {}

        for col in delta_bw_cols:
            val = float(row[col])
            # Z-normalize
            mu = df.loc[benign_mask, col].mean()
            sigma = df.loc[benign_mask, col].std()
            z_val = (val - mu) / sigma if sigma > 0 else 0

            name = col.replace("D_BW_", "")[:35]
            bw_deltas[name] = z_val

        sorted_bw = sorted(bw_deltas.items(), key=lambda x: abs(x[1]), reverse=True)[
            :top_k_each
        ]

        names = [x[0] for x in sorted_bw]
        vals = [x[1] for x in sorted_bw]
        colors = ["#d73027" if v > 0 else "#4575b4" for v in vals]

        ax2.barh(names, vals, color=colors)
        ax2.axvline(0, color="black", linewidth=0.8)
        ax2.set_xlabel("Z-Score", fontsize=11)
        ax2.set_title(
            "BigWig Tracks (Functional Signals)", fontsize=12, fontweight="bold"
        )
        ax2.invert_yaxis()

    label_str = "Pathogenic" if row["label"] else "Benign"
    fig.suptitle(
        f"Variant: {row['chrom']}:{row['pos']} {row['ref']}→{row['alt']} ({label_str})",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )

    plt.tight_layout()
    plt.show()


def plot_case_study_variant(df, chrom, pos, ref, alt, top_k=12, desc=""):
    """Safely plot a variant if it exists in the dataset."""
    target_variants = df[
        (df["chrom"] == chrom)
        & (df["pos"] == pos)
        & (df["ref"] == ref)
        & (df["alt"] == alt)
    ]
    if len(target_variants) > 0:
        target_idx = target_variants.index[0]
        plot_variant_fingerprint(df.loc[[target_idx]], variant_index=0, top_k=top_k)
    else:
        print(f"⚠️ Variant {chrom}:{pos} {ref}>{alt} not found in dataset")
        if desc:
            print(f"   ({desc})")


# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================


def compare_pathogenic_benign(df, score_type="Impact_Score_Normalized"):
    """
    Statistical comparison of pathogenic vs benign variants.

    Returns:
    --------
    dict with keys: 'mean_path', 'mean_ben', 'cohens_d', 'p_value', 'auc_proxy'
    """
    path = df[df["label"] == True][score_type]
    benign = df[df["label"] == False][score_type]

    mean_path = path.mean()
    mean_ben = benign.mean()

    # Cohen's d
    pooled_std = np.sqrt(
        ((len(path) - 1) * path.std() ** 2 + (len(benign) - 1) * benign.std() ** 2)
        / (len(path) + len(benign) - 2)
    )
    cohens_d = (mean_path - mean_ben) / pooled_std if pooled_std > 0 else 0

    # T-test
    t_stat, p_value = stats.ttest_ind(path, benign)

    # AUC proxy (based on Mann-Whitney U)
    u_stat, _ = stats.mannwhitneyu(path, benign, alternative="two-sided")
    auc_proxy = u_stat / (len(path) * len(benign))

    return {
        "mean_path": mean_path,
        "mean_ben": mean_ben,
        "cohens_d": cohens_d,
        "p_value": p_value,
        "auc_proxy": auc_proxy,
    }


def analyze_per_element_differences(df):
    """
    Analyze which BED elements show significant differences between pathogenic and benign.
    """
    delta_bed_cols = [c for c in df.columns if c.startswith("D_BED_")]

    if not delta_bed_cols:
        print("⚠️ No BED delta columns found")
        return None

    results = []
    path_df = df[df["label"] == True]
    benign_df = df[df["label"] == False]

    for col in delta_bed_cols:
        element = col.replace("D_BED_", "")

        path_vals = np.abs(path_df[col])
        benign_vals = np.abs(benign_df[col])

        path_mean = path_vals.mean()
        benign_mean = benign_vals.mean()

        # Cohen's d
        pooled_std = np.sqrt(
            (
                (len(path_vals) - 1) * path_vals.std() ** 2
                + (len(benign_vals) - 1) * benign_vals.std() ** 2
            )
            / (len(path_vals) + len(benign_vals) - 2)
        )
        cohens_d = (path_mean - benign_mean) / pooled_std if pooled_std > 0 else 0

        # T-test
        t_stat, p_val = stats.ttest_ind(path_vals, benign_vals)

        results.append(
            {
                "Element": element,
                "Pathogenic_Mean": path_mean,
                "Benign_Mean": benign_mean,
                "Difference_Benign_minus_Path": benign_mean - path_mean,
                "Cohens_d": cohens_d,
                "P_Value": p_val,
                "Significant": p_val
                < (0.05 / len(delta_bed_cols)),  # Bonferroni correction
            }
        )

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("Difference_Benign_minus_Path", ascending=False)

    return results_df


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================


def main():
    """Main analysis pipeline."""
    import argparse

    parser = argparse.ArgumentParser(description="Analyze NTv3 inference results")
    parser.add_argument(
        "--input",
        default="clinvar_new_deltas.parquet",
        help="Input parquet file with inference results",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Run in debug mode (first 70 variants)"
    )
    args = parser.parse_args()

    print("=" * 80)
    print("NTv3 Genomic Variant Analysis Pipeline")
    print("=" * 80)

    # Load data
    print(f"\n📂 Loading data from {args.input}...")
    df = pd.read_parquet(args.input)
    print(f"✅ Loaded {len(df)} variants")
    print(f"   Shape: {df.shape}")

    if args.debug:
        df = df.head(70)
        print(f"🐞 DEBUG MODE: Analyzing first {len(df)} variants")

    # Compute normalized scores
    print("\n🔢 Computing normalized impact scores...")
    df = compute_normalized_impact_scores(df, normalization="z_score")

    # Statistical comparison
    print("\n📊 Pathogenic vs Benign Comparison:")
    print("-" * 80)

    score_types = ["Impact_Score_BED", "Impact_Score_BW", "Impact_Score_Normalized"]
    for score_type in score_types:
        if score_type in df.columns:
            stats_dict = compare_pathogenic_benign(df, score_type)
            print(f"\n{score_type}:")
            print(f"   Pathogenic mean: {stats_dict['mean_path']:.4f}")
            print(f"   Benign mean:     {stats_dict['mean_ben']:.4f}")
            print(f"   Cohen's d:       {stats_dict['cohens_d']:.3f}")
            print(f"   P-value:         {stats_dict['p_value']:.2e}")
            print(f"   AUC (proxy):     {stats_dict['auc_proxy']:.3f}")

    # Per-element analysis
    print("\n\n📊 Per-Element Pathogenic vs Benign Comparison:")
    print("-" * 80)
    element_df = analyze_per_element_differences(df)

    if element_df is not None:
        print("\nTop elements with higher disruption in Benign vs Pathogenic variants:")
        bonferroni_alpha = 0.05 / len(element_df)
        print(
            f"(Bonferroni-corrected significance threshold: p < {bonferroni_alpha:.4f})"
        )
        print()
        print(element_df.head(10).to_string(index=False))

        print("\n📈 Summary Statistics:")
        print(f"   N Pathogenic: {len(df[df['label'] == True])}")
        print(f"   N Benign: {len(df[df['label'] == False])}")
        print(f"   Bonferroni-corrected α: {bonferroni_alpha:.4f}")
        print(
            f"   Significant elements: {element_df['Significant'].sum()}/{len(element_df)}"
        )

        cohens_d_vals = element_df["Cohens_d"].abs()
        print("\n   Effect sizes (Cohen's d):")
        print(f"   - Median |d|: {cohens_d_vals.median():.3f}")
        print(
            f"   - Range: [{element_df['Cohens_d'].min():.3f}, {element_df['Cohens_d'].max():.3f}]"
        )

    # Example visualization
    print("\n\n" + "=" * 80)
    print("Example usage:")
    print("=" * 80)

    if len(df) > 0:
        # Show first pathogenic variant
        path_variants = df[df["label"] == True]
        if len(path_variants) > 0:
            idx = path_variants.index[0]
            row = df.loc[idx]
            print(f"🔬 Variant: {row['chrom']}:{row['pos']} {row['ref']}→{row['alt']}")
            print(f"   Label: {'Pathogenic' if row['label'] else 'Benign'}")
            if "Impact_Score_Normalized" in row:
                print(f"   Impact Score: {row['Impact_Score_Normalized']:.4f}")

            # Show top disruptions
            delta_bed_cols = [c for c in df.columns if c.startswith("D_BED_")]
            if delta_bed_cols:
                deltas = [
                    (col.replace("D_BED_", ""), abs(row[col])) for col in delta_bed_cols
                ]
                deltas.sort(key=lambda x: x[1], reverse=True)

                print("\n   Top 3 Disruptions:")
                for i, (name, val) in enumerate(deltas[:3], 1):
                    sign = "+" if row[f"D_BED_{name}"] > 0 else "-"
                    direction = "Gain" if sign == "+" else "Loss"
                    print(f"   {i}. {name}: {sign}{val:.4f} ({direction})")

    print("\n" + "=" * 80)
    print("✅ Analysis complete!")
    print("=" * 80)

    return df


if __name__ == "__main__":
    main()
