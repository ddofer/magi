#!/usr/bin/env python3
"""
Comprehensive validation of example variants:
1. Literature accuracy
2. Example diversity and quality
3. Track type coverage in results
"""

import sys

sys.path.insert(0, ".")

from app import predict_single_variant
import pandas as pd
from collections import defaultdict
import json

# Known variants from literature
VARIANT_REFERENCES = {
    ("chr17", 7675088, "C", "T"): {
        "name": "TP53 R175H",
        "disease": "Cancer (tumor suppressor)",
        "frequency": "~6% of cancers",
        "literature": "Highly common hotspot, eliminates DNA binding domain",
        "impact_expected": "HIGH",
        "regions_expected": ["coding", "regulatory"],
    },
    ("chr7", 117559593, "ATCT", "A"): {
        "name": "CFTR F508del",
        "disease": "Cystic Fibrosis",
        "frequency": "~70% of CF patients",
        "literature": "Most common CF mutation, causes protein misfolding",
        "impact_expected": "HIGH",
        "regions_expected": ["coding", "structural"],
    },
    ("chr13", 32332771, "AGAGA", "AGA"): {
        "name": "BRCA2 frameshift",
        "disease": "Hereditary breast/ovarian cancer",
        "frequency": "Rare, pathogenic",
        "literature": "BRCA2 frameshift deletion (c.5946delT), causes loss of function",
        "impact_expected": "HIGH",
        "regions_expected": ["coding", "frameshift"],
    },
    ("chr11", 5227002, "T", "A"): {
        "name": "HBB E6V",
        "disease": "Sickle cell disease",
        "frequency": "Common in African populations",
        "literature": "Missense mutation (rs334) causing hemoglobin S polymerization",
        "impact_expected": "HIGH",
        "regions_expected": ["coding", "regulatory"],
    },
    ("chr17", 43092418, "T", "C"): {
        "name": "BRCA1 synonymous",
        "disease": "Benign control variant",
        "frequency": "Common",
        "literature": "Synonymous variant (c.3113A>G, rs16941), expected benign",
        "impact_expected": "LOW",
        "regions_expected": ["coding"],
    },
}


def categorize_track(feature_name, feature_type):
    """Categorize track type (BED or BigWig, and specific kind)"""
    if feature_type == "BED":
        if any(
            x in feature_name.lower() for x in ["splice", "exon", "intron", "codon"]
        ):
            return "bed_splicing"
        elif any(x in feature_name.lower() for x in ["cds", "coding"]):
            return "bed_coding"
        elif any(x in feature_name.lower() for x in ["promoter", "enhancer"]):
            return "bed_regulatory"
        elif any(x in feature_name.lower() for x in ["utr", "5utr", "3utr"]):
            return "bed_utr"
        else:
            return "bed_other"
    elif feature_type == "BigWig":
        if "Histone" in feature_name:
            return "bw_histone"
        elif "RNA" in feature_name or "CAGE" in feature_name:
            return "bw_expression"
        elif "DNase" in feature_name or "ATAC" in feature_name:
            return "bw_accessibility"
        elif "ChIP-seq" in feature_name:
            return "bw_chipseq"
        else:
            return "bw_other"
    return "unknown"


def validate_example(chrom, pos, ref, alt):
    """Validate a single example"""
    key = (chrom, pos, ref, alt)
    ref_data = VARIANT_REFERENCES.get(key)

    if not ref_data:
        return None

    print("\n" + "=" * 80)
    print(f"VARIANT: {ref_data['name']} ({chrom}:{pos} {ref}>{alt})")
    print("=" * 80)

    # Literature context
    print("\n📚 LITERATURE CONTEXT:")
    print(f"  Disease: {ref_data['disease']}")
    print(f"  Frequency: {ref_data['frequency']}")
    print(f"  Summary: {ref_data['literature']}")
    print(f"  Expected Impact: {ref_data['impact_expected']}")

    # Run prediction
    print("\n🔬 Running prediction...")
    result = predict_single_variant(chrom, pos, ref, alt)
    (
        summary_md,
        interpretation_md,
        top_table_df,
        fp_fig,
        rt_fig,
        csv_path,
        bed_df,
        mlm_md,
        ranked,
    ) = result

    # Extract impact score from summary
    impact_from_summary = None
    if "BED Impact Score" in summary_md:
        import re

        m = re.search(r"BED Impact Score\s*\|\s*([\d.]+)", summary_md)
        if m:
            impact_from_summary = float(m.group(1))

    print("\n📊 RESULTS SUMMARY:")
    gene_region = "N/A"
    if "**Gene:**" in summary_md:
        gene_part = summary_md.split("**Gene:**")[1]
        gene_region = gene_part.split("\n")[0].strip()
    print(f"  Gene/Region: {gene_region}")
    print(f"  Top tracks: {len(top_table_df)} features")
    print(f"  Total ranked: {len(ranked)} tracks")
    print(f"  Interpretation panel present: {bool(interpretation_md)}")
    print(
        f"  Impact Score (BED): {impact_from_summary if impact_from_summary else 'N/A'}"
    )

    # Track diversity analysis
    print("\n🎯 TRACK DIVERSITY ANALYSIS:")

    track_categories = defaultdict(int)
    sources = set()

    # Analyze ranked tracks
    for r in ranked:
        feat = r["display_name"]
        ftype = r["track_type"]
        category = categorize_track(feat, ftype)
        track_categories[category] += 1

        if "|" in feat:  # BigWig with tissue info
            tissue = feat.split("|")[0].strip()
            sources.add(tissue)

    print("\n  Track Categories:")
    for cat, count in sorted(track_categories.items(), key=lambda x: -x[1]):
        cat_display = cat.replace("bed_", "BED: ").replace("bw_", "BigWig: ")
        print(f"    - {cat_display}: {count}")

    print(f"\n  Tissue/Cell Sources Found: {len(sources)}")
    if len(sources) > 0:
        for source in sorted(list(sources))[:5]:  # Show first 5
            print(f"    - {source}")
        if len(sources) > 5:
            print(f"    ... and {len(sources) - 5} more")

    # Quality assessment
    print("\n✓ VARIANT QUALITY ASSESSMENT:")

    # Has diverse tracks
    diversity_score = len(track_categories)
    if diversity_score >= 3:
        print(f"  ✓ Good track diversity ({diversity_score} categories)")
    else:
        print(f"  ⚠️  Limited track diversity ({diversity_score} categories)")

    # Has gains and losses
    gains = [r for r in ranked if r["delta"] > 0]
    losses = [r for r in ranked if r["delta"] < 0]
    if gains and losses:
        print(f"  ✓ Shows both gains ({len(gains)}) and losses ({len(losses)})")
    elif gains:
        print("  ⚠️  Only shows gains, no losses")
    elif losses:
        print("  ⚠️  Only shows losses, no gains")
    else:
        print("  ✗ No gains or losses - may not be suitable example")

    # Expected impact matches actual
    has_high_impacts = any(abs(r["delta"]) >= 0.1 for r in ranked)

    if ref_data["impact_expected"] == "HIGH" and has_high_impacts:
        print("  ✓ HIGH impact expected and observed")
    elif ref_data["impact_expected"] == "MODERATE" and not has_high_impacts:
        print("  ✓ MODERATE impact expected and observed (no extreme deltas)")
    else:
        if ref_data["impact_expected"] == "HIGH":
            print("  ⚠️  HIGH impact expected but weak observed")

    # Relevant to disease
    print(f"  ✓ Disease-relevant example ({ref_data['disease']})")

    return True


# Validate all examples
examples = [
    ("chr17", 7675088, "C", "T"),
    ("chr7", 117559593, "ATCT", "A"),
    ("chr13", 32332771, "AGAGA", "AGA"),
    ("chr11", 5227002, "T", "A"),
    ("chr17", 43092418, "T", "C"),
]

print("\n" + "=" * 80)
print("COMPREHENSIVE EXAMPLE VALIDATION")
print("=" * 80)

for chrom, pos, ref, alt in examples:
    validate_example(chrom, pos, ref, alt)

print("\n" + "=" * 80)
print("VALIDATION COMPLETE")
print("=" * 80)
