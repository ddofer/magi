#!/usr/bin/env python3
"""Compare refactored pipeline outputs against legacy notebook parquets."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import PATHS
from clinvar.prompts import build_prompts_table

# (label, legacy_path_key, new_path_key, key_cols, stage_kind)
STAGE_PAIRS = [
    (
        "01 strict snp",
        "legacy_snp_strict",
        "snp_strict",
        ["#VariationID", "chrom", "pos", "ref", "alt"],
        "parquet",
    ),
    (
        "01 strict indel",
        "legacy_indel_strict",
        "indel_strict",
        ["#VariationID", "chrom", "pos", "ref", "alt"],
        "parquet",
    ),
    (
        "01 annotated snp",
        "legacy_snp_annotated",
        "snp_annotated",
        ["#VariationID"],
        "parquet",
    ),
    (
        "01 annotated indel",
        "legacy_indel_annotated",
        "indel_annotated",
        ["#VariationID"],
        "parquet",
    ),
    (
        "02 signaled snp",
        "legacy_snp_signaled",
        "snp_annotated_signaled",
        ["#VariationID"],
        "parquet",
    ),
    (
        "02 signaled indel",
        "legacy_indel_signaled",
        "indel_annotated_signaled",
        ["#VariationID"],
        "parquet",
    ),
    (
        "03 prompts snp",
        "legacy_snp_signaled",
        "snp_with_prompts",
        ["#VariationID"],
        "prompts",
    ),
    (
        "03 prompts indel",
        "legacy_indel_signaled",
        "indel_with_prompts",
        ["#VariationID"],
        "prompts",
    ),
]


def _load(path: str) -> pd.DataFrame | None:
    p = Path(path)
    if not p.exists():
        return None
    return pd.read_parquet(p)


def compare_parquet_pair(
    label: str,
    legacy_path: str,
    new_path: str,
    key_cols: list[str],
) -> dict:
    legacy = _load(legacy_path)
    new = _load(new_path)

    if legacy is None and new is None:
        return {"label": label, "status": "SKIP", "detail": "both missing"}
    if legacy is None:
        return {"label": label, "status": "MISSING_LEGACY", "detail": legacy_path}
    if new is None:
        return {"label": label, "status": "MISSING_NEW", "detail": new_path}

    report = {
        "label": label,
        "legacy_rows": len(legacy),
        "new_rows": len(new),
        "legacy_cols": len(legacy.columns),
        "new_cols": len(new.columns),
    }

    keys = [c for c in key_cols if c in legacy.columns and c in new.columns]
    if not keys:
        report["status"] = "NO_KEY"
        return report

    legacy_keys = legacy[keys].astype(str).drop_duplicates()
    new_keys = new[keys].astype(str).drop_duplicates()
    merged_keys = legacy_keys.merge(new_keys, on=keys, how="outer", indicator=True)

    only_legacy = int((merged_keys["_merge"] == "left_only").sum())
    only_new = int((merged_keys["_merge"] == "right_only").sum())
    shared = int((merged_keys["_merge"] == "both").sum())

    report.update(
        {
            "shared_keys": shared,
            "only_legacy": only_legacy,
            "only_new": only_new,
            "status": (
                "OK"
                if only_legacy == 0 and only_new == 0 and len(legacy) == len(new)
                else "DIFF"
            ),
        }
    )

    skip_cols = {"transcript_set", "promoter_transcript_set"}
    overlap_cols = sorted(set(legacy.columns) & set(new.columns) - set(keys) - skip_cols)
    if shared and overlap_cols:
        l = legacy.copy()
        n = new.copy()
        for col in keys:
            l[col] = l[col].astype(str)
            n[col] = n[col].astype(str)
        l = l.merge(new_keys.astype(str), on=keys).set_index(keys)
        n = n.merge(new_keys.astype(str), on=keys).set_index(keys)
        common_idx = l.index.intersection(n.index)
        mismatches = []
        for col in overlap_cols:
            a = l.loc[common_idx, col]
            b = n.loc[common_idx, col]
            if a.dtype == object or b.dtype == object:
                ne = int((a.astype(str) != b.astype(str)).sum())
            else:
                ne = int(((a != b) & ~(a.isna() & b.isna())).sum())
            if ne:
                mismatches.append((col, ne))
        mismatches.sort(key=lambda x: -x[1])
        report["column_mismatches"] = mismatches[:15]
        report["total_mismatch_cells"] = int(sum(m for _, m in mismatches))

    return report


def compare_prompts_pair(
    label: str,
    legacy_signaled_path: str,
    new_prompts_path: str,
    variant_type: str,
) -> dict:
    legacy_signaled = _load(legacy_signaled_path)
    new_prompts = _load(new_prompts_path)

    if legacy_signaled is None:
        return {"label": label, "status": "MISSING_LEGACY", "detail": legacy_signaled_path}
    if new_prompts is None:
        return {"label": label, "status": "MISSING_NEW", "detail": new_prompts_path}

    legacy_prompts = build_prompts_table(legacy_signaled, variant_type=variant_type)

    key = "#VariationID"
    legacy_prompts[key] = legacy_prompts[key].astype(str)
    new_prompts[key] = new_prompts[key].astype(str)

    merged = legacy_prompts[[key, "llm_prompt"]].merge(
        new_prompts[[key, "llm_prompt"]],
        on=key,
        suffixes=("_legacy", "_new"),
    )

    prompt_mismatch = int((merged["llm_prompt_legacy"] != merged["llm_prompt_new"]).sum())

    return {
        "label": label,
        "status": "OK" if prompt_mismatch == 0 else "DIFF",
        "variants": len(merged),
        "prompt_mismatch": prompt_mismatch,
        "match_rate_pct": round(100 * (1 - prompt_mismatch / len(merged)), 2) if len(merged) else 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=[s[0] for s in STAGE_PAIRS],
        default=None,
        help="Compare one stage only",
    )
    parser.add_argument(
        "--prefix",
        choices=["01", "02", "03", "all"],
        default="all",
        help="Compare only stages 01, 02, or 03",
    )
    args = parser.parse_args()

    pairs = STAGE_PAIRS
    if args.stage:
        pairs = [p for p in STAGE_PAIRS if p[0] == args.stage]
    elif args.prefix != "all":
        pairs = [p for p in STAGE_PAIRS if p[0].startswith(args.prefix)]

    print("Comparing legacy (parent repo) vs refactored (clinvar_pipeline/output)\n")
    for label, legacy_key, new_key, keys, kind in pairs:
        if kind == "prompts":
            variant_type = "snp" if "snp" in label else "indel"
            rep = compare_prompts_pair(label, PATHS[legacy_key], PATHS[new_key], variant_type)
        else:
            rep = compare_parquet_pair(label, PATHS[legacy_key], PATHS[new_key], keys)
        print(f"── {label} ──")
        for k, v in rep.items():
            print(f"  {k}: {v}")
        print()


if __name__ == "__main__":
    main()
