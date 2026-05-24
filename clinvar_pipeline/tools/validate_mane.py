#!/usr/bin/env python3
"""Validate MANE/Promoter metadata: local build vs parent reference."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from clinvar.mane_build import build_mane_tables
from config import DATA_ROOT, PATHS, PROJECT_ROOT


def _feature_counts(mane: pd.DataFrame) -> dict[str, int]:
    return mane["Feature"].value_counts().to_dict()


def validate(*, rebuild: bool = False) -> dict:
    report: dict = {"status": "OK", "checks": []}

    def check(name: str, ok: bool, detail: str) -> None:
        report["checks"].append({"name": name, "ok": bool(ok), "detail": str(detail)})
        if not ok:
            report["status"] = "FAIL"

    gff = Path(PATHS["MANE_gff"])
    local_mane = Path(PATHS["MANE_processed"])
    local_prom = Path(PATHS["Promoter_processed"])
    parent_mane = PROJECT_ROOT / "metadata/MANE_processed.csv"
    parent_prom = PROJECT_ROOT / "metadata/Promoter_processed.csv"

    check("MANE_gff exists", gff.exists(), str(gff))
    check("local MANE_processed exists", local_mane.exists(), str(local_mane))
    check("local Promoter_processed exists", local_prom.exists(), str(local_prom))

    if not local_mane.exists() or not local_prom.exists():
        return report

    mane_local = pd.read_csv(local_mane)
    prom_local = pd.read_csv(local_prom)
    check("MANE rows > 0", len(mane_local) > 0, f"{len(mane_local):,} rows")
    check("Promoter rows > 0", len(prom_local) > 0, f"{len(prom_local):,} rows")
    check(
        "Promoter has coordinates",
        {"Promoter_Start", "Promoter_End"}.issubset(prom_local.columns),
        str(list(prom_local.columns[:8])),
    )

    if parent_mane.exists():
        mane_parent = pd.read_csv(parent_mane)
        check(
            "MANE row count vs parent",
            len(mane_local) == len(mane_parent),
            f"local={len(mane_local):,} parent={len(mane_parent):,}",
        )
        check(
            "MANE columns vs parent",
            set(mane_local.columns) == set(mane_parent.columns),
            f"Δ cols={set(mane_local.columns) ^ set(mane_parent.columns)}",
        )
        if len(mane_local) == len(mane_parent):
            # spot-check first/last transcript_id
            for col in ("transcript_id", "Feature", "Chromosome"):
                if col in mane_local.columns:
                    same = float(
                        (mane_local[col].astype(str).values == mane_parent[col].astype(str).values).mean()
                    )
                    check(f"MANE {col} match rate", same > 0.999, f"{same:.4%}")

    if parent_prom.exists():
        prom_parent = pd.read_csv(parent_prom)
        check(
            "Promoter row count vs parent",
            len(prom_local) == len(prom_parent),
            f"local={len(prom_local):,} parent={len(prom_parent):,}",
        )

    if rebuild and gff.exists():
        mane_built, prom_built = build_mane_tables(gff)
        check(
            "rebuild MANE rows vs local file",
            len(mane_built) == len(mane_local),
            f"built={len(mane_built):,} local={len(mane_local):,}",
        )
        check(
            "rebuild Promoter rows vs local file",
            len(prom_built) == len(prom_local),
            f"built={len(prom_built):,} local={len(prom_local):,}",
        )
        fc_built = _feature_counts(mane_built)
        fc_local = _feature_counts(mane_local)
        check(
            "rebuild MANE feature mix",
            fc_built.get("mRNA", 0) == fc_local.get("mRNA", 0),
            f"mRNA built={fc_built.get('mRNA')} local={fc_local.get('mRNA')}",
        )

    return report


def main() -> None:
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true", help="Also validate GFF → table rebuild")
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()

    report = validate(rebuild=args.rebuild)
    print(f"MANE validation: {report['status']}")
    for c in report["checks"]:
        mark = "✓" if c["ok"] else "✗"
        print(f"  {mark} {c['name']}: {c['detail']}")

    if args.write_report:
        out = Path(PATHS["validation_root"])
        out.mkdir(parents=True, exist_ok=True)
        (out / "mane_validation.json").write_text(json.dumps(report, indent=2))
        print(f"Report → {out / 'mane_validation.json'}")

    sys.exit(0 if report["status"] == "OK" else 1)


if __name__ == "__main__":
    main()
