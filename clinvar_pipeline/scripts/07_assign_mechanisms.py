#!/usr/bin/env python3
"""Assign NT signal mechanisms (uses vendor/ or repo-root nt_mechanism_assignment.py)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PIPELINE_ROOT.parent
VENDOR_ROOT = PIPELINE_ROOT / "vendor"
sys.path.insert(0, str(PIPELINE_ROOT))
if (VENDOR_ROOT / "nt_mechanism_assignment.py").exists():
    sys.path.insert(0, str(VENDOR_ROOT))
elif str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))

from config import PATHS
from nt_mechanism_assignment import run_assignment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    Path(PATHS["snp_mechanisms_csv"]).parent.mkdir(parents=True, exist_ok=True)
    run_assignment(
        PATHS["snp_annotated_signaled"],
        PATHS["snp_mechanisms_csv"],
        "snp",
    )
    run_assignment(
        PATHS["indel_annotated_signaled"],
        PATHS["indel_mechanisms_csv"],
        "indel",
    )
    print("\n✓ Mechanism assignment complete.")


if __name__ == "__main__":
    main()
