#!/usr/bin/env python3
"""Stage 1: NT deltas → ClinVar filter → MANE annotation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import PATHS
from clinvar.pipeline import prepare_variants


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snp", action="store_true", help="Process SNPs")
    parser.add_argument("--indel", action="store_true", help="Process indels")
    args = parser.parse_args()

    run_snp = args.snp or not (args.snp or args.indel)
    run_indel = args.indel or not (args.snp or args.indel)
    prepare_variants(PATHS, run_snp=run_snp, run_indel=run_indel)


if __name__ == "__main__":
    main()
