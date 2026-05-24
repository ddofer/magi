#!/usr/bin/env python3
"""Stage 2: extract top-k BED/BW/MLM signals from annotated parquets."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import DEFAULT_TOP_K_SIGNALS, PATHS
from clinvar.pipeline import extract_signals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant-type",
        choices=["snp", "indel"],
        required=True,
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=DEFAULT_TOP_K_SIGNALS,
        help=f"Top signals per category (default {DEFAULT_TOP_K_SIGNALS})",
    )
    args = parser.parse_args()
    extract_signals(PATHS, args.variant_type, k_signals=args.top_k)


if __name__ == "__main__":
    main()
