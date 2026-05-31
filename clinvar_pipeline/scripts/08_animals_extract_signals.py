#!/usr/bin/env python3
"""Stage 8: OMIA animal deltas → top-k signal summaries (signaled parquets)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import DEFAULT_ANIMALS_TOP_K_SIGNALS, PATHS
from clinvar.animals_pipeline import extract_animals_signals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant-type", choices=["snp", "indel"], required=True)
    parser.add_argument(
        "--top-k",
        type=int,
        default=DEFAULT_ANIMALS_TOP_K_SIGNALS,
        help=f"Top signals per category (default {DEFAULT_ANIMALS_TOP_K_SIGNALS})",
    )
    args = parser.parse_args()
    extract_animals_signals(PATHS, args.variant_type, k_signals=args.top_k)


if __name__ == "__main__":
    main()
