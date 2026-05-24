#!/usr/bin/env python3
"""Compute quality vs rationaled cohort counts → output/validation/cohort_funnel.json."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import PATHS
from clinvar.variant_prep import measure_cohort_funnel, write_cohort_funnel


def main() -> None:
    funnel = measure_cohort_funnel(PATHS)
    write_cohort_funnel(PATHS, funnel)
    for vt, stats in funnel.items():
        print(f"  {vt}: quality={stats['quality']:,}  rationaled={stats['rationaled']:,}")


if __name__ == "__main__":
    main()
