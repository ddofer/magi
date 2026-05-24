#!/usr/bin/env python3
"""Build MANE_processed.csv and Promoter_processed.csv from MANE RefSeq GFF."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import DATA_ROOT, PATHS
from clinvar.mane_build import write_mane_metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gff",
        default=PATHS["MANE_gff"],
        help="MANE GRCh38 RefSeq genomic GFF (.gz)",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DATA_ROOT / "metadata"),
        help="Output directory for processed CSVs",
    )
    args = parser.parse_args()

    out = Path(args.out_dir)
    write_mane_metadata(
        args.gff,
        out / "MANE_processed.csv",
        out / "Promoter_processed.csv",
    )
    print("\n✓ MANE metadata complete.")


if __name__ == "__main__":
    main()
