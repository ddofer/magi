#!/usr/bin/env python3
"""Copy pipeline input files from an external ClinVar data tree into clinvar_pipeline/data/."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import DATA_ROOT, PROJECT_ROOT

# (parent relative, data/ relative)
INPUT_FILES = [
    ("metadata/MANE_processed.csv", "metadata/MANE_processed.csv"),
    ("metadata/Promoter_processed.csv", "metadata/Promoter_processed.csv"),
    ("metadata/functional_tracks_metadata_human.csv", "metadata/functional_tracks_metadata_human.csv"),
    ("metadata/MANE.GRCh38.v1.5.refseq_genomic.gff.gz", "metadata/MANE.GRCh38.v1.5.refseq_genomic.gff.gz"),
    ("clinvar_tables/variant_summary.txt.gz", "clinvar_tables/variant_summary.txt.gz"),
    ("clinvar_tables/submission_summary.txt.gz", "clinvar_tables/submission_summary.txt.gz"),
    ("thresholds_snps.csv", "thresholds_snps.csv"),
    ("thresholds_indels.csv", "thresholds_indels.csv"),
    ("AF/snps_AF_sub.csv", "AF/snps_AF_sub.csv"),
    ("AF/indel_AF_sub.csv", "AF/indel_AF_sub.csv"),
    ("uncertain_ids.csv", "uncertain_ids.csv"),
]

LARGE_PARQUETS = [
    ("parquet/clinvar_new_deltas.parquet", "parquet/clinvar_new_deltas.parquet"),
    ("parquet/clinvar_indel_deltas.parquet", "parquet/clinvar_indel_deltas.parquet"),
]

LARGE_BYTES = 200 * 1024 * 1024  # 200 MB


def _link_or_copy(src: Path, dst: Path, *, symlink_large: bool) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    size = src.stat().st_size
    if symlink_large and size >= LARGE_BYTES:
        os.symlink(src.resolve(), dst)
        return "symlink"
    shutil.copy2(src, dst)
    return "copy"


def sync(
    *,
    source_root: Path | None = None,
    include_parquets: bool = True,
    symlink_large: bool = True,
    force: bool = False,
) -> None:
    root = source_root or PROJECT_ROOT
    todo = list(INPUT_FILES)
    if include_parquets:
        todo.extend(LARGE_PARQUETS)

    copied, linked, skipped, missing = [], [], [], []
    for src_rel, dst_rel in todo:
        src = root / src_rel
        dst = DATA_ROOT / dst_rel
        if not src.exists():
            missing.append(str(src))
            continue
        if dst.exists() and not force:
            skipped.append(str(dst))
            continue
        mode = _link_or_copy(src, dst, symlink_large=symlink_large)
        (linked if mode == "symlink" else copied).append(f"{src_rel} → data/{dst_rel} ({mode})")

    print(f"Synced inputs into {DATA_ROOT}")
    for line in copied + linked:
        print(f"  {line}")
    if skipped:
        print(f"Skipped {len(skipped)} existing (use --force)")
    if missing:
        print(f"Missing {len(missing)}:")
        for m in missing:
            print(f"  {m}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help="ClinVar data root (default: parent of clinvar_pipeline/)",
    )
    parser.add_argument("--no-parquets", action="store_true", help="Skip large delta parquets")
    parser.add_argument("--copy-large", action="store_true", help="Copy large files instead of symlinking")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    sync(
        source_root=args.source,
        include_parquets=not args.no_parquets,
        symlink_large=not args.copy_large,
        force=args.force,
    )


if __name__ == "__main__":
    main()
