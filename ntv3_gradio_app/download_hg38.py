#!/usr/bin/env python3
"""
Download and extract hg38 reference genome
==========================================
Downloads hg38.fa.gz from UCSC and extracts it to data/hg38.fa
Only runs if hg38.fa doesn't already exist (3GB download).

Usage:
    python download_hg38.py
"""
import os
import urllib.request
import gzip
import shutil
from pathlib import Path

# Define paths relative to this script
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
FASTA_PATH = DATA_DIR / "hg38.fa"
GZ_PATH = DATA_DIR / "hg38.fa.gz"
URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"


def download_and_extract():
    """Download and extract hg38 reference genome."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if FASTA_PATH.exists():
        size_mb = FASTA_PATH.stat().st_size / (1024 * 1024)
        print(f"✅ '{FASTA_PATH}' already exists ({size_mb:.0f} MB). Skipping download.")
        return

    print("⬇️  Downloading hg38.fa.gz from UCSC (~3 GB)...")
    print(f"   URL: {URL}")
    print("   This may take 10-30 minutes depending on your connection.\n")

    # Download with progress reporting
    def reporthook(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            percent = min(100, downloaded * 100 / total_size)
            mb_downloaded = downloaded / (1024 * 1024)
            mb_total = total_size / (1024 * 1024)
            print(
                f"\r   Progress: {percent:.1f}% ({mb_downloaded:.0f}/{mb_total:.0f} MB)",
                end="",
            )

    urllib.request.urlretrieve(URL, GZ_PATH, reporthook=reporthook)
    print("\n")

    print(f"📦 Extracting to {FASTA_PATH}...")
    with gzip.open(GZ_PATH, "rb") as f_in:
        with open(FASTA_PATH, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    print("🧹 Cleaning up compressed file...")
    GZ_PATH.unlink()

    size_mb = FASTA_PATH.stat().st_size / (1024 * 1024)
    print(f"✅ Done! hg38.fa extracted ({size_mb:.0f} MB)\n")


if __name__ == "__main__":
    download_and_extract()
