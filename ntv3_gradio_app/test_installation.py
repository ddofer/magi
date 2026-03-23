#!/usr/bin/env python3
"""
Quick validation script for the MAGI Gradio app.

The local genome is optional. If `data/hg38.fa` is absent, the app can fall back
to UCSC sequence retrieval.
"""

import sys
import os
import requests


def check_dependencies():
    """Check if all required packages are importable"""
    print("🔍 Checking dependencies...")
    required = [
        ("torch", "PyTorch"),
        ("transformers", "Transformers"),
        ("gradio", "Gradio"),
        ("pandas", "Pandas"),
        ("numpy", "NumPy"),
        ("matplotlib", "Matplotlib"),
        ("pyfaidx", "pyfaidx"),
    ]

    missing = []
    for module, name in required:
        try:
            __import__(module)
            print(f"  ✓ {name}")
        except ImportError:
            print(f"  ✗ {name} (missing)")
            missing.append(module)

    if missing:
        print(f"\n❌ Missing packages: {', '.join(missing)}")
        print(f"   Install with: pip install {' '.join(missing)}")
        return False
    else:
        print("\n✅ All dependencies installed\n")
        return True


def check_data_files():
    """Check if required data files exist (local genome optional)."""
    print("🔍 Checking data files...")

    required_files = [
        ("data/MANE_processed.csv", "MANE transcripts", 200_000_000),  # ~235 MB
        ("data/Promoter_processed.csv", "Promoter annotations", 10_000_000),  # ~11 MB
        (
            "data/functional_tracks_metadata_human.csv",
            "BigWig metadata",
            500_000,
        ),  # ~519 KB
    ]

    all_present = True
    for file_path, description, min_size in required_files:
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            size_mb = size / 1_000_000
            if size >= min_size:
                print(f"  ✓ {description} ({size_mb:.1f} MB)")
            else:
                print(
                    f"  ⚠ {description} exists but may be corrupted ({size_mb:.1f} MB < {min_size / 1_000_000:.1f} MB expected)"
                )
                all_present = False
        else:
            print(f"  ✗ {description} (missing: {file_path})")
            all_present = False

    local_fa = os.path.exists("data/hg38.fa")
    local_fai = os.path.exists("data/hg38.fa.fai")
    local_fa_gz = os.path.exists("data/hg38.fa.gz")

    if local_fa:
        size_mb = os.path.getsize("data/hg38.fa") / 1_000_000
        print(f"  ✓ Local reference genome available ({size_mb:.1f} MB)")
        if local_fai:
            idx_mb = os.path.getsize("data/hg38.fa.fai") / 1_000_000
            print(f"  ✓ Genome index available ({idx_mb:.2f} MB)")
        else:
            print(
                "  ⚠ Local genome index missing (will auto-build on first local lookup)"
            )
    elif local_fa_gz:
        gz_mb = os.path.getsize("data/hg38.fa.gz") / 1_000_000
        print(
            f"  ⚠ Compressed genome found ({gz_mb:.1f} MB); app will use UCSC API unless decompressed"
        )
    else:
        print("  ⚠ Local genome not found; app will use UCSC API fallback")

    if not all_present:
        print("\n❌ Some data files missing or incomplete")
        return False
    else:
        print("\n✅ Required data files present\n")
        return True


def test_sequence_access():
    """Test sequence retrieval path: local hg38.fa if available, otherwise UCSC API."""
    print("🔍 Testing sequence access...")

    def _test_ucsc_api():
        response = requests.get(
            "https://api.genome.ucsc.edu/getData/sequence",
            params={
                "genome": "hg38",
                "chrom": "chr17",
                "start": 7675087,
                "end": 7675088,
            },
            timeout=10,
        )
        response.raise_for_status()
        dna = str(response.json().get("dna", "")).upper()
        if dna and dna[0] in {"A", "C", "G", "T", "N"}:
            print(f"  ✓ UCSC API reachable (chr17:7675088 -> {dna[0]})")
            print("✅ Sequence access test passed\n")
            return True
        return False

    if not os.path.exists("data/hg38.fa") or not os.path.exists("data/hg38.fa.fai"):
        if os.path.exists("data/hg38.fa") and not os.path.exists("data/hg38.fa.fai"):
            print(
                "  ⚠ Local hg38.fa found but index missing; validating UCSC API path to avoid long index build"
            )
        try:
            if _test_ucsc_api():
                return True
        except Exception as e:
            print(f"  ✗ UCSC API test failed: {e}")
            return False

    try:
        import pyfaidx

        genome = pyfaidx.Fasta("data/hg38.fa")

        # Test accessing chr17 (TP53 region), with the same fallback behavior as the app.
        seq = None
        for test_chrom in ("chr17", "17"):
            if test_chrom in genome.keys():
                seq = str(genome[test_chrom][7675087:7675088]).upper()
                print(f"  ✓ Successfully accessed {test_chrom} locally (TP53 position: {seq})")
                break

        if seq is None:
            print("  ⚠ Local genome could not serve chr17; testing UCSC fallback instead")
            genome.close()
            return _test_ucsc_api()

        genome.close()
        print("✅ Genome access test passed\n")
        return True

    except Exception as e:
        print(f"  ✗ Error accessing genome: {e}")
        return False


def test_model_loading():
    """Test that we can load the configured NTv3 model."""
    print("🔍 Testing model loading (this may take 30-60 seconds)...")

    try:
        from transformers import AutoModel, AutoTokenizer

        model_name = "InstaDeepAI/NTv3_650M_post"
        hf_token = os.environ.get("HF_TOKEN")
        print("  ⏳ Loading tokenizer...")
        tokenizer = AutoTokenizer.from_pretrained(
            model_name,
            trust_remote_code=True,
            token=hf_token,
        )

        print("  ⏳ Loading model (650M configuration)...")
        model = AutoModel.from_pretrained(
            model_name,
            trust_remote_code=True,
            token=hf_token,
        )

        print("  ✓ Model loaded successfully")
        print(f"  ✓ Tokenizer vocab size: {tokenizer.vocab_size}")
        print(
            f"  ✓ Model parameters: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M"
        )

        print("✅ Model loading test passed\n")
        return True

    except Exception as e:
        print(f"  ✗ Error loading model: {e}")
        if "gated" in str(e).lower() or "401" in str(e):
            print("  ⚠ Accept the model terms and set HF_TOKEN before retrying")
        return False


def main():
    """Run all validation tests"""
    print("\n" + "=" * 60)
    print("  MAGI Gradio App - Installation Validation")
    print("=" * 60 + "\n")

    # Change to app directory if needed
    if os.path.exists("ntv3_gradio_app"):
        os.chdir("ntv3_gradio_app")
        print("📁 Working directory: ntv3_gradio_app/\n")

    # Run checks
    results = []
    results.append(("Dependencies", check_dependencies()))
    results.append(("Data Files", check_data_files()))
    results.append(("Sequence Access", test_sequence_access()))

    # Model test is optional (downloads ~400MB if not cached)
    if all(r[1] for r in results):
        try:
            results.append(("Model Loading", test_model_loading()))
        except KeyboardInterrupt:
            print("\n⚠ Model test skipped (user interrupt)")

    # Summary
    print("\n" + "=" * 60)
    print("  Validation Summary")
    print("=" * 60)
    for name, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {name:20s} {status}")
    print("=" * 60 + "\n")

    if all(r[1] for r in results):
        print("🎉 All tests passed! Ready to run the app:")
        print("   python app.py")
        print("   or")
        print("   gradio app.py")
        return 0
    else:
        print(
            "⚠️  Some tests failed. Please fix the issues above before running the app."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
