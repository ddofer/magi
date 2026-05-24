#!/usr/bin/env bash
set -euo pipefail
cd /home/stavz/masters/Clinvar/clinvar_pipeline
source /home/stavz/miniforge3/etc/profile.d/conda.sh
conda activate vep_env
LOG=/home/stavz/masters/Clinvar/clinvar_pipeline/output/snp_pipeline_run.log
exec > >(tee "$LOG") 2>&1
echo "=== $(date -Is) SNP pipeline start ==="
python3 scripts/01_prepare_variants.py --snp
python3 scripts/02_extract_signals.py --variant-type snp
python3 scripts/03_build_prompts.py --variant-type snp
python3 tools/compare_outputs.py --prefix all
echo "=== $(date -Is) SNP pipeline done ==="
