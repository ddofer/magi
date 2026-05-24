#!/usr/bin/env bash
# Sync inputs, gather canonical LLM, validate MANE, run analysis pipeline, validate figures.
set -euo pipefail
cd "$(dirname "$0")/.."
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate vep_env 2>/dev/null || true

echo "══ 1/6 Sync inputs → data/ ══"
python3 tools/sync_inputs.py --force

echo "══ 2/6 Canonical LLM (snp + indel parquets) ══"
python3 tools/gather_canonical_llm.py --force

echo "══ 3/6 Validate MANE metadata ══"
python3 tools/validate_mane.py --rebuild --write-report

echo "══ 4/6 Impact scores ══"
python3 scripts/06_compute_impact_scores.py

echo "══ 5/6 Figures ══"
python3 scripts/07_assign_mechanisms.py
python3 scripts/run_figures.py --fig all

echo "══ 6/6 Validate figure outputs ══"
python3 tools/validate_figures.py --check-outputs --write-report

echo "Done."
