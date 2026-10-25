#!/bin/bash
# NTv3 variant inference on one GPU, preemptible. Extra arguments go to inference.py, e.g.
#   sbatch inference.sbatch --input my_variants.parquet --output my_deltas.parquet
#   sbatch inference.sbatch --region chr17:43044295-43045000 --output brca1_snvs.parquet
#   sbatch inference.sbatch --input clinvar_indel.parquet --output clinvar_indel_deltas.parquet --embeddings
# A preempted run is requeued and resumes from <output>.parts/ (finished chunks are kept).
# Submit from the repo root; logs/ must exist (mkdir -p logs) or SLURM cannot write the log.
#SBATCH --job-name=magi-inference
#SBATCH -p low
#SBATCH -w sblxa123c
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --open-mode=append
#SBATCH --output=logs/%x-%j.out
cd "$SLURM_SUBMIT_DIR"  # submit from the repo root
export HF_HOME=/storage/models/hf_home HF_HUB_OFFLINE=1
echo "start $(date)  restart=${SLURM_RESTART_COUNT:-0}"
srun .venv/bin/python inference.py "$@"
