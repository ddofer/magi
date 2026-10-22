# MAGI — agent notes

Public repo (github.com/ddofer/magi). Layout, usage and citation are in README.md; this file holds only what is not there.

## Environment

- Python 3.13 in `.venv`, managed by `uv` (torch cu128, transformers 5). Run scripts with `uv run --no-sync python ...` or `.venv/bin/python`; add packages with `uv pip install` (no `pyproject.toml`; `pip` is not in the venv). One venv for everything, including `clinvar_pipeline/` and `rebuttal/`.
- NTv3 (`InstaDeepAI/NTv3_650M_post`) is gated: needs `HF_TOKEN`. Use the shared cache `HF_HOME=/storage/models/hf_home`; never copy models into the repo.
- LLM stages read `GEMINI_API_KEY` from the environment. Never hardcode keys in notebooks or scripts.

## SLURM

Heavy work goes through SLURM, **low priority** for this repo: `-p low` on every `srun`/`sbatch` (preemptible; see the cluster policy for the SIGTERM/requeue pattern on long runs).
GPU (`--gres=gpu:1`, even for a quick check): `inference.py`, `inference_indel.py`, `inference_omia.py`, `smoke_test.py`, `gradio_app/app.py`, `paper_figures/regen_magi_app_panels.py`.
CPU only (still SLURM if over 4 cores, 16 GB or 10 min): `analysis.py`, `clinvar_pipeline/`, `paper_figures/regen_fig2_de_panels.py`, `rebuttal/analyses/`, notebooks in `analyses/`.
No `.sbatch` scripts exist yet; write one next to the script you run, logs under `logs/` (gitignored). Throwaway output goes in `scratch/` (gitignored).

Example:

```bash
srun -n 1 -p low --gres=gpu:1 -c 4 --mem=32G -t 0:20:00 uv run --no-sync python smoke_test.py --mode snp --model-size 100M
```

## Configuration

- `inference.py` and `inference_indel.py` have no CLI flags. Settings are module constants near the top (`MODEL_NAME`, `INPUT_DATA_FILE`, `CONTEXT_LEN`, `BATCH_SIZE`, `USE_*`). `smoke_test.py` overrides them by monkeypatching. Only `inference_omia.py` has argparse.
- `paper_figures/` scripts find the repo from their own path; set `MAGI_ROOT` only when data lives elsewhere.

## Private material (gitignored, present on disk)

`CONTEXT.md`, `docs/`, `rebuttal/`, `NAR_MAGI_Proof.PDF`, the delta parquets, `hg38.fa*`, `data/genomes/`, `data/results/` and the files listed at the top of `.gitignore`. Never `git add -f` any of them, and never paste their contents into tracked files. If `CONTEXT.md` exists, read it first: it is the glossary of canonical terms for prose and code. If `rebuttal/README.md` exists, follow its rules for the revision work.

## Gotchas

- `gradio_app/inference.py` and `gradio_app/analysis.py` are separate copies, not imports of the root files. A fix in one does not reach the other.
- `paper_figures/fig5/fig5b*`, `fig5c*`, `fig7/`, `fig8/` are finished vector files that cannot be rebuilt (upstream sources gone). Do not delete or regenerate them. See `paper_figures/RUNBOOK.md` before touching figures.
- Paths with spaces and odd nesting exist: `LLM pipeline inputs/LLM pipeline inputs/` and `analyses/analyses/`. Quote paths.

## Tests and sub-project docs

`clinvar_pipeline/README.md` and `ARTIFACTS.md`, `gradio_app/README.md` and `SETUP_GUIDE.md`, `paper_figures/RUNBOOK.md`.
- Unit tests: `uv run --no-sync python -m unittest discover tests` (no model load; CPU, light).
