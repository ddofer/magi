# MAGI — agent notes

Public repo (github.com/ddofer/magi). Layout, usage and citation are in README.md; this file holds only what is not there.

## Environment

- Python 3.13 in `.venv`, managed by `uv` (torch cu128, transformers 5). Run scripts with `uv run --no-sync python ...` or `.venv/bin/python`; add packages with `uv pip install` (no `pyproject.toml`; `pip` is not in the venv). One venv for everything, including `clinvar_pipeline/` and `rebuttal/`.
- NTv3 (`InstaDeepAI/NTv3_650M_post`) is gated: needs `HF_TOKEN`. Use the shared cache `HF_HOME=/storage/models/hf_home`; never copy models into the repo.
- LLM stages read `GEMINI_API_KEY` from the environment. Never hardcode keys in notebooks or scripts.

## SLURM

Heavy work goes through SLURM, **low priority** for this repo: `-p low` on every `srun`/`sbatch` (preemptible; see the cluster policy for the SIGTERM/requeue pattern on long runs).
At most 3 GPUs in total across this repo's jobs. No MLflow logging for this repo (user decision, overrides the cluster default).
The policy hook misreads `squeue -j <id>` as heavy work: use `squeue --me`. It also blocks even `--help` of GPU scripts outside SLURM: `srun -n 1 -p low -w sblxa123c -c 1 --mem=4G -t 0:05:00 ...`.
Network: login and compute nodes reach Hugging Face and PyPI; compute nodes also reach GitHub and Zenodo. NCBI, Ensembl and Google are blocked, so get annotations from HF mirrors (e.g. the vepyr VEP 116 cache).
**Always add `-w sblxa123c`**: this repo lives under `/opt/hpc`, which is not mounted on B300, so a job routed there dies before writing its log.
GPU (`--gres=gpu:1`, even for a quick check): `inference.py`, `inference_indel.py`, `inference_omia.py`, `smoke_test.py`, `gradio_app/app.py`, `paper_figures/regen_magi_app_panels.py`.
CPU only (still SLURM if over 4 cores, 16 GB or 10 min): `analysis.py`, `clinvar_pipeline/`, `paper_figures/regen_fig2_de_panels.py`, `rebuttal/analyses/`, notebooks in `analyses/`.
`inference.sbatch` is the template (srun inside, open-mode append, `HF_HUB_OFFLINE=1`: the cached model loads without a token); write others next to the script you run, logs under `logs/` (gitignored). Throwaway output goes in `scratch/` (gitignored).

Example:

```bash
srun -n 1 -p low -w sblxa123c --gres=gpu:1 -c 4 --mem=32G -t 0:20:00 uv run --no-sync python smoke_test.py --mode snp --model-size 100M
```

## Configuration

- `inference.py` settings are module constants near the top (`MODEL_NAME`, `INPUT_DATA_FILE`, `CONTEXT_LEN`, `BATCH_SIZE`, `USE_*`, `SHARE_WINDOWS`, `COMPILE`); its CLI flags only set those constants. `smoke_test.py` and `inference_indel.py` monkeypatch them. Run new variant sets with `sbatch inference.sbatch --input ... --output ...`.
- `paper_figures/` scripts find the repo from their own path; set `MAGI_ROOT` only when data lives elsewhere.
- GPU speed (A100, 64 kb, measured 2026-10-22). Tried and dropped for no gain: SDPA attention, compiling single blocks, max-autotune. bf16 needs the config `*_compute_dtype` fields (an outer autocast does nothing) and gave ~8% error on large BigWig deltas.
- The NTv3 stock forward runs all 24 species' BigWig heads; root `forward_trunk` skips them, `gradio_app/inference.py` patches the head to run only the batch's species.

## Private material (gitignored, present on disk)

`CONTEXT.md`, `docs/`, `rebuttal/`, `NAR_MAGI_Proof.PDF`, the delta parquets, `hg38.fa*`, `data/genomes/`, `data/results/` and the files listed at the top of `.gitignore`. Never `git add -f` any of them, and never paste their contents into tracked files. If `CONTEXT.md` exists, read it first: it is the glossary of canonical terms for prose and code. If `rebuttal/README.md` exists, follow its rules for the revision work.

## Gotchas

- `gradio_app/inference.py` and `gradio_app/analysis.py` are separate copies, not imports of the root files. A fix in one does not reach the other.
- `paper_figures/fig5/fig5b*`, `fig5c*`, `fig7/`, `fig8/` are finished vector files that cannot be rebuilt (upstream sources gone). Do not delete or regenerate them. See `paper_figures/RUNBOOK.md` before touching figures.
- `inference_omia.py` writes to a relative `data/`: the tracked `data/*_variants.parquet` files and `data/results/`. For a rerun, cd into a separate folder that has symlinks to `data/omia_raw_variants.csv` and `data/genomes`, and set `PYTHONPATH=<repo>`. `--validate-only` parses and checks REF letters without the model. Animal runs use a 16 kb window; human runs use 64 kb.
- `MLM_logprob_delta` is z-scored before the rules use it (`clinvar_pipeline/clinvar/variant_prep.py`). Label it as a z-score.
- Paths with spaces and odd nesting exist: `LLM pipeline inputs/LLM pipeline inputs/` and `analyses/analyses/`. Quote paths.

## For a future full rerun (decided 2026-10-22: not before)

- Track value scale. `D_BED_*`/`D_BW_*` are sigmoid of the head output: sigmoid(class 1) for BED instead of the softmax probability, and sigmoid(softplus) for BigWig, which squashes high-signal tracks. Kept for comparability with the paper. Switching to softmax (BED) and raw signal (BigWig) means re-deriving the rule cut-offs (±0.1/0.15/0.2, 0.08), rebuilding `gradio_app/data/magi_baseline_stats.csv` and `clinvar_pipeline/data/thresholds_*.csv`, and fixing the LLM prompt text, all in one versioned change. New columns already use the standard scales (`LFC_*` raw signal, `WIN_*` softmax).
- Window placement. NTv3 deltas move ~37% when the window shifts by 1 letter (the variant moves within its 128-letter token); every variant sits at the same place, so this is noise rather than bias. The user chose one pass per allele; treat it as a stated limitation.

## Tests and sub-project docs

`clinvar_pipeline/README.md` and `ARTIFACTS.md`, `gradio_app/README.md` and `SETUP_GUIDE.md`, `paper_figures/RUNBOOK.md`.
- Unit tests: `uv run --no-sync python -m unittest discover tests` (no model load; CPU, light).
