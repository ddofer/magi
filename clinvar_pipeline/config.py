"""Paths and settings for the refactored ClinVar pipeline.

``data/`` — synced inputs and gathered canonical artifacts (LLM results).
``output/`` — everything produced by running scripts (parquet stages, impact,
mechanisms, figures, validation reports). Kept separate so ``data/`` can be
refreshed from the parent repo without overwriting pipeline results.
"""

from __future__ import annotations

from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_ROOT.parent
OUTPUT_ROOT = PACKAGE_ROOT / "output"
DATA_ROOT = PACKAGE_ROOT / "data"


def _in_repo(rel: str) -> str:
    return str(PROJECT_ROOT / rel)


def _out(rel: str) -> str:
    return str(OUTPUT_ROOT / rel)


def _data(rel: str) -> str:
    return str(DATA_ROOT / rel)


def _data_or_repo(local_rel: str, repo_rel: str) -> str:
    local = DATA_ROOT / local_rel
    if local.exists():
        return str(local)
    return _in_repo(repo_rel)


PATHS = {
    # ── inputs (prefer data/ copies) ──
    "snp_deltas": _data_or_repo("parquet/clinvar_new_deltas.parquet", "parquet/clinvar_new_deltas.parquet"),
    "indel_deltas": _data_or_repo("parquet/clinvar_indel_deltas.parquet", "parquet/clinvar_indel_deltas.parquet"),
    "variant_summary": _data_or_repo("clinvar_tables/variant_summary.txt.gz", "clinvar_tables/variant_summary.txt.gz"),
    "submission_summary": _data_or_repo(
        "clinvar_tables/submission_summary.txt.gz", "clinvar_tables/submission_summary.txt.gz"
    ),
    "MANE_gff": _data_or_repo(
        "metadata/MANE.GRCh38.v1.5.refseq_genomic.gff.gz", "metadata/MANE.GRCh38.v1.5.refseq_genomic.gff.gz"
    ),
    "MANE_processed": _data_or_repo("metadata/MANE_processed.csv", "metadata/MANE_processed.csv"),
    "Promoter_processed": _data_or_repo("metadata/Promoter_processed.csv", "metadata/Promoter_processed.csv"),
    "tracks_metadata": _data_or_repo(
        "metadata/functional_tracks_metadata_human.csv", "metadata/functional_tracks_metadata_human.csv"
    ),
    "snp_z_thresholds": _data_or_repo("thresholds_snps.csv", "thresholds_snps.csv"),
    "indel_z_thresholds": _data_or_repo("thresholds_indels.csv", "thresholds_indels.csv"),
    "uncertain_ids": _data_or_repo("uncertain_ids.csv", "uncertain_ids.csv"),
    "af_snps": _data_or_repo("AF/snps_AF_sub.csv", "AF/snps_AF_sub.csv"),
    "af_indels": _data_or_repo("AF/indel_AF_sub.csv", "AF/indel_AF_sub.csv"),
    # canonical LLM judge results (one file per variant type — built by tools/gather_canonical_llm.py)
    "snp_llm_results": _data("llm/snp_llm_results.parquet"),
    "indel_llm_results": _data("llm/indel_llm_results.parquet"),
    # parent-repo fallbacks when gathered parquets are missing
    "snp_eval_csv": _in_repo("parquet/snp_evaluation_results_v2.csv"),
    "indel_eval_csv": _in_repo("parquet/indel_evaluation_results_v2.csv"),
    # legacy reference outputs (read-only comparison)
    "legacy_snp_strict": _in_repo("parquet/snps_strict_test.parquet"),
    "legacy_indel_strict": _in_repo("parquet/indels_strict_test.parquet"),
    "legacy_snp_annotated": _in_repo("parquet/snps_annotated_test.parquet"),
    "legacy_indel_annotated": _in_repo("parquet/indels_annotated_test.parquet"),
    "legacy_snp_signaled": _in_repo("parquet/annotated_snps_signaled.parquet"),
    "legacy_indel_signaled": _in_repo("parquet/annotated_indels_signaled.parquet"),
    # ── pipeline outputs ──
    "snp_strict": _out("parquet/snps_strict.parquet"),
    "indel_strict": _out("parquet/indels_strict.parquet"),
    "snp_annotated": _out("parquet/snps_annotated.parquet"),
    "indel_annotated": _out("parquet/indels_annotated.parquet"),
    "snp_annotated_signaled": _out("parquet/snps_signaled.parquet"),
    "indel_annotated_signaled": _out("parquet/indels_signaled.parquet"),
    "snp_with_prompts": _out("parquet/snps_with_prompts.parquet"),
    "indel_with_prompts": _out("parquet/indels_with_prompts.parquet"),
    "snp_deltas_impact": _out("impact/snp_deltas_impact.parquet"),
    "indel_deltas_impact": _out("impact/indel_deltas_impact.parquet"),
    # derived mechanism tables (computed from signaled parquets, not copied from parent)
    "snp_mechanisms_csv": _out("mechanisms/snps_assigned_mechanisms.csv"),
    "indel_mechanisms_csv": _out("mechanisms/indels_assigned_mechanisms.csv"),
    # figure outputs (flat names under output/figures/)
    "figures_root": _out("figures"),
    "fig1c": _out("figures/fig1c.png"),
    "fig2a": _out("figures/fig2a.png"),
    "fig2b": _out("figures/fig2b.png"),
    "fig3a": _out("figures/fig3a.png"),
    "fig3b": _out("figures/fig3b.png"),
    "fig3c": _out("figures/fig3c.png"),
    "fig3d": _out("figures/fig3d.png"),
    "fig3e": _out("figures/fig3e.png"),
    "fig3e_indel": _out("figures/fig3e_indel.png"),
    "figs1": _out("figures/figs1.png"),
    # OMIA animal variants (notebook: notebooks/analysis_v2 animals.ipynb)
    "animals_snp_deltas": _data("parquet/deltas_animals_snp.parquet"),
    "animals_indel_deltas": _data("parquet/deltas_animals_indel.parquet"),
    "animals_snp_signaled": _out("parquet/animals_snp_signaled.parquet"),
    "animals_indel_signaled": _out("parquet/animals_indel_signaled.parquet"),
    "animals_snp_with_prompts": _out("parquet/animals_snp_with_prompts.parquet"),
    "animals_indel_with_prompts": _out("parquet/animals_indel_with_prompts.parquet"),
    "animals_snp_eval": _data("llm/animals_snp_evaluation_results.parquet"),
    "animals_indel_eval": _data("llm/animals_indel_evaluation_results_2.parquet"),
    "omia_all_species_inferred_full": _data_or_repo(
        "parquet/omia_all_species_inferred_full.parquet",
        "data/results/omia_all_species_inferred_full.parquet",
    ),
    "omia_all_species_inferred_relevant": _data_or_repo(
        "parquet/omia_all_species_inferred_relevant.parquet",
        "data/results/omia_all_species_inferred_relevant.parquet",
    ),
    "validation_root": _out("validation"),
    "cohort_funnel": _out("validation/cohort_funnel.json"),
    "batch_work_dir": _out("batch_eval"),
}

IMPACT_MERGE_KEYS = ["chrom", "pos", "ref", "alt"]
IMPACT_METADATA_COLS = ["label", "variant_type", "#VariationID", "variant_id"]
IMPACT_SCORE_COLS = [
    "gnomADe_AF",
    "BED_mean_abs",
    "BED_max_abs",
    "BED_sum_log",
    "BED_top3_mean",
    "BW_z_mean_abs",
    "BW_z_max_abs",
    "BW_z_sum_log",
    "BW_z_top10_mean",
    "Global_z_mean_abs",
    "Global_z_max_abs",
    "Global_z_sum_log",
    "Composite_mean",
    "Composite_top",
]

MIN_GOLD_STARS = 2
DEFAULT_TOP_K_SIGNALS = 5
DEFAULT_ANIMALS_TOP_K_SIGNALS = 10
DEFAULT_LLM_MODEL = "gemini-3-flash-preview"
