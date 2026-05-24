"""Canonical LLM result columns used by figure scripts."""

from __future__ import annotations

# Core judge outputs (fig3abc + fig3e)
SNP_LLM_CORE = [
    "#VariationID",
    "concordance",
    "signal_category",
    "MLM_category",
    "primary_signal_mechanism",
    "rationale_mechanism",
    "explanation",
    "key_signals",
    "nt_missed",
    "notes",
    "label",
]

INDEL_LLM_CORE = [
    "#VariationID",
    "concordance",
    "signal_category",
    "MLM_category",
    "primary_signal_mechanism",
    "rationale_mechanism",
    "explanation",
    "key_signals",
    "nt_missed",
    "notes",
    "embedding_impact",
    "variant_type",
    "indel_size",
]
