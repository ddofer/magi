"""Regex canonicalizers for LLM mechanism labels (from fig3_a_b_c.ipynb)."""

from __future__ import annotations

import re

import pandas as pd


def build_primary_mechanism_canonicalizer(
    merge_mlm_like: bool = True,
    mlm_bucket_name: str = "MLM / Sequence-constraint signal",
    model_favors_bucket_name: str = "Model favors variant",
):
    CANON = {
        model_favors_bucket_name: [
            r"\bmodel favou?rs variant\b",
            r"\bmodel favou?rs the variant\b",
        ],
        "Benign / Neutral / None": [
            r"\bnone detected\b",
            r"\bnone identified\b",
            r"\bno mechanism detected\b",
            r"\bneutral\b",
            r"\bbenign\b",
            r"\btolerated\b",
            r"\bsilent\b",
            r"\bsynonymous\b",
        ],
        "Uncertain / Unknown": [
            r"\bunknown\b",
            r"\buncertain\b",
            r"\bundetermined\b",
            r"\binconclusive\b",
            r"\bunclear\b",
            r"\bna\b",
        ],
        "Pharmacogenomic": [r"\bpharmacogen", r"\bdrug response\b", r"\bpharmacogenetic\b"],
        "Splicing": [
            r"\bsplice\b",
            r"\bsplicing\b",
            r"\bdonor\b",
            r"\bacceptor\b",
            r"\bexon skipping\b",
            r"\bintron retention\b",
            r"\bcryptic splice\b",
            r"\bpseudoexon\b",
            r"\bexon\b",
            r"\bintron\b",
        ],
        "Start codon / Initiation / ORF": [
            r"\bstart codon\b",
            r"\bstart[- ]?loss\b",
            r"\binitiation\b",
            r"\btranslation initiation\b",
            r"\borf\b",
        ],
        "Nonsense / Truncation / NMD": [
            r"\bnonsense\b",
            r"\bstop[- ]?gain\b",
            r"\bpremature stop\b",
            r"\btruncat",
            r"\bnmd\b",
            r"\btermination codon\b",
        ],
        "Stop-loss / Protein extension": [r"\bstop[- ]?loss\b", r"\bprotein extension\b"],
        "Regulatory / Expression": [
            r"\bregulatory\b",
            r"\bexpression\b",
            r"\btranscription\b",
            r"\bcage\b",
            r"\b3'\s*utr\b",
            r"\b5'\s*utr\b",
        ],
        mlm_bucket_name: [
            r"\bmlm\b",
            r"\bllr\b",
            r"\blogprob\b",
            r"\bconstraint\b",
            r"\bconservation\b",
            r"\bmissense\b",
            r"\bprotein\b",
            r"\bcoding\b",
            r"\bnon[- ]?synonymous\b",
        ],
    }
    PRIORITY = [
        model_favors_bucket_name,
        "Benign / Neutral / None",
        "Uncertain / Unknown",
        "Pharmacogenomic",
        "Splicing",
        "Start codon / Initiation / ORF",
        "Nonsense / Truncation / NMD",
        "Stop-loss / Protein extension",
        "Regulatory / Expression",
        mlm_bucket_name,
    ]
    CANON_RE = {k: [re.compile(p, re.I) for p in pats] for k, pats in CANON.items()}
    EFFECT_GUARD_RE = re.compile(
        r"\b(model favou?rs variant|missense|protein|coding|mlm|llr|"
        r"constraint|conservation|splice|nmd|nonsense|truncat)\b",
        re.I,
    )

    def _canonicalize(label):
        if label is None or (isinstance(label, float) and pd.isna(label)):
            return "Uncertain / Unknown"
        s = str(label).strip()
        if not s:
            return "Uncertain / Unknown"
        s_norm = s.replace("_", " ").replace("|", " ").replace(",", " ").lower()
        benign_hit = any(rx.search(s_norm) for rx in CANON_RE["Benign / Neutral / None"])
        if benign_hit and EFFECT_GUARD_RE.search(s_norm):
            for bucket in [b for b in PRIORITY if b != "Benign / Neutral / None"]:
                for rx in CANON_RE[bucket]:
                    if rx.search(s_norm):
                        return bucket
            return "Other / Unclassified"
        for bucket in PRIORITY:
            for rx in CANON_RE[bucket]:
                if rx.search(s_norm):
                    return bucket
        return "Other / Unclassified"

    return _canonicalize


def build_rationale_mechanism_canonicalizer():
    CANON = {
        "Splicing": [r"\bsplice\b", r"\bsplicing\b", r"\bdonor\b", r"\bacceptor\b"],
        "Start codon / Initiation / ORF": [r"\bstart codon\b", r"\borf\b", r"\binitiation\b"],
        "Nonsense / Truncation / NMD": [r"\bnonsense\b", r"\btruncat", r"\bnmd\b"],
        "Frameshift": [r"\bframeshift\b"],
        "Missense / Protein function": [r"\bmissense\b", r"\bprotein function\b"],
        "Sequence constraint / Conservation": [r"\bconstraint\b", r"\bconservation\b", r"\bconserved\b"],
        "Benign / Neutral / None": [r"\bbenign\b", r"\bneutral\b", r"\bsynonymous\b"],
        "Uncertain / Unknown": [r"\bunknown\b", r"\buncertain\b"],
    }
    PRIORITY = list(CANON.keys())
    CANON_RE = {k: [re.compile(p, re.I) for p in pats] for k, pats in CANON.items()}

    def _canonicalize(label):
        if label is None or (isinstance(label, float) and pd.isna(label)):
            return "Uncertain / Unknown"
        s_norm = str(label).replace("_", " ").lower()
        for bucket in PRIORITY:
            for rx in CANON_RE[bucket]:
                if rx.search(s_norm):
                    return bucket
        return "Other / Unclassified"

    return _canonicalize
