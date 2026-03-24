"""
describe_assignment_rules.py
=============================
Programmatically describes the EXACT rules used by NT mechanism assignment.
Reads from the hypothesis logic dicts and prints human-readable rule descriptions.

Run standalone or import describe_all_rules() / describe_mechanism().
"""


# =============================================================================
# 1. HYPOTHESIS LOGIC — FLAT THRESHOLDS, NO STRENGTH TIERS
# =============================================================================
# 
# Design principles:
#   - Single low threshold per feature (captures more variants)
#   - No WEAK_* categories — every mechanism is either matched or not
#   - Priority still determines which mechanism "wins" as primary
#   - OR logic for primary conditions, AND logic for secondary
#   - Secondary conditions use even lower thresholds (supportive evidence)
#   - Tiers 0–4: shared across SNP and INDEL (D_BED_* features)
#   - Tier 5: variant-type-specific MLM / embedding signals

# =============================================================================
# 1a. SHARED FEATURE-LEVEL LOGIC (Tiers 0–4)
# =============================================================================

_SHARED_LOGIC = {
    # ==========================================================================
    # TIER 0: DEFINITIVE LOSS-OF-FUNCTION (Priority 0)
    # ==========================================================================
    
    "SPLICE_SITE_DESTROYED": {
        "primary": [
            {"feature": "D_BED_splice_donor", "threshold": -0.2, "direction": "negative"},
            {"feature": "D_BED_splice_acceptor", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Loss of canonical splice donor or acceptor site",
        "priority": 0
    },
    
    "START_LOSS": {
        "primary": [
            {"feature": "D_BED_start_codon", "threshold": -0.15, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Loss of translation initiation codon",
        "priority": 0
    },
    
    "STOP_CODON_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_stop_codon", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Disruption of stop codon causing read-through",
        "priority": 0
    },

    # ==========================================================================
    # TIER 1: SPLICE & EXON MECHANISMS (Priority 1)
    # ==========================================================================
    
    "CRYPTIC_SPLICE_GAIN": {
        "primary": [
            {"feature": "D_BED_splice_donor", "threshold": 0.2, "direction": "positive"},
            {"feature": "D_BED_splice_acceptor", "threshold": 0.2, "direction": "positive"}
        ],
        "secondary": [],
        "description": "Creation of new cryptic splice site",
        "priority": 1
    },
    
    "SPLICE_INDUCED_EXONIZATION": {
        "primary": [
            {"feature": "D_BED_exon", "threshold": 0.15, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_intron", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "Intronic sequence included as coding exon via splice disruption",
        "priority": 1
    },
    
    "SPLICE_INDUCED_INTRON_RETENTION": {
        "primary": [
            {"feature": "D_BED_intron", "threshold": 0.2, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_splice_donor", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "Intron retention due to splice site weakening",
        "priority": 1
    },
    
    "EXON_SKIPPING_INDUCED": {
        "primary": [
            {"feature": "D_BED_skipped_exon", "threshold": 0.15, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_always_on_exon", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "Induction of exon skipping",
        "priority": 1
    },
    
    "CONSTITUTIVE_EXON_LOSS": {
        "primary": [
            {"feature": "D_BED_always_on_exon", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Loss of constitutively included exon",
        "priority": 1
    },
    
    "EXON_LOSS": {
        "primary": [
            {"feature": "D_BED_exon", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Loss of exon identity",
        "priority": 1
    },
    
    "CODING_LOSS": {
        "primary": [
            {"feature": "D_BED_exon", "threshold": -0.15, "direction": "negative"}
        ],
        "secondary": [
            {"feature": "D_BED_ORF", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "Loss of coding exon identity with ORF disruption",
        "priority": 1
    },
    
    "ORF_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_ORF", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Disruption of open reading frame",
        "priority": 1
    },
    
    "CRYPTIC_EXON_INCLUSION": {
        "primary": [
            {"feature": "D_BED_skipped_exon", "threshold": 0.2, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_exon", "threshold": 0.1, "direction": "positive"}
        ],
        "description": "Inclusion of cryptic/normally-skipped exon into transcript",
        "priority": 1
    },

    # ==========================================================================
    # TIER 2: REGULATORY & UTR (Priority 2)
    # ==========================================================================
    
    "START_LOSS_WITH_UTR_EXTENSION": {
        "primary": [
            {"feature": "D_BED_start_codon", "threshold": -0.1, "direction": "negative"}
        ],
        "secondary": [
            {"feature": "D_BED_5UTR+", "threshold": 0.1, "direction": "positive"}
        ],
        "description": "Start codon loss with consequent 5'UTR extension",
        "priority": 2
    },
    
    "START_LOSS_WITH_ORF_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_start_codon", "threshold": -0.1, "direction": "negative"}
        ],
        "secondary": [
            {"feature": "D_BED_ORF", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "Start codon loss with downstream ORF disruption",
        "priority": 2
    },
    
    "PROMOTER_DESTRUCTION": {
        "primary": [
            {"feature": "D_BED_promoter_Tissue_specific", "threshold": -0.2, "direction": "negative"},
            {"feature": "D_BED_promoter_Tissue_invariant", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Destruction of promoter element",
        "priority": 2
    },
    
    "ENHANCER_DESTRUCTION": {
        "primary": [
            {"feature": "D_BED_enhancer_Tissue_specific", "threshold": -0.2, "direction": "negative"},
            {"feature": "D_BED_enhancer_Tissue_invariant", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Destruction of enhancer element",
        "priority": 2
    },
    
    "REGULATORY_GAIN": {
        "primary": [
            {"feature": "D_BED_enhancer_Tissue_specific", "threshold": 0.2, "direction": "positive"},
            {"feature": "D_BED_promoter_Tissue_specific", "threshold": 0.2, "direction": "positive"},
            {"feature": "D_BED_enhancer_Tissue_invariant", "threshold": 0.2, "direction": "positive"},
            {"feature": "D_BED_promoter_Tissue_invariant", "threshold": 0.2, "direction": "positive"}
        ],
        "secondary": [],
        "description": "Gain of regulatory element (potential ectopic expression)",
        "priority": 2
    },
    
    "UTR5_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_5UTR+", "threshold": -0.15, "direction": "negative"},
            {"feature": "D_BED_5UTR-", "threshold": 0.15, "direction": "positive"}
        ],
        "secondary": [],
        "description": "Disruption of 5' UTR structure or regulation",
        "priority": 2
    },
    
    "UTR3_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_3UTR+", "threshold": -0.15, "direction": "negative"},
            {"feature": "D_BED_3UTR-", "threshold": 0.15, "direction": "positive"}
        ],
        "secondary": [],
        "description": "Disruption of 3' UTR structure or regulation",
        "priority": 2
    },
    
    "POLYA_SIGNAL_LOSS": {
        "primary": [
            {"feature": "D_BED_polyA_signal", "threshold": -0.15, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Loss of polyadenylation signal",
        "priority": 2
    },
    
    "CTCF_BOUNDARY_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_CTCF-bound", "threshold": -0.15, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Disruption of CTCF binding / chromatin boundary element",
        "priority": 2
    },

    # ==========================================================================
    # TIER 3: COMPLEX / RECIPROCAL CHANGES (Priority 3)
    # ==========================================================================
    
    "CRYPTIC_ORF_GAIN": {
        "primary": [
            {"feature": "D_BED_ORF", "threshold": 0.2, "direction": "positive"}
        ],
        "secondary": [],
        "description": "Creation of cryptic open reading frame",
        "priority": 3
    },
    
    "CRYPTIC_ORF_GAIN_WITH_INTRON_LOSS": {
        "primary": [
            {"feature": "D_BED_ORF", "threshold": 0.15, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_intron", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "ORF gain with intron loss — likely splice-induced exonization",
        "priority": 3
    },
    
    "ABERRANT_EXON_INCLUSION": {
        "primary": [
            {"feature": "D_BED_always_on_exon", "threshold": 0.2, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_intron", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "Aberrant inclusion of intronic sequence as constitutive exon",
        "priority": 3
    },
    
    "INTRON_RETENTION": {
        "primary": [
            {"feature": "D_BED_intron", "threshold": 0.2, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_exon", "threshold": -0.1, "direction": "negative"}
        ],
        "description": "Retention of intron in mature transcript",
        "priority": 3
    },
    
    "RECIPROCAL_EXON_GAIN_INTRON_LOSS": {
        "primary": [
            {"feature": "D_BED_exon", "threshold": 0.15, "direction": "positive"}
        ],
        "secondary": [
            {"feature": "D_BED_intron", "threshold": -0.15, "direction": "negative"}
        ],
        "description": "Exon gain with reciprocal intron loss (exonization)",
        "priority": 3
    },
    
    "RECIPROCAL_EXON_LOSS_INTRON_GAIN": {
        "primary": [
            {"feature": "D_BED_exon", "threshold": -0.15, "direction": "negative"}
        ],
        "secondary": [
            {"feature": "D_BED_intron", "threshold": 0.15, "direction": "positive"}
        ],
        "description": "Exon loss with reciprocal intron gain (intronization)",
        "priority": 3
    },
    
    "GENE_ARCHITECTURE_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_protein_coding_gene", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Disruption of overall gene architecture",
        "priority": 3
    },
    
    "LNCRNA_DISRUPTION": {
        "primary": [
            {"feature": "D_BED_lncRNA", "threshold": -0.2, "direction": "negative"}
        ],
        "secondary": [],
        "description": "Disruption of lncRNA",
        "priority": 3
    },

    # ==========================================================================
    # TIER 4: CATCH-ALL — any D_BED signal above noise floor (Priority 4)
    # ==========================================================================
    
    "SPLICE_CHANGE": {
        "primary": [
            {"feature": "D_BED_splice_donor", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_splice_acceptor", "threshold": 0.08, "direction": "any"}
        ],
        "secondary": [],
        "description": "Detectable change in splice site signal",
        "priority": 4
    },
    
    "REGULATORY_CHANGE": {
        "primary": [
            {"feature": "D_BED_enhancer_Tissue_specific", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_promoter_Tissue_specific", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_enhancer_Tissue_invariant", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_promoter_Tissue_invariant", "threshold": 0.08, "direction": "any"}
        ],
        "secondary": [],
        "description": "Detectable change in regulatory element",
        "priority": 4
    },
    
    "EXON_CHANGE": {
        "primary": [
            {"feature": "D_BED_exon", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_always_on_exon", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_skipped_exon", "threshold": 0.08, "direction": "any"}
        ],
        "secondary": [],
        "description": "Detectable change in exon identity",
        "priority": 4
    },
    
    "ORF_CHANGE": {
        "primary": [
            {"feature": "D_BED_ORF", "threshold": 0.08, "direction": "any"}
        ],
        "secondary": [],
        "description": "Detectable change in ORF",
        "priority": 4
    },
    
    "INTRON_CHANGE": {
        "primary": [
            {"feature": "D_BED_intron", "threshold": 0.08, "direction": "any"}
        ],
        "secondary": [],
        "description": "Detectable change in intron signal",
        "priority": 4
    },
    
    "UTR_CHANGE": {
        "primary": [
            {"feature": "D_BED_5UTR+", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_5UTR-", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_3UTR+", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_3UTR-", "threshold": 0.08, "direction": "any"}
        ],
        "secondary": [],
        "description": "Detectable change in UTR signal",
        "priority": 4
    },
    
    "GENE_CHANGE": {
        "primary": [
            {"feature": "D_BED_protein_coding_gene", "threshold": 0.08, "direction": "any"},
            {"feature": "D_BED_lncRNA", "threshold": 0.08, "direction": "any"}
        ],
        "secondary": [],
        "description": "Detectable change in gene-level signal",
        "priority": 4
    },
}


# =============================================================================
# 1b. VARIANT-TYPE-SPECIFIC TIER 5 (MLM / Embedding fallbacks)
# =============================================================================
#
# These fire only when no D_BED feature-level annotation matched (Tiers 0–4).
#
# SNP signals:
#   LLR = log(P(ALT)/P(REF))          — primary constraint signal
#   MLM_Prior = P(REF)                 — model confidence in reference allele
#
# INDEL signals:
#   MLM_logprob_delta = logP(ALT) - logP(REF)  — analogous to SNP LLR
#   MLM_KL_max = max KL divergence across window — distributional disruption
#   EMB_cosine_dist = cosine distance between ref/alt embeddings — representation disruption

_TIER5_SNP = {
    # Primary: LLR ≤ -1.0 (ALT ~2.7x less likely than REF)
    # Secondary: MLM_Prior ≥ 0.3 (model was reasonably confident about REF)
    # Together: "the model was confident about REF and strongly prefers it over ALT"
    "SEQUENCE_CONSTRAINED": {
        "primary": [
            {"feature": "LLR", "threshold": -1.0, "direction": "negative"}
        ],
        "secondary": [
            {"feature": "MLM_Prior", "threshold": 0.3, "direction": "positive"}
        ],
        "description": "Position under sequence constraint — MLM confident in REF and strongly prefers it over ALT",
        "priority": 5
    },
}

_TIER5_INDEL = {
    # Primary (OR): logprob delta or KL divergence catches disruption
    #   MLM_logprob_delta ≤ -1.0 : direct analog of SNP LLR
    #   MLM_KL_max ≥ 2.0 : distributional disruption even if per-position logprobs look OK
    "SEQUENCE_CONSTRAINED": {
        "primary": [
            {"feature": "MLM_logprob_delta", "threshold": -1.0, "direction": "negative"},
            {"feature": "MLM_KL_max", "threshold": 2.0, "direction": "positive"}
        ],
        "secondary": [],
        "description": "Sequence under constraint — MLM probability or distributional disruption detected",
        "priority": 5
    },
    
    # Separate mechanism: embedding-space disruption
    # Catches cases where probability-space metrics look fine but the internal
    # representation is substantially altered — a different modality of signal.
    "EMBEDDING_DISRUPTED": {
        "primary": [
            {"feature": "EMB_cosine_dist", "threshold": 0.15, "direction": "positive"}
        ],
        "secondary": [],
        "description": "Embedding-space disruption — MLM internal representation substantially altered by variant",
        "priority": 5
    },
}


def get_hypothesis_logic(variant_type="snp"):
    """
    Return the full hypothesis logic dict for a given variant type.
    
    Parameters
    ----------
    variant_type : str
        "snp" or "indel" (case-insensitive).
    
    Returns
    -------
    dict : Complete hypothesis logic (Tiers 0–5).
    """
    vt = variant_type.strip().lower()
    logic = dict(_SHARED_LOGIC)  # copy Tiers 0–4
    
    if vt == "snp":
        logic.update(_TIER5_SNP)
    elif vt == "indel":
        logic.update(_TIER5_INDEL)
    else:
        raise ValueError(f"Unknown variant_type '{variant_type}'. Use 'snp' or 'indel'.")
    
    return logic


def _direction_to_operator(direction: str) -> str:
    return {
        "negative": "≤",
        "positive": "≥",
        "any": "|abs| ≥"
    }[direction]


def _format_condition(cond: dict) -> str:
    feat = cond['feature'].replace('D_BED_', '')
    op = _direction_to_operator(cond['direction'])
    thresh = cond['threshold']
    return f"{feat} {op} {thresh}"


def describe_mechanism(name: str, rules: dict) -> str:
    lines = []
    lines.append(f"MECHANISM: {name}")
    lines.append(f"  Priority tier: {rules['priority']}")
    lines.append(f"  Description: {rules['description']}")
    
    primary = rules['primary']
    if len(primary) == 1:
        lines.append(f"  Rule: {_format_condition(primary[0])}")
    else:
        cond_strs = [_format_condition(c) for c in primary]
        lines.append(f"  Rule (ANY of): {' OR '.join(cond_strs)}")
    
    secondary = rules.get('secondary', [])
    if secondary:
        cond_strs = [_format_condition(c) for c in secondary]
        lines.append(f"  AND also requires (ALL of): {' AND '.join(cond_strs)}")
    
    return "\n".join(lines)


def describe_all_rules(variant_type="snp") -> str:
    """
    Generate a complete description of the assignment logic for a variant type.
    
    Parameters
    ----------
    variant_type : str
        "snp" or "indel".
    """
    vt = variant_type.strip().lower()
    logic = get_hypothesis_logic(vt)
    
    sections = []
    
    sections.append("=" * 72)
    sections.append(f"NT MECHANISM ASSIGNMENT RULES — {vt.upper()} — COMPLETE SPECIFICATION")
    sections.append("=" * 72)
    sections.append("")
    
    sections.append("ASSIGNMENT ALGORITHM:")
    sections.append("  1. Mechanisms are evaluated in priority order (0 = highest).")
    sections.append("  2. Within each mechanism:")
    sections.append("     a. PRIMARY conditions use OR logic: at least one must pass.")
    sections.append("     b. SECONDARY conditions use AND logic: all must pass (if any exist).")
    sections.append("  3. A variant can match MULTIPLE mechanisms.")
    sections.append("     - NT_Mechanism = the first (highest priority) match.")
    sections.append("     - NT_All_Mechanisms = all matches, listed with their priority.")
    sections.append("  4. If no mechanism matches, the variant is UNCERTAIN_SIGNIFICANCE.")
    sections.append("  5. NaN / missing feature values are treated as 0.")
    sections.append("")
    sections.append("DIRECTION SEMANTICS:")
    sections.append("  'negative'  → value ≤ threshold  (signal LOSS)")
    sections.append("  'positive'  → value ≥ threshold  (signal GAIN)")
    sections.append("  'any'       → |value| ≥ |threshold|  (any change)")
    sections.append("")
    sections.append("FEATURE SCALES:")
    sections.append("  D_BED_*          delta probabilities in [-1, 1].  D = P(variant) - P(reference).")
    if vt == "snp":
        sections.append("  LLR              log(P(ALT) / P(REF)).  Negative = REF preferred = constrained.")
        sections.append("  MLM_Prior        P(REF) at the position.  Higher = model more confident about REF.")
    elif vt == "indel":
        sections.append("  MLM_logprob_delta  logP(ALT) - logP(REF) across window.  Negative = constrained.")
        sections.append("  MLM_KL_max       Max KL divergence across affected positions.  Higher = more disrupted.")
        sections.append("  EMB_cosine_dist  Cosine distance between ref/alt embeddings.  Higher = more disrupted.")
    sections.append("")
    
    # Group by priority tier
    by_tier = {}
    for name, rules in logic.items():
        tier = rules['priority']
        by_tier.setdefault(tier, []).append((name, rules))
    
    tier_labels = {
        0: "TIER 0 — DEFINITIVE LOSS-OF-FUNCTION",
        1: "TIER 1 — SPLICE & EXON MECHANISMS",
        2: "TIER 2 — REGULATORY & UTR",
        3: "TIER 3 — COMPLEX / RECIPROCAL",
        4: "TIER 4 — CATCH-ALL (any D_BED signal above noise floor)",
        5: f"TIER 5 — MLM / SEQUENCE CONSTRAINT ({vt.upper()}-specific)",
    }
    
    for tier in sorted(by_tier.keys()):
        sections.append("-" * 72)
        sections.append(tier_labels.get(tier, f"TIER {tier}"))
        sections.append("-" * 72)
        for name, rules in by_tier[tier]:
            sections.append("")
            sections.append(describe_mechanism(name, rules))
        sections.append("")
    
    # Summary table
    sections.append("=" * 72)
    sections.append("SUMMARY TABLE")
    sections.append("=" * 72)
    sections.append(f"{'Mechanism':<42} {'Tier':>4}  {'#Primary':>8}  {'#Secondary':>10}  Primary Features")
    sections.append("-" * 120)
    
    for tier in sorted(by_tier.keys()):
        for name, rules in by_tier[tier]:
            n_p = len(rules['primary'])
            n_s = len(rules.get('secondary', []))
            feats = ", ".join(set(c['feature'].replace('D_BED_', '') for c in rules['primary']))
            sections.append(f"{name:<42} {tier:>4}  {n_p:>8}  {n_s:>10}  {feats}")
    
    total = len(logic)
    sections.append(f"\nTotal mechanisms defined: {total}")
    
    return "\n".join(sections)


def describe_feature_usage(variant_type="snp") -> str:
    logic = get_hypothesis_logic(variant_type)
    
    feature_to_mechs = {}
    for name, rules in logic.items():
        for cond in rules['primary'] + rules.get('secondary', []):
            feat = cond['feature']
            feature_to_mechs.setdefault(feat, []).append(name)
    
    lines = []
    lines.append(f"FEATURE USAGE MAP — {variant_type.upper()}")
    lines.append("=" * 60)
    for feat in sorted(feature_to_mechs.keys()):
        mechs = feature_to_mechs[feat]
        lines.append(f"  {feat.replace('D_BED_', ''):<35} → {len(mechs)} mechanisms: {', '.join(mechs)}")
    
    lines.append(f"\nTotal unique features referenced: {len(feature_to_mechs)}")
    return "\n".join(lines)


def describe_tier5_differences() -> str:
    """Show what differs between SNP and INDEL Tier 5."""
    lines = []
    lines.append("=" * 72)
    lines.append("TIER 5 COMPARISON: SNP vs INDEL")
    lines.append("=" * 72)
    
    lines.append("")
    lines.append("SNP Tier 5:")
    for name, rules in _TIER5_SNP.items():
        lines.append(f"  {describe_mechanism(name, rules)}")
    
    lines.append("")
    lines.append("INDEL Tier 5:")
    for name, rules in _TIER5_INDEL.items():
        lines.append(f"  {describe_mechanism(name, rules)}")
    
    lines.append("")
    lines.append("KEY DIFFERENCES:")
    lines.append("  SNP:   Single-position LLR + MLM_Prior confidence check")
    lines.append("  INDEL: Window-level logprob_delta OR KL divergence (probability space)")
    lines.append("         + separate EMBEDDING_DISRUPTED (representation space)")
    
    return "\n".join(lines)


if __name__ == "__main__":
    for vt in ["snp", "indel"]:
        print(describe_all_rules(vt))
        print("\n")
        print(describe_feature_usage(vt))
        print("\n\n")
    
    print(describe_tier5_differences())
