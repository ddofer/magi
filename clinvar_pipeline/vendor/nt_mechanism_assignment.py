import pandas as pd
import numpy as np

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


# =============================================================================
# 1c. BUILD COMPLETE LOGIC PER VARIANT TYPE
# =============================================================================

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


# =============================================================================
# 2. CORE FUNCTIONS — VECTORIZED (no df.apply, no row-by-row iteration)
# =============================================================================

def _check_condition_vec(series, threshold, direction):
    """Vectorized condition check. Returns boolean Series."""
    s = series.fillna(0)
    if direction == "negative":
        return s <= threshold
    elif direction == "positive":
        return s >= threshold
    elif direction == "any":
        return s.abs() >= abs(threshold)
    return pd.Series(False, index=series.index)


def _get_column_safe(df, feature):
    """Get a column from df, returning zeros if missing."""
    if feature in df.columns:
        return df[feature].fillna(0)
    return pd.Series(0, index=df.index, dtype=float)


def assign_mechanisms_vectorized(df, hypothesis_logic):
    """
    Vectorized mechanism assignment across the entire dataframe.
    
    For each mechanism, evaluates primary (OR) and secondary (AND) conditions
    using numpy/pandas vectorized ops. No row-by-row iteration.
    
    Adds columns: NT_Mechanism, NT_Description, NT_Trigger, 
                  NT_All_Mechanisms, NT_Mechanism_Count
    """
    n = len(df)
    sorted_hypotheses = sorted(hypothesis_logic.items(), key=lambda x: x[1]['priority'])
    
    # Pre-compute: for each mechanism, a boolean mask of matching rows
    # and the trigger info (feature name + value for the first matching primary)
    mech_names = []
    mech_masks = []
    mech_priorities = []
    mech_descriptions = []
    mech_trigger_features = []  # Series of trigger feature name per row
    mech_trigger_values = []    # Series of trigger value per row
    
    for mech_name, rules in sorted_hypotheses:
        # PRIMARY (OR): build mask for each primary condition, then OR them
        primary_masks = []
        primary_features = []
        for cond in rules['primary']:
            col = _get_column_safe(df, cond['feature'])
            mask = _check_condition_vec(col, cond['threshold'], cond['direction'])
            primary_masks.append(mask)
            primary_features.append((cond['feature'], col))
        
        if not primary_masks:
            continue
            
        primary_or = primary_masks[0]
        for m in primary_masks[1:]:
            primary_or = primary_or | m
        
        # SECONDARY (AND): all must pass
        secondary_and = pd.Series(True, index=df.index)
        for cond in rules.get('secondary', []):
            col = _get_column_safe(df, cond['feature'])
            mask = _check_condition_vec(col, cond['threshold'], cond['direction'])
            secondary_and = secondary_and & mask
        
        combined = primary_or & secondary_and
        
        if combined.sum() == 0:
            continue
        
        # Determine which primary condition triggered (first match)
        # Initialize with first condition
        trigger_feat = pd.Series("", index=df.index)
        trigger_val = pd.Series(0.0, index=df.index)
        assigned = pd.Series(False, index=df.index)
        
        for mask_i, (feat_name, feat_col) in zip(primary_masks, primary_features):
            new_assign = mask_i & combined & ~assigned
            trigger_feat = trigger_feat.where(~new_assign, feat_name.replace('D_BED_', ''))
            trigger_val = trigger_val.where(~new_assign, feat_col)
            assigned = assigned | new_assign
        
        mech_names.append(mech_name)
        mech_masks.append(combined)
        mech_priorities.append(rules['priority'])
        mech_descriptions.append(rules.get('description', ''))
        mech_trigger_features.append(trigger_feat)
        mech_trigger_values.append(trigger_val)
    
    # Now assign: for each row, collect all matching mechanisms
    # Primary mechanism = first match (sorted by priority already)
    
    # Initialize output columns
    nt_mechanism = pd.Series('UNCERTAIN_SIGNIFICANCE', index=df.index)
    nt_description = pd.Series('No NT signal above noise floor', index=df.index)
    nt_trigger = pd.Series('', index=df.index)
    nt_all_mechs = pd.Series('', index=df.index)
    nt_mech_count = pd.Series(0, index=df.index, dtype=int)
    
    # Track which rows already have a primary assignment
    has_primary = pd.Series(False, index=df.index)
    
    # Build all_mechanisms strings and assign primary
    # We need to iterate mechanisms (not rows!) — this is O(mechanisms), not O(rows)
    for i, (mname, mask, prio, desc, tfeat, tval) in enumerate(
        zip(mech_names, mech_masks, mech_priorities, mech_descriptions,
            mech_trigger_features, mech_trigger_values)
    ):
        # Count: every matching row gets +1
        nt_mech_count = nt_mech_count + mask.astype(int)
        
        # Append to all_mechanisms string
        mech_tag = f"{mname}(P{prio})"
        current_all = nt_all_mechs[mask]
        nt_all_mechs[mask] = current_all.where(
            current_all == '', current_all + '; '
        ) + mech_tag
        
        # Primary assignment: first mechanism to match wins
        new_primary = mask & ~has_primary
        if new_primary.sum() > 0:
            nt_mechanism[new_primary] = mname
            nt_description[new_primary] = desc
            nt_trigger[new_primary] = tfeat[new_primary] + '=' + tval[new_primary].map(lambda v: f"{v:.3f}")
            has_primary = has_primary | new_primary
    
    df['NT_Mechanism'] = nt_mechanism
    df['NT_Description'] = nt_description
    df['NT_Trigger'] = nt_trigger
    df['NT_All_Mechanisms'] = nt_all_mechs
    df['NT_Mechanism_Count'] = nt_mech_count
    
    return df


def extract_top_signals_vectorized(df, delta_cols, k=3):
    """
    Vectorized extraction of top K signals by magnitude, gain, and loss.
    
    Uses numpy argsort on the delta matrix instead of row-by-row iteration.
    """
    mat = df[delta_cols].fillna(0).values  # (n_variants, n_features)
    col_names = np.array([c.replace('D_BED_', '') for c in delta_cols])
    n = mat.shape[0]
    
    def _format_topk(indices, values):
        """Format top-k as string column."""
        results = []
        for i in range(n):
            parts = []
            for j in range(min(k, len(indices[i]))):
                idx = indices[i][j]
                val = values[i][idx] if isinstance(idx, (int, np.integer)) else values[i][j]
                parts.append(f"{col_names[indices[i][j]]}={mat[i, indices[i][j]]:.3f}")
            results.append("; ".join(parts) if parts else "None")
        return results
    
    # Top by absolute magnitude
    abs_mat = np.abs(mat)
    abs_topk = np.argsort(-abs_mat, axis=1)[:, :k]
    
    abs_strs = []
    for i in range(n):
        parts = [f"{col_names[j]}={mat[i, j]:.3f}" for j in abs_topk[i]]
        abs_strs.append("; ".join(parts))
    
    # Top gains (positive values, descending)
    gain_strs = []
    for i in range(n):
        pos_idx = np.where(mat[i] > 0.01)[0]
        if len(pos_idx) == 0:
            gain_strs.append("None")
        else:
            sorted_pos = pos_idx[np.argsort(-mat[i, pos_idx])][:k]
            gain_strs.append("; ".join(f"{col_names[j]}={mat[i,j]:.3f}" for j in sorted_pos))
    
    # Top losses (negative values, ascending)
    loss_strs = []
    for i in range(n):
        neg_idx = np.where(mat[i] < -0.01)[0]
        if len(neg_idx) == 0:
            loss_strs.append("None")
        else:
            sorted_neg = neg_idx[np.argsort(mat[i, neg_idx])][:k]
            loss_strs.append("; ".join(f"{col_names[j]}={mat[i,j]:.3f}" for j in sorted_neg))
    
    df['Top_Abs_Signals'] = abs_strs
    df['Top_Gain_Signals'] = gain_strs
    df['Top_Loss_Signals'] = loss_strs
    
    return df


# =============================================================================
# 3. MAIN EXECUTION
# =============================================================================

def run_assignment(input_file, output_file, variant_type):
    """
    Run mechanism assignment pipeline for a given variant type.
    
    Parameters
    ----------
    input_file : str
        Path to input parquet.
    output_file : str
        Path to output CSV.
    variant_type : str
        "snp" or "indel".
    """
    vt = variant_type.strip().lower()
    logic = get_hypothesis_logic(vt)
    
    print(f"{'='*60}")
    print(f"  MECHANISM ASSIGNMENT — {vt.upper()}")
    print(f"{'='*60}")
    
    print(f"Loading {input_file}...")
    df = pd.read_parquet(input_file)
    print(f"Loaded {len(df)} variants")
    
    delta_cols = [c for c in df.columns if c.startswith("D_BED_")]
    print(f"Found {len(delta_cols)} delta columns")
    
    df[delta_cols] = df[delta_cols].apply(pd.to_numeric, errors='coerce')
    
    # --- ASSIGN MECHANISMS (vectorized) ---
    print("\nAssigning NT mechanisms...")
    df = assign_mechanisms_vectorized(df, logic)
    
    # --- EXTRACT TOP SIGNALS (vectorized) ---
    print("Extracting top signals...")
    df = extract_top_signals_vectorized(df, delta_cols)
    
    # --- SUMMARY ---
    print(f"\n{'='*60}")
    print("MECHANISM ASSIGNMENT SUMMARY")
    print(f"{'='*60}")
    print(df['NT_Mechanism'].value_counts())
    
    print(f"\nPriority tier distribution:")
    tier_map = {name: rules['priority'] for name, rules in logic.items()}
    df['_tier'] = df['NT_Mechanism'].map(tier_map).fillna(-1).astype(int)
    print(df['_tier'].value_counts().sort_index())
    df.drop(columns=['_tier'], inplace=True)
    
    # Coverage
    n_assigned = (df['NT_Mechanism'] != 'UNCERTAIN_SIGNIFICANCE').sum()
    print(f"\nCoverage: {n_assigned}/{len(df)} ({100*n_assigned/len(df):.1f}%) variants assigned a mechanism")
    
    df.to_csv(output_file)
    print(f"Saved to: {output_file}")
    
    return df


def main():
    # --- SNPs ---
    df_snp = run_assignment(
        input_file="parquet/annotated_snps_signaled.parquet",
        output_file="snps_assigned_mechanisms.csv",
        variant_type="snp"
    )
    
    print("\n")
    
    # --- INDELs ---
    df_indel = run_assignment(
        input_file="parquet/annotated_indels_signaled.parquet",
        output_file="indels_assigned_mechanisms.csv",
        variant_type="indel"
    )
    
    return df_snp, df_indel


if __name__ == "__main__":
    df_snp, df_indel = main()
