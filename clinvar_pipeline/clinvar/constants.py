import numpy as np

CLINVAR_SUMMARY_COLS = [
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "VariationID",
    "Type",
    "GeneSymbol",
    "Assembly",
    "ReviewStatus",
]

SUBMISSION_COLS = [
    "#VariationID",
    "Description",
    "ExplanationOfInterpretation",
    "ReviewStatus",
]

ALLOWED_OTHER_STATUS = {
    "no assertion criteria provided",
    "no classification provided",
}
CRITERIA_SINGLE = "criteria provided, single submitter"

FEATURE_PRIORITY = [
    "start_codon",
    "stop_codon",
    "CDS",
    "UTR5",
    "UTR3",
    "UTR",
    "exon",
    "transcript",
    "gene",
]

REVIEW_STATUS_TO_GOLD_STARS = {
    "criteria provided, single submitter": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": np.nan,
    "reviewed by expert panel": 3,
    "no assertion provided": np.nan,
    "no interpretation for the single variant": np.nan,
    "practice guideline": 4,
}

Z_COLS_SNP = ["MLM_logprob_ref", "MLM_logprob_alt", "MLM_logprob_delta"]

Z_COLS_INDEL = [
    "MLM_logprob_ref",
    "MLM_logprob_alt",
    "MLM_logprob_delta",
    "EMB_l2_dist",
    "EMB_max_pos_dist",
    "EMB_mean_pos_dist",
]

PLACEHOLDER_RATIONALES = {
    "No rationale provided.",
    "No detailed rationale provided.",
}

ANNOTATION_COLUMNS = [
    "gene",
    "mRNA",
    "mRNA_promoter",
    "mRNA_exon",
    "coding_sequence",
    "start_codon",
    "stop_codon",
    "five_prime_UTR",
    "three_prime_UTR",
    "mRNA_intron",
    "mRNA_splice",
    "lncRNA",
    "lncRNA_promoter",
    "lncRNA_exon",
    "snRNA",
    "snRNA_promoter",
    "snRNA_exon",
    "antisenseRNA",
    "antisenseRNA_promoter",
    "antisenseRNA_exon",
    "telomeraseRNA",
    "telomeraseRNA_promoter",
    "telomeraseRNA_exon",
    "RNaseMRPRNA",
    "RNaseMRPRNA_promoter",
    "RNaseMRPRNA_exon",
    "snoRNA",
    "snoRNA_promoter",
    "snoRNA_exon",
    "other",
]

RNA_TYPES = [
    "lncRNA",
    "snRNA",
    "antisenseRNA",
    "telomeraseRNA",
    "RNaseMRPRNA",
    "snoRNA",
]
