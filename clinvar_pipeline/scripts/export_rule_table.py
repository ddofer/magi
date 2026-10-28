#!/usr/bin/env python3
"""Write the rule layer (SNV and indel rules) for the supplement, generated from the rule definitions.

rule_table.csv: one row per condition (machine-readable).
rule_table_paper.csv: one row per rule in checking order, conditions written out, rules shared by
SNVs and indels merged, and a last row for variants where no rule fires. Caption: clinvar_pipeline/README.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "vendor"))

from config import OUTPUT_ROOT
from nt_mechanism_assignment import rule_table

RULE_TABLE = OUTPUT_ROOT / "mechanisms" / "rule_table.csv"
NAMES = {  # feature labels; other annotation elements are written as "Δ <element>"
    "LLR": "LLR = ln(P alt / P ref)",
    "MLM_Prior": "Reference-letter probability",
    "MLM_KL_max": "Max KL(Alt ‖ Ref)",
    "MLM_logprob_delta": "Mean log-probability change (alt − ref, z-score)",
    "EMB_cosine_dist": "Embedding cosine distance",
    "D_BED_always_on_exon": "Δ constitutive exon",
    "D_BED_protein_coding_gene": "Δ protein-coding gene",
    "D_BED_promoter_Tissue_specific": "Δ tissue-specific promoter",
    "D_BED_promoter_Tissue_invariant": "Δ tissue-invariant promoter",
    "D_BED_enhancer_Tissue_specific": "Δ tissue-specific enhancer",
    "D_BED_enhancer_Tissue_invariant": "Δ tissue-invariant enhancer",
    "D_BED_polyA_signal": "Δ polyadenylation signal",
    "D_BED_CTCF-bound": "Δ CTCF-bound site",
}
COLUMNS = {
    "check": "Check order", "tier": "Tier", "mechanism": "Mechanism (NT_Mechanism)",
    "applies_to": "Applies to", "any_of": "Any of (OR)", "all_of": "All of (AND)",
    "description": "Description", "pre_empted_by": "Never assigned: always pre-empted by",
}


def _condition(row) -> str:
    """e.g. 'Δ splice donor ≤ −0.2' or '|Δ exon| ≥ 0.08'."""
    f = row.feature
    name = NAMES.get(f) or "Δ " + f.removeprefix("D_BED_").replace("_", " ")
    if row.direction == "any":
        return f"|{name}| ≥ {abs(row.threshold):g}"
    op = "≤" if row.direction == "negative" else "≥"
    return f"{name} {op} {row.threshold:g}".replace(" -", " −")  # typographic minus


def paper_table(table: pd.DataFrame) -> pd.DataFrame:
    """One row per rule in checking order; a rule identical for SNVs and indels appears once."""
    cond = table.apply(_condition, axis=1)
    primary = table["condition_role"].str.startswith("primary")
    table = table.assign(any_of=cond.where(primary), all_of=cond.mask(primary),
                         applies_to=table["variant_type"].map({"snp": "SNV", "indel": "Indel"}))
    keys = ["mechanism", "tier", "description", "pre_empted_by"]
    join = lambda c: "; ".join(c.dropna())  # noqa: E731
    rules = table.groupby(["applies_to", "order_within_tier", *keys], sort=False).agg(
        any_of=("any_of", join), all_of=("all_of", join)).reset_index()
    keys += ["any_of", "all_of"]
    rules["applies_to"] = rules.groupby(keys)["applies_to"].transform(lambda v: "Both" if len(v) > 1 else v.iloc[0])
    rules = rules.drop_duplicates(keys).sort_values(["tier", "order_within_tier"], kind="stable").reset_index(drop=True)
    rules["check"] = range(1, len(rules) + 1)
    rules.loc[len(rules)] = {"check": "", "tier": "—", "mechanism": "UNCERTAIN_SIGNIFICANCE", "applies_to": "Both",
                             "any_of": "No rule above fires", "description": "No rule fired (not a ClinVar VUS)"}
    return rules[list(COLUMNS)].rename(columns=COLUMNS).fillna("")


def write_rule_table(out=RULE_TABLE) -> None:
    table = pd.concat([rule_table("snp"), rule_table("indel")], ignore_index=True)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False)
    paper = paper_table(table)
    paper_out = out.with_name(f"{out.stem}_paper{out.suffix}")
    paper.to_csv(paper_out, index=False, encoding="utf-8-sig")  # utf-8-sig: Excel shows ≤ and Δ
    print(f"Saved {len(table)} conditions to {out} and {len(paper)} rows to {paper_out}")
    blocked = table.loc[~table["can_be_assigned"], ["variant_type", "mechanism", "pre_empted_by"]].drop_duplicates()
    if len(blocked):
        print("Rules that can never be the assigned mechanism:")
        print(blocked.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(RULE_TABLE))
    write_rule_table(parser.parse_args().out)


if __name__ == "__main__":
    main()
