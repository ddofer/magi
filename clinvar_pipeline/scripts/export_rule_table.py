#!/usr/bin/env python3
"""Write the rule layer (SNV and indel rules) for the supplement.

rule_table.csv: one row per condition (machine-readable).
rule_table_paper.csv: one row per rule, conditions written out, rules shared by SNVs and indels merged.
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
FEATURE_NAMES = {  # sequence-model metrics; annotation elements are written as "Δ <element>"
    "LLR": "LLR",
    "MLM_Prior": "Ref letter probability",
    "MLM_KL_max": "Max KL divergence",
    "MLM_logprob_delta": "Δ log-prob (z-score)",
    "EMB_cosine_dist": "Embedding cosine distance",
}
ELEMENT_NAMES = {  # annotation elements whose model names read badly
    "always_on_exon": "constitutive exon",
    "protein_coding_gene": "protein-coding gene",
    "promoter_Tissue_specific": "tissue-specific promoter",
    "promoter_Tissue_invariant": "tissue-invariant promoter",
    "enhancer_Tissue_specific": "tissue-specific enhancer",
    "enhancer_Tissue_invariant": "tissue-invariant enhancer",
}


def _condition(row) -> str:
    """e.g. 'Δ splice donor ≤ -0.20' or '|Δ exon| ≥ 0.08'."""
    f = row.feature
    element = f.removeprefix("D_BED_")
    name = FEATURE_NAMES.get(f, "Δ " + ELEMENT_NAMES.get(element, element.replace("_", " ")))
    if row.direction == "any":
        return f"|{name}| ≥ {abs(row.threshold):g}"
    return f"{name} {'≤' if row.direction == 'negative' else '≥'} {row.threshold:g}"


def paper_table(table: pd.DataFrame) -> pd.DataFrame:
    """One row per rule, in assignment order; a rule identical for SNVs and indels appears once.

    Caption facts: Δ is Alt minus Ref on MAGI's sigmoid scale; the lowest tier wins, then the
    lower order in tier; no rule firing gives UNCERTAIN_SIGNIFICANCE. Mechanism names are the
    NT_Mechanism values.
    """
    table = table.assign(condition=table.apply(_condition, axis=1))
    rules = table.groupby(["variant_type", "mechanism"], sort=False).agg(
        tier=("tier", "first"),
        order=("order_within_tier", "first"),
        any_of=("condition", lambda c: "; ".join(c[table.loc[c.index, "condition_role"].str.startswith("primary")])),
        all_of=("condition", lambda c: "; ".join(c[table.loc[c.index, "condition_role"].str.startswith("secondary")])),
        description=("description", "first"),
        pre_empted_by=("pre_empted_by", "first"),
    ).reset_index()
    keys = ["mechanism", "tier", "any_of", "all_of", "description", "pre_empted_by"]
    rules["applies_to"] = rules.groupby(keys)["variant_type"].transform(
        lambda v: ", ".join({"snp": "SNV", "indel": "Indel"}[x] for x in v)
    )
    rules = rules.drop_duplicates(keys).sort_values(["tier", "order"], kind="stable")
    rules["order"] += 1
    return rules.rename(columns={
        "tier": "Tier", "order": "Order in tier", "mechanism": "Mechanism", "applies_to": "Applies to",
        "any_of": "Fires when any of", "all_of": "And all of", "description": "Description",
        "pre_empted_by": "Never assigned: always pre-empted by",
    })[["Tier", "Order in tier", "Mechanism", "Applies to", "Fires when any of", "And all of",
        "Description", "Never assigned: always pre-empted by"]]


def write_rule_table(out=RULE_TABLE) -> None:
    table = pd.concat([rule_table("snp"), rule_table("indel")], ignore_index=True)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False)
    paper = paper_table(table)
    paper_out = out.with_name(f"{out.stem}_paper{out.suffix}")
    paper.to_csv(paper_out, index=False, encoding="utf-8-sig")  # utf-8-sig: Excel shows ≤ and Δ
    print(f"Saved {len(table)} conditions to {out} and {len(paper)} rules to {paper_out}")
    blocked = table.loc[~table["can_be_primary"], ["variant_type", "mechanism", "pre_empted_by"]].drop_duplicates()
    if len(blocked):
        print("Rules that can never be the assigned mechanism:")
        print(blocked.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(RULE_TABLE))
    write_rule_table(parser.parse_args().out)


if __name__ == "__main__":
    main()
