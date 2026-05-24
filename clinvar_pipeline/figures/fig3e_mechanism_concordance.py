"""Fig 3e — concordance by NT-assigned mechanism (signal_mechanism_concordance_fig.ipynb)."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(PIPELINE_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))

from config import PATHS
from figures.common import ensure_dirs
from figures.merge_eval import _read_llm_results
from figures.plotting import plot_concordance_by_mechanism
from nt_mechanism_assignment import run_assignment


def _attach_mechanisms(llm_df: pd.DataFrame, mech_csv: str) -> pd.DataFrame:
    mech = pd.read_csv(mech_csv, index_col=0)
    if "#VariationID" not in mech.columns:
        mech = mech.reset_index()
    lookup = mech.set_index("#VariationID")["NT_Mechanism"]
    out = llm_df.copy()
    out["NT_Mechanism"] = out["#VariationID"].map(lookup).fillna("UNCERTAIN_SIGNIFICANCE")
    return out


def run(*, assign: bool = False) -> None:
    ensure_dirs()
    if assign or not Path(PATHS["snp_mechanisms_csv"]).exists():
        run_assignment(PATHS["snp_annotated_signaled"], PATHS["snp_mechanisms_csv"], "snp")
        run_assignment(PATHS["indel_annotated_signaled"], PATHS["indel_mechanisms_csv"], "indel")

    snp_res = _attach_mechanisms(_read_llm_results("snp"), PATHS["snp_mechanisms_csv"])
    indel_res = _attach_mechanisms(_read_llm_results("indel"), PATHS["indel_mechanisms_csv"])

    plot_concordance_by_mechanism(
        snp_res,
        "snp",
        Path(PATHS["fig3e"]),
        Path(PATHS["fig3e"]).with_suffix(""),
        min_samples=20,
    )
    plot_concordance_by_mechanism(
        indel_res,
        "indel",
        Path(PATHS["fig3e_indel"]),
        Path(PATHS["fig3e_indel"]).with_suffix(""),
        min_samples=20,
    )
    print(f"  Fig3e → {PATHS['fig3e']} (+ indel: {PATHS['fig3e_indel']})")
