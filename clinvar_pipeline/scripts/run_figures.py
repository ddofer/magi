"""Run all publication figures into clinvar_pipeline/output/figures/."""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from figures import (
    fig1_pies,
    fig2_densities,
    fig2a_impact_vs_af,
    fig3_abc,
    fig3d_indel_frame,
    fig3e_mechanism_concordance,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fig",
        choices=["1", "2a", "2", "3abc", "3d", "3e", "all"],
        default="all",
    )
    parser.add_argument("--assign-mechanisms", action="store_true", help="Re-run NT mechanism assignment for fig3e")
    args = parser.parse_args()

    runners = {
        "1": fig1_pies.run,
        "2a": fig2a_impact_vs_af.run,
        "2": fig2_densities.run,
        "3abc": fig3_abc.run,
        "3d": fig3d_indel_frame.run,
        "3e": lambda: fig3e_mechanism_concordance.run(assign=args.assign_mechanisms),
    }
    if args.fig == "all":
        todo = ["1", "2a", "2", "3abc", "3d", "3e"]
    else:
        todo = [args.fig]
    print("══ Generating figures ══")
    failed: list[str] = []
    for key in todo:
        print(f"\n── Fig {key} ──")
        try:
            runners[key]()
        except Exception:
            failed.append(key)
            traceback.print_exc()
    if failed:
        print(f"\n✗ Failed figures: {', '.join(failed)}")
        sys.exit(1)
    print("\n✓ Figures complete.")


if __name__ == "__main__":
    main()
