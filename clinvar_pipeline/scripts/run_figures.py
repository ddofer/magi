"""Run all publication figures into clinvar_pipeline/output/figures/."""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _runners(selected: list[str], *, assign_mechanisms: bool = False) -> dict[str, callable]:
    runners: dict[str, callable] = {}

    if "1" in selected:
        from figures import fig1_pies

        runners["1"] = fig1_pies.run
    if "2a" in selected:
        from figures import fig2a_impact_vs_af

        runners["2a"] = fig2a_impact_vs_af.run
    if "2" in selected:
        from figures import fig2_densities

        runners["2"] = fig2_densities.run
    if "3abc" in selected:
        from figures import fig3_abc

        runners["3abc"] = fig3_abc.run
    if "3d" in selected:
        from figures import fig3d_indel_frame

        runners["3d"] = fig3d_indel_frame.run
    if "3e" in selected:
        from figures import fig3e_mechanism_concordance

        runners["3e"] = lambda: fig3e_mechanism_concordance.run(assign=assign_mechanisms)
    if "s1" in selected:
        from figures import figs1_animals_concordance

        runners["s1"] = figs1_animals_concordance.run

    return runners


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fig",
        choices=["1", "2a", "2", "3abc", "3d", "3e", "s1", "all"],
        default="all",
    )
    parser.add_argument("--assign-mechanisms", action="store_true", help="Re-run NT mechanism assignment for fig3e")
    args = parser.parse_args()

    if args.fig == "all":
        todo = ["1", "2a", "2", "3abc", "3d", "3e", "s1"]
    else:
        todo = [args.fig]
    runners = _runners(todo, assign_mechanisms=args.assign_mechanisms)
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
