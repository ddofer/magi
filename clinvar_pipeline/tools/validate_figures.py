#!/usr/bin/env python3
"""Validate figure inputs exist and (optionally) that figure PNGs were produced."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import PATHS

REQUIRED_INPUTS = {
    "snp deltas (benchmark)": PATHS["snp_deltas"],
    "indel deltas (benchmark)": PATHS["indel_deltas"],
    "pipeline snp annotated": PATHS["snp_annotated"],
    "pipeline indel annotated": PATHS["indel_annotated"],
    "pipeline snp signaled": PATHS["snp_annotated_signaled"],
    "pipeline indel signaled": PATHS["indel_annotated_signaled"],
    "snp impact": PATHS["snp_deltas_impact"],
    "indel impact": PATHS["indel_deltas_impact"],
    "snp llm canonical": PATHS["snp_llm_results"],
    "indel llm canonical": PATHS["indel_llm_results"],
}

EXPECTED_FIGURES = [
    PATHS["fig1c"],
    PATHS["fig2a"],
    PATHS["fig2b"],
    PATHS["fig3a"],
    PATHS["fig3b"],
    PATHS["fig3c"],
    PATHS["fig3d"],
    PATHS["fig3e"],
]


def validate(*, check_outputs: bool = False) -> dict:
    report: dict = {"status": "OK", "inputs": [], "outputs": []}

    for name, path in REQUIRED_INPUTS.items():
        ok = Path(path).exists()
        report["inputs"].append({"name": name, "path": path, "ok": bool(ok)})
        if not ok:
            report["status"] = "FAIL"

    if check_outputs:
        for path in EXPECTED_FIGURES:
            ok = Path(path).exists()
            report["outputs"].append({"path": path, "ok": bool(ok)})
            if not ok:
                report["status"] = "FAIL"

    return report


def main() -> None:
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-outputs", action="store_true")
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()

    report = validate(check_outputs=args.check_outputs)
    print(f"Figure validation: {report['status']}")
    for item in report["inputs"]:
        mark = "✓" if item["ok"] else "✗"
        print(f"  {mark} input  {item['name']}: {item['path']}")
    if args.check_outputs:
        for item in report["outputs"]:
            mark = "✓" if item["ok"] else "✗"
            print(f"  {mark} output {item['path']}")

    if args.write_report:
        out = Path(PATHS["validation_root"])
        out.mkdir(parents=True, exist_ok=True)
        (out / "figures_validation.json").write_text(json.dumps(report, indent=2))

    sys.exit(0 if report["status"] == "OK" else 1)


if __name__ == "__main__":
    main()
