#!/usr/bin/env python3
"""Bootstrap data/, validate MANE, run analysis steps, validate figures. Writes run.log."""

from __future__ import annotations

import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "output" / "validation" / "run.log"


def run(cmd: list[str]) -> int:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as f:
        f.write(f"\n$ {' '.join(cmd)}\n")
        p = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        return p.returncode


def main() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text(f"=== run_validation {datetime.now().isoformat()} ===\n")
    steps = [
        [sys.executable, "tools/sync_inputs.py", "--force"],
        [sys.executable, "tools/gather_canonical_llm.py", "--force"],
        [sys.executable, "tools/validate_mane.py", "--rebuild", "--write-report"],
        [sys.executable, "scripts/06_compute_impact_scores.py"],
        [sys.executable, "scripts/07_assign_mechanisms.py"],
        [sys.executable, "tools/compute_cohort_funnel.py"],
        [sys.executable, "scripts/run_figures.py", "--fig", "all"],
        [sys.executable, "tools/validate_figures.py", "--check-outputs", "--write-report"],
    ]
    for cmd in steps:
        rc = run(cmd)
        if rc != 0:
            print(f"FAILED: {' '.join(cmd)} (see {LOG})")
            sys.exit(rc)
    print(f"OK — log: {LOG}")


if __name__ == "__main__":
    main()
