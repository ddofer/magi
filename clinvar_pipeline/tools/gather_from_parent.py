#!/usr/bin/env python3
"""Deprecated: use tools/sync_inputs.py + tools/gather_canonical_llm.py."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.gather_canonical_llm import gather as gather_llm
from tools.sync_inputs import sync


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    print("Note: gather_from_parent.py is deprecated; syncing inputs + canonical LLM only.")
    sync(force=args.force)
    gather_llm(force=args.force)


if __name__ == "__main__":
    main()
