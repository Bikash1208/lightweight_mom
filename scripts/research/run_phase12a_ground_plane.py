#!/usr/bin/env python3
"""CLI for Phase 12A: optional PEC ground plane (grounded-substrate stack).

Writes `results/phase12a_ground_plane/`. Adds a capability (an opt-in `backing="pec"` code path
across the core modules) and measures what it gives; does not build a production library, does
not re-run Phases 8-11, and does not unlock geometry inversion or STEP export.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from src.phase12a.run import print_final_report, run_phase12a  # noqa: E402


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "results" / "phase12a_ground_plane"))
    args = parser.parse_args()

    result = run_phase12a(args.output_dir)
    print_final_report(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
