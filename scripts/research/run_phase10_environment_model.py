#!/usr/bin/env python3
"""CLI for Phase 10: environment-aware electromagnetic response model.

Writes `results/phase10_environment_model/`. Does not modify Phase 1-9; performs no geometry
inversion beyond the first-stage, MoM-verified test harness in `src/phase10/inverse.py`, and
generates no final MTS STEP.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from src.phase10.run import print_final_report, run_phase10  # noqa: E402


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "results" / "phase10_environment_model"))
    parser.add_argument("--hfss-reference", default=None, help="Path to an HFSS reference CSV, if one exists.")
    args = parser.parse_args()

    result = run_phase10(args.output_dir, hfss_reference_path=args.hfss_reference)
    print_final_report(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
