#!/usr/bin/env python3
"""CLI for Phase 9: physical control-variable + local-periodicity validation.

Writes `results/phase9_validation/`. Read-only with respect to Phase 7/8; performs no geometry
inversion and generates no final MTS STEP (validation-only diagnostic supercell geometries are
solved in memory, never exported).
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.phase9.run import print_final_report, run_phase9  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "results" / "phase9_validation"))
    parser.add_argument("--hfss-reference", default=None, help="Path to an HFSS reference CSV, if one exists.")
    args = parser.parse_args()

    result = run_phase9(args.output_dir, hfss_reference_path=args.hfss_reference)
    print_final_report(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
