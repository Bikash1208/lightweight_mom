#!/usr/bin/env python3
"""CLI for Phase 12B: mesh alignment, duplicate-geometry audit, grounded period study.

Writes `results/phase12b_mesh_period_study/`. Task A (mesh alignment / duplicate detection) is a
core, additive capability in `discretization.py`/`sweep.py`; Tasks B-D are measurement/audit only
-- no earlier phase's results are modified in place, and geometry inversion / STEP export stay
LOCKED regardless of this phase's own findings.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from src.phase12b.run import print_final_report, run_phase12b  # noqa: E402


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "results" / "phase12b_mesh_period_study"))
    args = parser.parse_args()

    result = run_phase12b(args.output_dir)
    print_final_report(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
