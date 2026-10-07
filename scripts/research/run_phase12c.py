#!/usr/bin/env python3
"""CLI for Phase 12C: does local periodicity hold for the grounded 8 mm reflectarray cell?

Writes `results/phase12c_grounded_locality/` (Tasks A-D) and `results/phase12c_hfss_export/`
(Task E). Task A adds a real, new grounded-supercell capability (`phase12c.grounded_supercell`,
wrapping Phase 9's frozen `SupercellUnitCell`); Tasks B-D are measurement/audit only -- Task D
never edits Phase 9/10/11's own results. Geometry inversion / STEP export stay LOCKED regardless
of this phase's own findings.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from src.phase12c.run import print_final_report, run_phase12c  # noqa: E402


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "results" / "phase12c_grounded_locality"))
    args = parser.parse_args()

    result = run_phase12c(args.output_dir)
    print_final_report(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
