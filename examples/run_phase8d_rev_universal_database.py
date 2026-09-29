#!/usr/bin/env python3
"""CLI for Phase 8D-REV: universal unit-cell database architecture.

Registers the existing, unmodified Phase 7 324-cell/4.8 mm library as the one benchmark
`LibraryEntry`, runs the Sec. 20 benchmark queries and Sec. 24-25 integrity checks, and writes
every required file under ``--output-dir`` (default `results/universal_database/`) plus
`UNIVERSAL_UNIT_CELL_DATABASE.md` at the project root.
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.phase8d_rev.run import print_final_report, run_phase8d_rev  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "results" / "universal_database"))
    parser.add_argument(
        "--regression-suite-passed",
        action="store_true",
        help="Set only after independently confirming `pytest` passes in full for this state of "
        "the repository; used solely to fill in Sec. 31's last acceptance-criteria line.",
    )
    args = parser.parse_args()

    result = run_phase8d_rev(
        args.output_dir,
        regression_suite_passed=args.regression_suite_passed,
    )
    print_final_report(result)
    return 0 if result["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
