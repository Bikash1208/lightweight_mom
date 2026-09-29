#!/usr/bin/env python3
"""CLI for Phase 8E: physical control-variable and holographic-scattering bridge.

Writes `results/phase8e_control_bridge/` (report.md, JSON/CSV outputs, plots/). Read-only with
respect to Phase 7 and Phase 8D-REV; performs no geometry inversion and generates no STEP file.
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.phase8e.run import print_final_report, run_phase8e  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(ROOT / "results" / "phase8e_control_bridge"))
    args = parser.parse_args()

    result = run_phase8e(args.output_dir)
    print_final_report(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
