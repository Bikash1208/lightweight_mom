#!/usr/bin/env python3
"""Phase 12D CLI: resolve the local-periodicity step limit, and ingest the real HFSS data.

Usage:
    .venv/bin/python scripts/research/run_phase12d.py [--output-dir DIR]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from src.phase12d.run import print_final_report, run_phase12d

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "results" / "phase12d_local_periodicity"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    result = run_phase12d(args.output_dir)
    print_final_report(result)
    print(f"\nOutputs written to {result.output_dir}")


if __name__ == "__main__":
    main()
