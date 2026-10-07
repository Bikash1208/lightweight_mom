#!/usr/bin/env python3
"""Phase 8C CLI: unit-cell scattering architecture and control-quantity validation, against the
real Phase 7 data (read from the internal solve checkpoint, never re-solved).

    python scripts/research/run_phase8c_validation.py [--output-dir DIR] [--checkpoint PATH]

Does not build a GUI, does not modify the solver or Phase 7 data, and does not generate a final
MTS STEP file or calculate an aperture (Sec. 2-4, 57-60).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import argparse

from src.phase8c import data_source
from src.phase8c.run import print_final_report, run_phase8c


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="results/phase8c_scattering")
    parser.add_argument("--checkpoint", default=str(data_source.DEFAULT_CHECKPOINT_JSON))
    args = parser.parse_args(argv)

    result = run_phase8c(args.output_dir, checkpoint_path=args.checkpoint, progress=print)
    print()
    print_final_report(result)
    print()
    print(f"Output: {result['output_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
