#!/usr/bin/env python3
"""Phase 8B CLI: evaluate the holographic formulation and physical control-variable question
against the real Phase 7 library, and write `results/phase8b_formulation/`.

    python examples/run_phase8b_formulation.py [--library-dir DIR] [--output-dir DIR]
                                                [--design-frequency-ghz F]

Does not build a GUI, does not modify the solver or Phase 7 data, and does not generate a final
MTS STEP file (Sec. 2-4, 49).
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.feasibility_schema import DEFAULT_LIBRARY_DIR
from src.phase8b.run import print_final_report, run_phase8b


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library-dir", default=str(DEFAULT_LIBRARY_DIR))
    parser.add_argument("--output-dir", default="results/phase8b_formulation")
    parser.add_argument("--design-frequency-ghz", type=float, default=12.5)
    args = parser.parse_args(argv)

    result = run_phase8b(args.output_dir, library_dir=args.library_dir, design_frequency_ghz=args.design_frequency_ghz)
    print_final_report(result)
    print()
    print(f"Output: {result['output_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
