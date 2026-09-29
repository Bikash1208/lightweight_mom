"""Write the Phase 6A HFSS comparison export (5 plots + hfss_comparison_reference.csv +
HFSS_COMPARISON.md) for an already-exported Phase 6 validation case.

Run from the project root, after examples/export_validation_cases.py has produced the case
folder:
    python examples/export_hfss_comparison.py [results/validation_cases/case01_non_resonant_patch]

Reads only that folder's mom_results.json / geometry.json / solver_config.json; never solves.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.hfss_comparison import export_hfss_comparison

ROOT = Path(__file__).resolve().parent.parent


def main(argv):
    case_dir = Path(argv[0]) if argv else ROOT / "results" / "validation_cases" / "case01_non_resonant_patch"
    result = export_hfss_comparison(case_dir)

    print(f"case folder: {result.case_dir}")
    print(f"  {result.csv_path}")
    print(f"  {result.md_path}")
    for p in result.plot_paths.values():
        print(f"  {p}")
    print("\nThis is a manual comparison reference only; no HFSS comparison was performed.")


if __name__ == "__main__":
    main(sys.argv[1:])
