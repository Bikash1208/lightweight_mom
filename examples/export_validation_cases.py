"""Export the Phase 6 validation-case folders (STEP + JSON + CSV + plots + README, one folder
per case; see VALIDATION_EXPORT.md).

Run from the project root:
    python examples/export_validation_cases.py [results/validation_cases] [case_id ...]

With no case ids, exports all 10 cases from `src.validation_cases.CASES`. Needs the optional CAD
dependencies: `pip install -r requirements-cad.txt`.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.validation_cases import CASES, CASES_BY_ID
from src.validation_export import export_all_cases

ROOT = Path(__file__).resolve().parent.parent


def main(argv):
    args = list(argv)
    out_dir = ROOT / "results" / "validation_cases"
    if args and args[0] not in CASES_BY_ID:
        out_dir = Path(args.pop(0))
    unknown = [a for a in args if a not in CASES_BY_ID]
    if unknown:
        raise SystemExit(f"unknown case id(s) {unknown}; known: {sorted(CASES_BY_ID)}")
    cases = [CASES_BY_ID[i] for i in args] if args else list(CASES)

    print(f"Exporting {len(cases)} validation case(s) to {out_dir}")

    def progress(e):
        print(f"  [{e['done']}/{e['total']}] {e['case_id']}")

    t0 = time.time()
    results = export_all_cases(out_dir, cases=cases, progress=progress)
    dt = time.time() - t0

    print(f"\ndone in {dt:.2f}s")
    for case_id, r in results.items():
        n_metal = r.step_result.metal.n_solids if r.step_result.metal.n_solids else 0
        print(f"  {case_id}: {r.out_dir} "
              f"(metal volume {r.step_result.metal.volume_mm3:.4f} mm^3, "
              f"{len(r.database.records)} S-matrix record(s))")

    print(f"\nSee {out_dir / cases[0].case_id / 'README.md'} (per-case) and VALIDATION_EXPORT.md "
          f"(pipeline-wide) for how to reproduce any of these in HFSS. This project does not "
          f"perform that comparison automatically.")


if __name__ == "__main__":
    main(sys.argv[1:])
