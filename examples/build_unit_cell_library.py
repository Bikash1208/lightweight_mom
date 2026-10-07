"""Build the Phase 7 production unit-cell library (18x18 L/W sweep, 201 frequencies, the locked
configuration in `src/library_schema.py::LibraryConfig.default`).

Run from the project root:
    python examples/build_unit_cell_library.py [results/unit_cell_library]

This is a long build (hours, not minutes -- see LIBRARY_REPORT.md's Sec. 17 and the project's
own measured mesh/Floquet cost trade-off); it writes an append-only checkpoint under
`<out>/_solve/checkpoint.jsonl` as it goes (`runner.run_sweep`, Phase 5, unmodified), so an
interrupted run resumes exactly where it stopped when re-run with the same command -- nothing is
re-solved, and nothing already written is overwritten.

Needs the optional CAD/plotting dependencies: `pip install -r requirements-cad.txt`.
"""

import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.filterwarnings("ignore", category=UserWarning)  # grating-order warning: none expected in this band

from src.library_schema import LibraryConfig  # noqa: E402
from src.unit_cell_library import build_full_library  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def main(argv):
    out_dir = Path(argv[0]) if argv else ROOT / "results" / "unit_cell_library"
    config = LibraryConfig.default()

    n_freq = int(round((config.f_max_hz - config.f_min_hz) / config.f_step_hz)) + 1
    total = config.n_geometries * n_freq
    print(
        f"Building unit-cell library: {config.n_geometries} geometries x {n_freq} frequencies "
        f"= {total} solves -> {out_dir}"
    )

    t0 = time.time()
    last = [0]

    def progress(e):
        if e["done"] - last[0] >= 1000 or e["done"] == e["total"]:
            last[0] = e["done"]
            elapsed = time.time() - t0
            rate = e["done"] / elapsed if elapsed > 0 else 0.0
            eta = (e["total"] - e["done"]) / rate if rate > 0 else float("nan")
            print(
                f"  {e['done']}/{e['total']} ({100 * e['done'] / e['total']:.1f}%) "
                f"elapsed={elapsed:.0f}s eta={eta:.0f}s",
                flush=True,
            )

    result = build_full_library(config, out_dir, progress=progress)

    print(f"\ndone in {time.time() - t0:.0f}s")
    print(f"  {result.build.n_ok_cells} solved OK, {result.build.n_invalid_cells} invalid/skipped")
    print(f"  computed this run: {result.build.computed}, reused from checkpoint: {result.build.reused}")
    print(
        f"  phase coverage: {result.coverage.coverage_deg:.2f} deg of 360 deg "
        f"(largest gap {result.coverage.largest_gap_deg:.2f} deg)"
    )
    print(f"  {len(result.step_paths)} STEP files, {len(result.plot_paths)} plots")
    print(f"  manifest: {result.manifest_path}")
    print(f"  report:   {result.library_report_path}")
    print(f"  validation: {result.library_validation_path}")
    print("\nNo HFSS comparison has been performed; see reports/LIBRARY_VALIDATION.md.")


if __name__ == "__main__":
    main(sys.argv[1:])
