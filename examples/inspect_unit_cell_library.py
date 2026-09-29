"""Query an already-built Phase 7 unit-cell library (req. 24-25): exact-geometry lookup, an
exact-frequency response, phase-candidate search and phase inversion. Read-only; never re-solves.

Run from the project root, after examples/build_unit_cell_library.py has produced the library:
    python examples/inspect_unit_cell_library.py [results/unit_cell_library] [target_phase_deg]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.unit_cell_library import UnitCellLibrary

ROOT = Path(__file__).resolve().parent.parent


def main(argv):
    base_dir = Path(argv[0]) if argv else ROOT / "results" / "unit_cell_library"
    target_phase_deg = float(argv[1]) if len(argv) > 1 else 45.0

    json_path = base_dir / "database" / "unit_cell_library.json"
    library = UnitCellLibrary.from_json(json_path)
    n_cells = len(library._doc["cells"])
    print(f"Loaded {json_path} ({n_cells} geometries)")

    manifest_path = base_dir / "manifest.json"
    if manifest_path.exists():
        import json

        manifest = json.loads(manifest_path.read_text())
        print(
            f"design frequency: {manifest['design_frequency_GHz']} GHz, "
            f"phase_source: {manifest['phase_source']} ({manifest['phase_source_status']})"
        )
        f_design = manifest["design_frequency_GHz"]
    else:
        f_design = 12.5

    example_l, example_w = 2.4, 2.4
    cell = library.get_cell(example_l, example_w)
    if cell is not None:
        print(
            f"\nget_cell(L={example_l}, W={example_w}): {cell['cell_id']} "
            f"geometry_status={cell['geometry_status']} solver_status={cell['solver_status']}"
        )
        try:
            resp = library.response(example_l, example_w, f_design)
            print(
                f"response(L={example_l}, W={example_w}, f={f_design} GHz): "
                f"|S11xx|={resp['mag']:.4f} ({resp['mag_db']:.2f} dB), phase={resp['phase_deg']:.2f} deg"
            )
        except (KeyError, ValueError) as e:
            print(f"response(...) unavailable: {e}")

    print(f"\nfind_phase_candidates(target={target_phase_deg} deg, f={f_design} GHz):")
    candidates = library.find_phase_candidates(target_phase_deg, f_design)
    if not candidates:
        print("  no candidates within the default tolerance")
    for c in candidates[:10]:
        print(
            f"  L={c.l_mm:g} W={c.w_mm:g} mm  phase={c.phase_deg:.2f} deg  "
            f"error={c.phase_error_deg:+.2f} deg  |S11xx|={c.magnitude_db:.2f} dB"
        )

    inversion = library.invert_phase(target_phase_deg, f_design)
    print(f"\ninvert_phase(target={target_phase_deg} deg, f={f_design} GHz): {inversion}")


if __name__ == "__main__":
    main(sys.argv[1:])
