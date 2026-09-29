"""Phase 8A Sec. 32-34: library integrity and cross-file consistency, verified from the actual
delivered Phase 7 outputs (`results/unit_cell_library/`), never taken on the manifest's word
alone -- every manifest field this module reports on is *recomputed* independently from the JSON/
CSV/geometry files and compared, so a manifest that lies about its own build is caught, not
trusted.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from .geometry_export import step_filename
from .library_validation import check_json_csv_consistency


@dataclass(frozen=True)
class LibraryBundle:
    base_dir: Path
    manifest: dict
    library_doc: dict            # unit_cell_library.json, parsed from disk
    json_path: Path
    csv_path: Path
    step_files: frozenset        # file names present in geometry/


def load_library_bundle(library_dir):
    """Loads the four things Sec. 33 asks to be cross-checked: database (JSON), CSV, manifest,
    geometry (STEP file listing). Raises FileNotFoundError, never substitutes a missing file."""
    library_dir = Path(library_dir)
    manifest_path = library_dir / "manifest.json"
    json_path = library_dir / "database" / "unit_cell_library.json"
    csv_path = library_dir / "database" / "unit_cell_library.csv"
    geom_dir = library_dir / "geometry"
    for p in (manifest_path, json_path, csv_path):
        if not p.exists():
            raise FileNotFoundError(f"Phase 8A: expected library file is missing: {p}")
    manifest = json.loads(manifest_path.read_text())
    library_doc = json.loads(json_path.read_text())
    step_files = frozenset(p.name for p in geom_dir.glob("*.step")) if geom_dir.exists() else frozenset()
    return LibraryBundle(library_dir, manifest, library_doc, json_path, csv_path, step_files)


@dataclass(frozen=True)
class LibraryIntegrityResult:
    ok: bool
    errors: tuple
    warnings: tuple
    n_cells: int
    n_ok_cells: int
    n_failed_cells: int
    n_step_files: int
    consistency: object           # library_validation.ConsistencyReport


def check_library_integrity(bundle, expected_n_geometries=324, expected_n_frequency_points=201):
    """Sec. 32-33's checklist, computed from ``bundle`` (already loaded independently from disk
    by `load_library_bundle`). Any mismatch is a `LIBRARY_INTEGRITY FAILURE` (an ``errors`` entry
    and ``ok=False``), never silently ignored."""
    errors, warnings = [], []
    cells = bundle.library_doc.get("cells", [])
    n_cells = len(cells)
    if n_cells != expected_n_geometries:
        errors.append(f"unit_cell_library.json has {n_cells} cell entries, expected {expected_n_geometries}")

    cell_ids = [c["cell_id"] for c in cells]
    dupes = sorted({cid for cid in cell_ids if cell_ids.count(cid) > 1})
    if dupes:
        errors.append(f"duplicate cell_id(s) in unit_cell_library.json: {dupes}")

    ok_cells = [c for c in cells if c.get("solver_status") == "OK"]
    failed_cells = [c for c in cells if c.get("solver_status") != "OK"]

    for c in ok_cells:
        n_f = len(c.get("frequency_GHz", []))
        if n_f != expected_n_frequency_points:
            errors.append(f"cell {c['cell_id']}: {n_f} frequency points, expected {expected_n_frequency_points}")
        for fam in ("S11", "S21", "S12", "S22"):
            n_fam = len(c.get(fam, []))
            if n_fam != n_f:
                errors.append(f"cell {c['cell_id']}: {fam} has {n_fam} entries, frequency_GHz has {n_f}")
        expected_step = step_filename(c["L_mm"], c["W_mm"])
        if expected_step not in bundle.step_files:
            errors.append(f"cell {c['cell_id']}: expected STEP file {expected_step} not found in geometry/")

    n_step_files = len(bundle.step_files)
    if n_step_files != len(ok_cells):
        errors.append(
            f"{n_step_files} STEP file(s) present in geometry/, expected exactly one per "
            f"successful cell ({len(ok_cells)})"
        )

    m = bundle.manifest
    if m.get("number_of_cells") != n_cells:
        errors.append(f"manifest number_of_cells={m.get('number_of_cells')!r} disagrees with actual count {n_cells}")
    if m.get("number_of_successful_cells") != len(ok_cells):
        errors.append(
            f"manifest number_of_successful_cells={m.get('number_of_successful_cells')!r} "
            f"disagrees with actual count {len(ok_cells)}"
        )
    if m.get("number_of_failed_cells") != len(failed_cells):
        errors.append(
            f"manifest number_of_failed_cells={m.get('number_of_failed_cells')!r} disagrees "
            f"with actual count {len(failed_cells)}"
        )
    if m.get("frequency_count") != expected_n_frequency_points:
        errors.append(
            f"manifest frequency_count={m.get('frequency_count')!r} disagrees with expected "
            f"{expected_n_frequency_points}"
        )

    consistency = check_json_csv_consistency(bundle.json_path, bundle.csv_path)
    if not consistency.ok:
        errors.append(f"unit_cell_library.json/.csv are inconsistent: {len(consistency.mismatches)} mismatch(es)")

    return LibraryIntegrityResult(
        ok=len(errors) == 0, errors=tuple(errors), warnings=tuple(warnings),
        n_cells=n_cells, n_ok_cells=len(ok_cells), n_failed_cells=len(failed_cells),
        n_step_files=n_step_files, consistency=consistency,
    )
