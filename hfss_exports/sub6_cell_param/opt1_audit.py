"""Final independent row, build, VERIFY, and readback audit for Option 1."""
import json
from pathlib import Path

import cell_param_pyaedt as cp
from opt1_cell import configure, values
from opt1_geometry import Geometry, all_geometries
from opt1_queue import cases

configure()
here = Path(__file__).resolve().parent
root = here.parents[1] / "results" / "sub6_cell_param" / "lab_B_opt1"
geometries = [Geometry(**item["geometry"]) for item in all_geometries(0.5)]
expected = {row["*"]: (g, row) for g in geometries for row in cases(g)}
found = {}
for file in (root / "rows").glob("*.json"):
    row = json.loads(file.read_text())
    run_id = row["run_id"]
    if run_id != file.stem or run_id not in expected or run_id in found:
        raise RuntimeError(f"Unexpected/duplicate row {run_id}")
    g, source = expected[run_id]
    intended = values(g, source["state"], source["direction"], int(source["phase_deg"]))
    good, why = cp.readback_ok(row)
    if row["project_variant"] != g.name or row["variables"] != intended or len(row["modes"]) != 8 or not good:
        raise RuntimeError(f"Invalid row {run_id}: {why}")
    found[run_id] = row
if set(found) != set(expected):
    raise RuntimeError(f"Missing {len(set(expected)-set(found))} rows")
failed = json.loads((root / "FAILED_ROWS.json").read_text())
if failed:
    raise RuntimeError(f"FAILED rows: {len(failed)}")
for g in geometries:
    build = json.loads((root / "builds" / f"{g.name}.json").read_text())
    stack = build["stack"]
    if build["status"] != "BUILT" or build["aedt"] != "2025.2" or build["mesh_operations"]:
        raise RuntimeError(f"Bad build {g.name}")
    if not stack["solid_ground"] or stack["fr4"] or stack["vias"] or stack["antipads"] or stack["stubs"] or stack["air_below_ground"]:
        raise RuntimeError(f"Wrong stack {g.name}")
    if not {"DISC", "ISLAND"}.issubset(build["objects"]) or build["geometry"]["checks"]["island_local_disc_copper_dc_contact"]:
        raise RuntimeError(f"Local DC check failed {g.name}")
    if build["geometry"]["geometry"]["line_w"] != g.line_w or build["geometry"]["geometry"]["channel_w"] != g.channel_w:
        raise RuntimeError(f"Line/channel width mismatch {g.name}")
verify = json.loads((root / "verify" / "determinism.json").read_text())
if verify["status"] != "PASS" or verify["external_reference_used"]:
    raise RuntimeError("Determinism VERIFY failed")
audit = {"status": "PASS", "aedt": "2025.2", "geometries": len(geometries),
         "expected_rows": len(expected), "solved_rows": len(found), "failed_rows": 0,
         "modes_per_row": 8, "readback_validated_rows": len(found),
         "geometry_dc_passed": True, "solid_ground_no_via_no_fr4": True,
         "determinism_verify": "PASS", "external_reference_used": False}
(root / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
print(f"PASS: {len(found)}/{len(expected)} rows, readback and build/DC gates, 0 FAILED")
