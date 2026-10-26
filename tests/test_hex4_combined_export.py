"""Phase 26: `hfss_exports/hex4_tensor_cp_combined/` (board B, combined beam) reproduces the committed design `results/phase26/fourhex/fourhex_B_cells.json`."""

import ast
import contextlib
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "hfss_exports" / "hex4_tensor_cp_combined"
DESIGN = ROOT / "results" / "phase26" / "fourhex" / "fourhex_B_cells.json"


def _script_assignments():
    tree = ast.parse((PKG / "hfss_build_board.py").read_text())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            with contextlib.suppress(ValueError):
                out[node.targets[0].id] = ast.literal_eval(node.value)
    return out


def test_script_cells_match_the_design_json():
    cells = json.loads(DESIGN.read_text())["cells"]
    script = np.array(_script_assignments()["CELLS"], dtype=float)
    ref = np.array([[c["x_mm"], c["y_mm"], c["slot_L_mm"], c["slot_W_mm"], c["slot_angle_deg"]] for c in cells])
    assert script.shape == ref.shape == (5684, 5)
    np.testing.assert_allclose(script[:, :4], ref[:, :4], atol=1e-3 + 1e-9)  # printed to 3 decimals
    assert np.abs(script[:, 4] - ref[:, 4]).max() < 1e-4 + 1e-9  # printed to 4 decimals


def test_layout_csv_matches_the_design_json_and_ports():
    cells = json.loads(DESIGN.read_text())["cells"]
    with open(PKG / "layout.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == len(cells)
    for r, c in zip(rows, cells, strict=True):
        assert abs(float(r["x_mm"]) - c["x_mm"]) < 1e-3 + 1e-9
        assert abs(float(r["y_mm"]) - c["y_mm"]) < 1e-3 + 1e-9
        assert abs(float(r["slot_angle_deg"]) - c["slot_angle_deg"]) < 1e-3 + 1e-9
        assert int(r["port"]) == c["port"]


def test_combined_beam_drive_and_gitignored_step():
    assert _script_assignments()["DRIVE_PORT"] == 5  # all four ports 0.25 W / 0 deg = combined beam
    gi = (ROOT / ".gitignore").read_text()
    assert "hfss_exports/**/*.step" in gi  # board.step (27.9 MB) stays out of git
