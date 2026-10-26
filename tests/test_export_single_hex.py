"""Phase 26 D: slot polygon corners of the single-hexagon export are exact rotations of the rectangle."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

_p = Path(__file__).resolve().parents[1] / "scripts/research/phase26_export_single_hex.py"
_s = importlib.util.spec_from_file_location("p26_export", _p)
_m = importlib.util.module_from_spec(_s)
sys.modules["p26_export"] = _m
_s.loader.exec_module(_m)


def test_slot_corners_are_the_rotated_rectangle():
    pts = np.array(_m.slot_corners(10.0, -5.0, 4.8, 2.0, 30.0))
    c = pts.mean(axis=0)
    assert np.allclose(c, [10.0, -5.0], atol=1e-3)
    sides = np.linalg.norm(np.diff(np.vstack([pts, pts[:1]]), axis=0), axis=1)
    assert np.allclose(sides, [4.8, 2.0, 4.8, 2.0], atol=1e-3)
    ax = pts[1] - pts[0]
    np.testing.assert_allclose(np.degrees(np.arctan2(ax[1], ax[0])), 30.0, atol=0.01)


def test_generated_script_is_valid_python_with_every_cell():
    cells = [{"x_mm": 0.0, "y_mm": 6.4, "slot_L_mm": 4.8, "slot_W_mm": 2.0, "slot_angle_deg": 12.5}] * 3
    src = _m.aedt_script(cells)
    compile(src, "hfss_build_board.py", "exec")
    assert src.count("(0.000, 6.400, 4.800, 2.000, 12.5000, 0)") == 3
