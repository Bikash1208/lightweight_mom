"""Phase 29 compact 7-tile package (MODEL_ONLY): stub tests of the generated AEDT scripts and the committed package files (no AEDT needed)."""

import ast
import csv
import json
import re
from pathlib import Path

import numpy as np
import pytest

from src.export import compact7_export as ex
from src.hologram import seven_hex_cluster as c

PKG = Path(__file__).resolve().parents[1] / "hfss_exports/phase29_compact7"


def _cells(n=3, r=40.0):
    return [
        {"x_mm": 6.4 * k, "y_mm": 0.0, "slot_L_mm": 4.8, "slot_W_mm": 1.6, "slot_angle_deg": 10.0 * k} for k in range(n)
    ]


def test_board_covers_all_seven_hexagons_with_margin():
    wx, wy = ex.board_size_mm(40.0)
    for q in c.tile_centres_mm(40.0):
        for x, y in ex.tile_hex_vertices(q[:2], 40.0):
            assert abs(x) + ex.MARGIN_MM <= wx / 2 + 1e-9 and abs(y) + ex.MARGIN_MM <= wy / 2 + 1e-9


def test_build_script_is_valid_python_with_no_hardcoded_design_names_and_guards():
    s = ex.aedt_build_script(_cells(), 40.0)
    ast.parse(s)
    assert (
        s.count('"phase29_compact7"') == 1
    )  # the design name appears once (DESIGN_NAME); later calls read oDesign.GetName()
    assert "oDesign.GetName()" in s and "def note(" in s and "except Exception" in s
    assert s.count("oBnd.AssignLumpedPort(") == 1 and "for k, (fx, fy) in enumerate(FEEDS)" in s
    assert re.search(r"FEEDS = \[\((0\.000, 0\.000)\)", s) and s.count("(69.282, 0.000)") == 1
    assert 'SetActiveDesign("phase29' not in s  # no hard-coded design names after creation


def test_weights_script_is_valid_python_and_prints_a_fallback():
    rows = [{"theta_deg": 0.0, "phi_deg": 0.0, "weights": [[0.8, 0.0], [0.0, 0.6]] + [[0.0, 0.0]] * 5}]
    s = ex.aedt_weights_script(rows)
    ast.parse(s)
    assert "EditSources" in s and "except Exception as exc" in s and "(0.64, 0.00)" in s and "(0.36, 90.00)" in s


@pytest.mark.skipif(not PKG.exists(), reason="package not generated")
def test_committed_package_is_consistent():
    pred = json.loads((PKG / "predictions.json").read_text())
    assert "MODEL_ONLY" in pred["labels"] and "INTER_TILE_COUPLING_NOT_MODELLED" in pred["labels"]
    assert len(pred["per_port_alone"]) == 7 and set(pred["weighted_examples"]) == {
        "broadside",
        "theta15_phi0",
        "theta30_phi30",
    }
    assert (
        "NOT_AVAILABLE_FROM_CURRENT_DATA" in pred["s11"] and "NOT_AVAILABLE_FROM_CURRENT_DATA" in pred["port_coupling"]
    )
    rows = list(csv.DictReader((PKG / "weights_table.csv").open(encoding="utf-8")))
    assert len(rows) == 361  # broadside + 15 theta rings x 24 phi
    for r in rows:
        p = np.array([float(r[f"P{k}_power_W"]) for k in range(1, 8)])
        assert p.sum() == pytest.approx(1.0, abs=2e-4) and np.count_nonzero(p > 0) <= 3
    n_cells = len(list(csv.DictReader((PKG / "layout.csv").open(encoding="utf-8"))))
    assert n_cells == pred["n_cells"] and "CELLS = [" in (PKG / "build_board.py").read_text(encoding="utf-8")
    assert (PKG / "build_board.py").read_text(encoding="utf-8").count("\n    (") >= n_cells
    cuts = list(csv.DictReader((PKG / "pattern_cuts.csv").open(encoding="utf-8")))
    assert {r["case"] for r in cuts} == {f"port{k}_alone" for k in range(1, 8)} | {
        "weighted_broadside",
        "weighted_theta15_phi0",
        "weighted_theta30_phi30",
    }
    for f in ("build_board.py", "set_weights.py"):
        ast.parse((PKG / f).read_text(encoding="utf-8"))
