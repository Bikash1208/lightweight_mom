"""Phase 26: 4-hexagon boards -- shared lattice, edge-power bookkeeping, export script, global phase reference."""

import ast
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
import phase26_export_fourhex as ex  # noqa: E402
import phase26_fourhex_design as fd  # noqa: E402
import phase26_tensor_cp_hexagon as cp  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable  # noqa: E402


@pytest.fixture(scope="module")
def boards():
    table = SupercellTable.from_files(cp.TABLE_GLOB)
    return {b: fd.design_board(b, table) for b in ("A", "B")}


def test_shared_lattice_cells_are_distinct_sites_and_at_least_one_pitch_apart(boards):
    xy = np.concatenate([np.c_[boards["A"][j]["x_mm"], boards["A"][j]["y_mm"]] for j in fd.FEEDS])
    from scipy.spatial import cKDTree

    d, _ = cKDTree(xy).query(xy, k=2)
    assert d[:, 1].min() >= 6.4 - 1e-6  # no overlapping D 6.0 circles, no duplicated site


def test_unconverted_power_is_fully_assigned_to_free_or_shared_edges(boards):
    edges = fd.crossing_edges(boards["A"])
    for row in edges.values():
        parts = sum(v for k, v in row.items() if k != "total_unconverted")
        assert parts == pytest.approx(row["total_unconverted"], rel=1e-9)
        assert 0.3 < row["total_unconverted"] < 0.6  # eta_conv ~ 0.54


def test_global_phase_reference_only_rotates_slots_of_a_hexagon_with_same_beam(boards):
    """Port 3 aims at phi 225 in both boards: B differs from A by a uniform slot rotation -1/2 kt u.r_feed (mod 180 deg)."""
    a, b = boards["A"][3]["design"], boards["B"][3]["design"]
    da = (np.degrees(b["slot_angle_rad"] - a["slot_angle_rad"]) + 90.0) % 180.0 - 90.0
    kt = cp.K0 * np.sin(np.radians(30.0))
    f = fd.FEEDS[3]
    expect = (np.degrees(-0.5 * kt * (f[0] * 1e-3 * np.cos(np.radians(225)) + f[1] * 1e-3 * np.sin(np.radians(225)))) + 90.0) % 180.0 - 90.0
    assert np.allclose(da, expect, atol=1e-6)


def test_export_script_is_valid_python_with_four_ports_and_drive_option(tmp_path):
    cells = [{"x_mm": 10.0, "y_mm": 12.8, "slot_L_mm": 4.8, "slot_W_mm": 2.0, "slot_angle_deg": 33.0, "port": 1}]
    src = ex.aedt_script(cells, beams="test")
    ast.parse(src)
    assert "DRIVE_PORT = 0" in src and "def set_sources" in src and src.count("Port%d") >= 2 and "FEEDS = [(64.952, 112.500)" in src
    assert "(\"mono_h\", \"7.575mm\")" in src
