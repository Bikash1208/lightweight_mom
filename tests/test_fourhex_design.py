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
    src = ex.aedt_script_fourhex(cells, beams="test")
    ast.parse(src)
    assert "DRIVE_PORT = 0" in src and "def set_sources" in src and src.count("Port%d") >= 2 and "FEEDS = [(64.952, 112.500)" in src
    assert "(\"mono_h\", \"7.575mm\")" in src


def test_march_density_fn_default_is_bit_identical_and_factor_acts(boards):
    from src.hologram.tensor_leaky_cp import march_rays

    table = SupercellTable.from_files(cp.TABLE_GLOB)
    a = march_rays(table, cp.hexagon(), (0.0, 0.0), cp.K0, np.radians(30), np.radians(225), 0.9, taper_c=0.5, n_psi=90)
    b = march_rays(table, cp.hexagon(), (0.0, 0.0), cp.K0, np.radians(30), np.radians(225), 0.9, taper_c=0.5, n_psi=90, density_fn=None)
    assert np.array_equal(a["w"], b["w"]) and np.array_equal(a["slot"], b["slot"])
    c = march_rays(table, cp.hexagon(), (0.0, 0.0), cp.K0, np.radians(30), np.radians(225), 0.9, taper_c=0.5, n_psi=90,
                   density_fn=lambda psi, rho, r_edge: np.broadcast_to(1.0 + 4.0 * rho[None, :] / rho.max(), (psi.size, rho.size)))
    assert not np.array_equal(a["w"], c["w"])  # the target-density shape is honoured (a constant is renormalised away)
    k = np.full((psi_n := 90, a["rho"].size), 0.5)
    d = march_rays(table, cp.hexagon(), (0.0, 0.0), cp.K0, np.radians(30), np.radians(225), 0.9, taper_c=0.5, n_psi=psi_n, density_fn=lambda *_: k)
    assert np.array_equal(a["w"], d["w"])


def test_compare_fourhex_on_synthetic_exports(tmp_path):
    import importlib.util

    p = ROOT / "scripts/research/phase26_compare_fourhex.py"
    s = importlib.util.spec_from_file_location("p26_cmp4", p)
    m = importlib.util.module_from_spec(s)
    sys.modules["p26_cmp4"] = m
    s.loader.exec_module(m)
    th, ph = np.arange(0.0, 180.0 + 1e-9, 2.0), np.arange(-180.0, 180.0 + 1e-9, 2.0)
    lines = ['"Freq [GHz]","Phi [deg]","Theta [deg]","dB(RealizedGainLHCP) []","dB(RealizedGainRHCP) []"']
    for q in ph:
        for t in th:
            d = m.single._angle_to(t, np.mod(q, 360.0), 30.0, 225.0)
            r = 23.0 - 12.0 * (d / 8.0) ** 2
            lh = 23.0 - 25.0 - 12.0 * (d / 8.0) ** 2
            dl = m.single._angle_to(t, np.mod(q, 360.0), 30.0, 255.0)
            lh = max(lh, 23.0 - 21.0 - 12.0 * (dl / 5.0) ** 2)
            lines.append(f"10.54,{q:g},{t:g},{lh:.3f},{max(r, -30.0):.3f}")
    (tmp_path / "g3.csv").write_text("\n".join(lines) + "\n")
    # 4-port Touchstone with S_ij = -55 dB everywhere off-diagonal, S_ii = -15 dB
    rows = ["# GHZ S DB R 50"]
    for f in (10.0, 10.54, 11.0):
        vals = []
        for i in range(4):
            for j in range(4):
                vals += [(-15.0 if i == j else -55.0), 0.0]
        rows.append(f"{f} " + " ".join(f"{v}" for v in vals[:8]))
        for k in range(1, 4):
            rows.append(" ".join(f"{v}" for v in vals[8 * k:8 * k + 8]))
    (tmp_path / "b.s4p").write_text("\n".join(rows) + "\n")
    rep = m.main(["--s4p", str(tmp_path / "b.s4p"), "--gain", f"3={tmp_path / 'g3.csv'}", "--pred", str(ROOT / "hfss_exports/hex4_tensor_cp_static/predictions.json")])
    c = rep["ports"]["3"]["checks"]
    assert c["peak_gain"]["state"] == "GAIN_WITHIN_PREDICTED_RANGE"
    assert c["neighbour_lhcp_lobe"]["hfss_at_predicted_direction_db"] == pytest.approx(-21.0, abs=1.0)
    assert rep["coupling"]["pairs"]["1-2"]["s_ij_db"] == pytest.approx(-55.0, abs=0.1)
    assert rep["coupling"]["pairs"]["1-2"]["state"] == "COUPLING_BELOW_RAW_ESTIMATE"
