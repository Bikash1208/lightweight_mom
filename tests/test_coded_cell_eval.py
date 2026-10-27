"""Phase 28 item 3: the tensor-sheet evaluator of diode-switched cells (`coded_cell_eval`) and the candidate geometries (`coded_cell_candidates`)."""

import math
import warnings

import numpy as np
import pytest

from src.dispersion.leakage_model import leakage_alpha_first_order_tensor, state_surface_reactances
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.dispersion.slotted_cells import SlottedCircularPatch
from src.dispersion.tensor_sheet import rotate_tensor
from src.geometry import UnitCell
from src.hologram import coded_cell_candidates as cc
from src.hologram import coded_cell_eval as ev
from src.materials import Metal
from src.periodic_bc import RectangularLattice
from src.ris.cells import OFF, ON, RT5880, ideal_kit, madp_kit

F10, F1054 = 10.0e9, 10.54e9


def _slotted(d, length, width, angle_deg):
    geo = SlottedCircularPatch(d * 1e-3, length * 1e-3, width * 1e-3, math.radians(angle_deg))
    return UnitCell(RectangularLattice(6.4e-3, 6.4e-3), geo, RT5880, Metal(), backing="pec")


# ------------------------------------------------------------------------------------------------ the metrics on a known tensor
def test_principal_metrics_on_a_synthetic_lossless_tensor():
    b, d = 7.0e-3, 2.0e-3
    y = 1j * np.array([[b + d, 0.0], [0.0, b - d]])
    m = ev.principal_metrics(y, RT5880, F10)
    assert m["b_bar_ms"] == pytest.approx(7.0) and m["delta_ms"] == pytest.approx(2.0)
    assert m["g_max_ms"] == pytest.approx(0.0, abs=1e-12) and m["asymmetry_ms"] == pytest.approx(0.0, abs=1e-12)
    beta = find_beta_from_sheet_admittance(1j * b, RT5880, F10).beta_rad_per_m
    xs = state_surface_reactances([1j * (b + d), 1j * (b - d)], beta, RT5880, F10)
    assert m["dx_aniso_ohm"] == pytest.approx(abs(xs[0] - xs[1]), rel=1e-12)
    alpha = leakage_alpha_first_order_tensor(beta, 1j * b, d, math.radians(30.0), RT5880, F10)
    assert m["alpha_beam_np_per_m"] == pytest.approx(alpha, rel=1e-12)
    # the design-loading contrast is the same anisotropy on the beta = 1.06 mean loading
    b_des = ev.design_loading(F10, RT5880, 1.06)
    beta_des = find_beta_from_sheet_admittance(1j * b_des, RT5880, F10).beta_rad_per_m
    assert beta_des / (2 * math.pi * F10 / 299792458.0) == pytest.approx(1.06, abs=1e-6)


def test_rotating_the_tensor_rotates_the_axis_and_keeps_delta_and_the_contrast():
    y0 = 1j * np.array([[9.0e-3, 0.0], [0.0, 5.0e-3]])
    m0 = ev.principal_metrics(y0, RT5880, F10)
    m1 = ev.principal_metrics(rotate_tensor(y0, math.radians(45.0)), RT5880, F10)
    for k in ("b_bar_ms", "delta_ms", "dx_aniso_ohm", "alpha_beam_np_per_m"):
        assert m1[k] == pytest.approx(m0[k], rel=1e-12)
    assert (m1["axis_deg"] - m0["axis_deg"]) % 180.0 == pytest.approx(45.0, abs=1e-9)


# ------------------------------------------------------------------------------------------------ the extraction against Phase 26 (calibration)
def test_evaluator_reproduces_the_phase26_slotted_cell_scale_and_the_90_degree_equivalence():
    """Plain slotted circular patch (Phase 26 B2 family) at 10.54 GHz, mesh 32: b_bar and delta in the range of the committed supercell table
    (delta 2.3-2.5 mS, ~100 ohm); a 90-degree rotated cell is the same tensor rotated (exact on the square mesh)."""
    y0 = ev.sheet_tensor(_slotted(6.0, 4.8, 1.6, 0.0), F1054, 32)
    m0 = ev.principal_metrics(y0, RT5880, F1054)
    assert m0["b_bar_ms"] == pytest.approx(7.953, abs=0.01) and m0["delta_ms"] == pytest.approx(2.487, abs=0.01)
    assert m0["dx_aniso_ohm"] == pytest.approx(110.3, abs=0.5)
    y90 = ev.sheet_tensor(_slotted(6.0, 4.8, 1.6, 90.0), F1054, 32)
    assert np.linalg.norm(y90 - rotate_tensor(y0, math.pi / 2)) / np.linalg.norm(y0) < 1e-9


# ------------------------------------------------------------------------------------------------ geometries
def test_sector_hub_geometry_is_symmetric_clipped_and_the_necks_are_metal():
    cell = cc.SectorHubCell(
        period=6.4e-3,
        radius=3.0e-3,
        radius_diag=3.4e-3,
        hub_radius=0.8e-3,
        radial_gap=0.8e-3,
        neck_width=0.4e-3,
        sector_gap=0.4e-3,
    )
    geo = cell.geometry()
    g = (np.arange(128) + 0.5) * 0.05e-3 - 3.2e-3  # symmetric grid, no point on a boundary
    x, y = np.meshgrid(g, g, indexing="ij")
    m = geo.in_feature(x, y)
    assert np.array_equal(m, geo.in_feature(-y, x))  # 90-degree rotation (the grid is symmetric)
    assert np.array_equal(m, geo.in_feature(-x, y))  # mirror
    assert not np.any(geo.in_feature(np.array([3.1e-3, -3.15e-3]), np.array([0.0, 0.0])))  # clipped inside the cell
    for (sx, sy, length, width, _direction), a in zip(cell.diode_sites(), cell.arms_deg(), strict=True):
        ang = math.radians(a)
        u = np.linspace(-length / 2, length / 2, 5) * 0.999
        v = np.linspace(-width / 2, width / 2, 5) * 0.999
        uu, vv = np.meshgrid(u, v, indexing="ij")
        px = sx + uu * math.cos(ang) - vv * math.sin(ang)
        py = sy + uu * math.sin(ang) + vv * math.cos(ang)
        assert np.all(geo.in_feature(px, py)), a


def test_sector_hub_states_are_one_hot_pairs_on_the_axis():
    cell = cc.SectorHubCell()
    sts = cell.states()
    assert len(sts) == 4 and all(len(s) == 8 for s in sts)
    for k, s in enumerate(sts):
        on = [a for a, v in zip(cell.arms_deg(), s, strict=True) if v == ON]
        assert sorted(on) == sorted([45.0 * k, 45.0 * k + 180.0])
    assert all(sum(v == ON for v in s) == 1 for s in cc.SectorHubCell(single_sided=True).states())
    assert cell.n_diodes == 8 and all(v in (ON, OFF) for s in sts for v in s)


def test_star_cell_states_shorten_the_perpendicular_slot():
    cell = cc.ShortenableStarCell()
    assert cell.n_diodes == 8
    for k, s in enumerate(cell.states()):
        on = {a % 180.0 for a, v in zip(cell.arms_deg(), s, strict=True) if v == ON}
        assert on == {(cc.ANGLES_DEG[k] + 90.0) % 180.0}


# ------------------------------------------------------------------------------------------------ the cell physics (ideal switch, mesh 32)
@pytest.fixture(scope="module")
def ideal_states():
    cell = cc.SectorHubCell(
        period=6.4e-3,
        radius=3.0e-3,
        radius_diag=3.4e-3,
        hub_radius=0.8e-3,
        radial_gap=0.8e-3,
        neck_width=0.4e-3,
        sector_gap=0.4e-3,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ys = [ev.sheet_tensor(cell.unit_cell(st, ideal_kit()), F10, 32) for st in cell.states()]
    return cell, ys


def test_opposite_axis_states_are_exact_90_degree_rotations_and_reciprocal(ideal_states):
    _, ys = ideal_states
    assert np.linalg.norm(ys[2] - rotate_tensor(ys[0], math.pi / 2)) / np.linalg.norm(ys[0]) < 1e-6
    assert np.linalg.norm(ys[3] - rotate_tensor(ys[1], math.pi / 2)) / np.linalg.norm(ys[1]) < 1e-6
    for y in ys:
        assert abs(y[0, 1] - y[1, 0]) / max(abs(y[0, 0]), abs(y[1, 1])) < 1e-6


def test_the_on_axis_sets_the_principal_axis_and_a_real_switch_is_not_the_ideal_one(ideal_states):
    cell, ys = ideal_states
    m = [ev.principal_metrics(y, RT5880, F10) for y in ys]
    assert abs(((m[0]["axis_deg"] - m[2]["axis_deg"]) % 180.0) - 90.0) < 1.0
    assert all(x["delta_ms"] > 0.5 for x in m)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y_madp = ev.sheet_tensor(cell.unit_cell(cell.states()[0], madp_kit()), F10, 32)
    assert not np.allclose(y_madp, ys[0], rtol=1e-3)  # the diode models change the tensor (the loads act)


# ------------------------------------------------------------------------------------------------ the rational-slope and Manhattan hub cells
def _grid(n=160, half=4.0e-3):
    g = (np.arange(n) + 0.5) * (2 * half / n) - half  # symmetric, no point on a multiple of 0.05 mm
    return np.meshgrid(g, g, indexing="ij")


@pytest.mark.parametrize(
    "cell",
    [
        cc.RationalSectorCell(
            period=8.0e-3,
            hub_half=1.6e-3,
            neck_gap=0.4e-3,
            neck_width=0.8e-3,
            gap_t=0.4e-3,
            arm_outer=3.6e-3,
            diag_outer=3.2e-3,
            clip=3.8e-3,
        ),
        cc.RationalSectorCell(
            period=6.4e-3,
            hub_half=1.2e-3,
            neck_gap=0.4e-3,
            neck_width=0.4e-3,
            gap_t=0.2e-3,
            arm_outer=3.0e-3,
            diag_outer=2.8e-3,
            clip=3.0e-3,
        ),
    ],
)
def test_rational_sector_geometry_symmetry_clip_and_loads_on_metal(cell):
    geo = cell.geometry()
    x, y = _grid(160, cell.period / 2)
    m = geo.in_feature(x, y)
    assert np.array_equal(m, geo.in_feature(-y, x))  # C4: the pinwheel is invariant under 90 degrees
    assert np.array_equal(m, geo.in_feature(-x, -y))
    assert not np.any(
        geo.in_feature(np.array([cell.period / 2 - 1e-5]), np.array([0.0]))
    )  # inside the cell with a margin
    assert m.mean() > 0.3  # a substantial metal fill
    for (sx, sy, length, width, d), a in zip(cell.diode_sites(), range(8), strict=True):
        u = np.linspace(-length / 2, length / 2, 5) * 0.999
        v = np.linspace(-width / 2, width / 2, 5) * 0.999
        uu, vv = np.meshgrid(u, v, indexing="ij")
        px, py = (sx + uu, sy + vv) if d == "x" else (sx + vv, sy + uu)
        assert np.all(geo.in_feature(px, py)), a
    sts = cell.states()
    assert len(sts) == 4 and all(sum(v == ON for v in s) == 2 for s in sts)
    notch = cc.RationalSectorCell(**{**cell.__dict__, "mode": "notch"})
    assert all(sum(v == ON for v in s) == 6 for s in notch.states())


def test_rational_sector_necks_have_the_right_conduction_direction_and_no_load_warning():
    cell = cc.RationalSectorCell(
        period=8.0e-3,
        hub_half=1.6e-3,
        neck_gap=0.4e-3,
        neck_width=0.8e-3,
        gap_t=0.4e-3,
        arm_outer=3.6e-3,
        diag_outer=3.2e-3,
        clip=3.8e-3,
    )
    dirs = [d for *_, d in cell.diode_sites()]
    assert (
        dirs == ["x", "y", "x", "y"] * 2
    )  # the base necks are horizontal and each 90-degree rotation swaps the direction
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # a load whose region is not metal raises a UserWarning in check_load_geometry
        cell.unit_cell(cell.states()[0], ideal_kit())
        from src.loads import check_load_geometry

        uc = cell.unit_cell(cell.states()[0], ideal_kit())
        check_load_geometry(uc.loads, uc, 40, 40)


def test_manhattan_hub_geometry_is_c4_symmetric_aligned_and_has_all_necks_on_metal():
    cell = cc.ManhattanHubCell()
    geo = cell.geometry()
    x, y = _grid(128, cell.period / 2)
    m = geo.in_feature(x, y)
    assert np.array_equal(m, geo.in_feature(-y, x))
    assert [d for *_, d in cell.diode_sites()] == ["x", "y", "x", "y"] * 2
    for (sx, sy, length, _width, d), _ in zip(cell.diode_sites(), range(8), strict=True):
        px, py = (
            (np.array([sx - length / 2 * 0.99, sx + length / 2 * 0.99]), np.array([sy, sy]))
            if d == "x"
            else (np.array([sx, sx]), np.array([sy - length / 2 * 0.99, sy + length / 2 * 0.99]))
        )
        assert np.all(geo.in_feature(px, py))
    assert all(sum(v == ON for v in s) == 2 for s in cell.states())
    assert all(sum(v == ON for v in s) == 6 for s in cc.ManhattanHubCell(mode="notch").states())


def test_rational_sector_arm_state_is_an_exact_90_degree_rotation_of_the_opposite_arm_state():
    cell = cc.RationalSectorCell(
        period=6.4e-3,
        hub_half=1.2e-3,
        neck_gap=0.4e-3,
        neck_width=0.4e-3,
        gap_t=0.2e-3,
        arm_outer=3.0e-3,
        diag_outer=2.8e-3,
        clip=3.0e-3,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y0 = ev.sheet_tensor(cell.unit_cell(cell.states()[0], ideal_kit()), F10, 32)
        y2 = ev.sheet_tensor(cell.unit_cell(cell.states()[2], ideal_kit()), F10, 32)
    assert np.linalg.norm(y2 - rotate_tensor(y0, math.pi / 2)) / np.linalg.norm(y0) < 1e-6


def test_ideal_switch_loads_agree_with_a_plain_pattern_without_loads():
    """The lumped-diode treatment in the candidate cells: an ideal-switch cell (ON necks = loads of 1 mohm, OFF necks = 1e-18 F) has the same sheet tensor as the same cell
    built with the ON necks as plain metal and the OFF necks removed from the pattern (no loads at all)."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
    import phase28_neck_load_check as chk

    r = chk.compare_state(0, 24)
    assert r["n_on"] == 6 and r["n_off"] == 2
    assert r["relative_tensor_difference"] < 5e-3
    assert abs(r["alpha_with_loads"] - r["alpha_plain_pattern"]) / r["alpha_plain_pattern"] < 5e-2
