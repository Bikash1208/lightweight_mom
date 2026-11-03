"""Phase 30: refined-screen loss estimates and test-board design invariants (scripts/research/phase30_*.py)."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase30_refine as R  # noqa: E402
import phase30_testboard as T  # noqa: E402
from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.sheet_model import find_beta_from_sheet_admittance  # noqa: E402

MESH = 24
SPEC = dict(major=3.05e-3, ratio=1.0, slot_l=1.9e-3, slot_w=0.5e-3, angle=0.0)


@pytest.fixture(scope="module")
def solved():
    mesh = MeshSettings(MESH, MESH, MESH)
    cell = R.make_cell(**SPEC)
    basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    yn = ts.sheet_admittance_tensor_for_cell(cell, R.F0, mesh)
    beta0 = find_beta_from_sheet_admittance(1j * 0.5 * (yn[0, 0].imag + yn[1, 1].imag), R.SUB, R.F0).beta_rad_per_m
    return cell, mesh, basis, beta0, 0.5 * (yn[0, 0].imag + yn[1, 1].imag)


def test_dielectric_alpha_is_linear_in_loss_tangent(solved):
    _, _, _, beta0, y_bar = solved
    a1 = R.dielectric_alpha(y_bar, beta0, R.F0)
    old = R.SUB
    try:
        R.SUB = type(old)(old.name, old.eps_r, 2.0 * R.TAND, old.thickness)
        a2 = R.dielectric_alpha(y_bar, beta0, R.F0)
    finally:
        R.SUB = old
    assert a1 > 0.0
    assert a2 / a1 == pytest.approx(2.0, rel=0.03)


def test_conductor_alpha_scales_as_inverse_sqrt_conductivity_and_splits(solved):
    cell, mesh, basis, beta0, _ = solved
    a, a_patch, a_gnd = R.conductor_alpha(cell, R.F0, mesh, beta0, basis)
    assert a == pytest.approx(a_patch + a_gnd, rel=1e-12)
    assert a_patch > 0.0 and a_gnd > 0.0
    old = R.SIGMA_CU
    try:
        R.SIGMA_CU = 4.0 * old
        b, _, _ = R.conductor_alpha(cell, R.F0, mesh, beta0, basis)
    finally:
        R.SIGMA_CU = old
    assert b / a == pytest.approx(0.5, rel=1e-9)


def test_grid_has_no_unphysical_slot_lengths():
    """Slot ends keep >= 0.3 mm of metal bridge (a bridge below one mesh cell is a split patch, not a slot)."""
    for s in R.grid_specs():
        if s["slot_l"] > 0:
            assert (s["major"] - s["slot_l"]) / 2 >= 0.3e-3


def test_testboard_lattice_is_fourfold_symmetric_and_clear_of_the_feed():
    x, y = T.lattice(T.R_BOARD)
    assert np.hypot(x, y).min() >= T.RHO_CLEAR
    assert np.hypot(x, y).max() <= T.R_BOARD
    pts = {(round(a, 9), round(b, 9)) for a, b in zip(x, y, strict=True)}
    assert all(((round(-a, 9), round(b, 9)) in pts) and ((round(a, 9), round(-b, 9)) in pts) for a, b in pts)


def test_eigenmode_step_geometry_roundtrip(tmp_path):
    """The STEP of E1 re-reads with 4 named bodies, the 0.15 mm gap and the L x W slot (cadquery/OCP)."""
    pytest.importorskip("cadquery")
    import phase30_eigenmode_step as S

    path = tmp_path / "E1.step"
    S.build("E1_strong_M0p64", path)
    out = S.check("E1_strong_M0p64", path)
    assert out["all_pass"], out["checks"]
    assert out["slot_LxW_mm_from_area_and_bbox"] == [1.89, 0.58]
