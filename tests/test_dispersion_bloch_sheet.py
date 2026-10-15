"""Phase 26 B2: slotted cells, Bloch (kt > k0) sheet tensor, Bloch dispersion, supercell rows."""

import numpy as np
import pytest

from src.basis_functions import MeshSettings, RooftopBasis
from src.dispersion import tensor_sheet as ts
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.dispersion.slotted_cells import CellRow, SlottedCircularPatch, SlottedSquarePatch
from src.frequency import free_space_wavenumber
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.materials import Metal, Substrate
from src.mom_matrix import assemble_impedance_matrix
from src.periodic_bc import RectangularLattice

F = 10.54e9
P = 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
K0 = free_space_wavenumber(F)
MESH = MeshSettings(24, 24, 24)


def _cell(geo):
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")


def test_circle_mask_area_and_rotation_invariance_of_the_outline():
    c = SlottedCircularPatch(5.6e-3, 0.0, 0.0)
    x = np.linspace(-3.2e-3, 3.2e-3, 641)
    xx, yy = np.meshgrid(x, x)
    area = c.in_feature(xx, yy).sum() * (x[1] - x[0]) ** 2
    assert area == pytest.approx(np.pi * 2.8e-3**2, rel=2e-3)
    slotted = SlottedCircularPatch(5.6e-3, 4.0e-3, 1.0e-3, 0.7)
    a1 = slotted.in_feature(xx, yy).sum() * (x[1] - x[0]) ** 2
    assert a1 == pytest.approx(np.pi * 2.8e-3**2 - 4.0e-3 * 1.0e-3, rel=5e-3)
    rot = SlottedCircularPatch(5.6e-3, 4.0e-3, 1.0e-3, 0.7 + np.pi / 2)
    assert (rot.in_feature(xx, yy) & ~slotted.in_feature(xx, yy)).sum() > 0  # the slot really turns
    assert slotted.feature_area == pytest.approx(a1, rel=5e-3)


def test_bloch_tensor_at_zero_kt_is_the_normal_incidence_extraction():
    cell = _cell(Geometry("rectangle", 4.0e-3, 2.4e-3))
    a = ts.sheet_admittance_tensor_for_cell(cell, F, MESH)
    b = ts.bloch_sheet_tensor(cell, F, MESH, 0.0, 0.0)
    assert np.allclose(a, b, rtol=1e-8, atol=1e-12)


def test_bare_circle_is_direction_independent_but_longitudinal_differs_from_transverse():
    cell = _cell(SlottedCircularPatch(5.0e-3, 0.0, 0.0))
    y0 = ts.bloch_sheet_tensor(cell, F, MESH, 1.05 * K0, 0.0)
    y90 = ts.bloch_sheet_tensor(cell, F, MESH, 0.0, 1.05 * K0)
    assert y0[0, 0] == pytest.approx(y90[1, 1], rel=1e-9)  # the longitudinal component, 90-degree mask symmetry
    assert y0[1, 1] == pytest.approx(y90[0, 0], rel=1e-9)  # the transverse component
    assert abs(y0[0, 1]) < 1e-9 * abs(y0[0, 0])
    # spatial dispersion: at kt > k0 even an isotropic patch has Y_longitudinal != Y_transverse
    assert abs(y0[0, 0].imag - y0[1, 1].imag) > 0.05 * abs(y0[0, 0].imag)


def test_slotted_circle_quarter_turn_swaps_axes_exactly():
    a = ts.bloch_sheet_tensor(_cell(SlottedCircularPatch(5.4e-3, 4.0e-3, 1.2e-3, 0.0)), F, MESH, 1.05 * K0, 0.0)
    b = ts.bloch_sheet_tensor(
        _cell(SlottedCircularPatch(5.4e-3, 4.0e-3, 1.2e-3, np.pi / 2)), F, MESH, 0.0, 1.05 * K0
    )
    # slot turned by 90 deg AND wave turned by 90 deg: the same physical situation, transposed axes
    assert b[1, 1] == pytest.approx(a[0, 0], rel=1e-8)
    assert b[0, 0] == pytest.approx(a[1, 1], rel=1e-8)
    assert a[0, 0].imag != pytest.approx(a[1, 1].imag, rel=0.05)  # and the slot does make it anisotropic


def test_bloch_dispersion_of_a_small_circle_matches_the_scalar_sheet_model():
    cell = _cell(SlottedCircularPatch(3.0e-3, 0.0, 0.0))
    beta, _y = ts.bloch_beta(cell, F, MESH, 0.0, n_grid=5, refine=1)
    y_normal = ts.sheet_admittance_tensor_for_cell(cell, F, MESH)[0, 0]
    ref = find_beta_from_sheet_admittance(y_normal, SUB, F).beta_rad_per_m
    assert beta == pytest.approx(ref, rel=2e-3)  # weak spatial dispersion for a small element
    assert K0 < beta < np.sqrt(2.2) * K0


def test_slot_makes_beta_depend_on_direction_and_normal_incidence_is_a_good_model_along_the_slot():
    cell = _cell(SlottedCircularPatch(5.6e-3, 4.6e-3, 1.2e-3, 0.0))
    b0, _ = ts.bloch_beta(cell, F, MESH, 0.0, n_grid=5, refine=1)
    b90, _ = ts.bloch_beta(cell, F, MESH, np.pi / 2, n_grid=5, refine=1)
    assert b90 > b0 * 1.003


def test_orders_assembly_equals_the_standard_assembly_for_equal_orders():
    cell = _cell(SlottedCircularPatch(4.0e-3, 0.0, 0.0))
    basis = RooftopBasis(cell, 16, 16)
    green = LayeredGreenFunction(SUB, F, backing="pec")
    z1 = assemble_impedance_matrix(basis, green, 0.3 * K0, 0.1 * K0, 12)
    z2 = ts.assemble_impedance_matrix_orders(basis, green, 0.3 * K0, 0.1 * K0, 12, 12)
    assert np.allclose(z1, z2, rtol=1e-12, atol=1e-18)


def test_cell_row_of_identical_cells_reproduces_the_single_cell_mean_sheet():
    one = SlottedCircularPatch(4.6e-3, 0.0, 0.0)
    row = UnitCell(RectangularLattice(2 * P, P), CellRow((one, one), P), SUB, Metal(), backing="pec")
    mesh_row = MeshSettings(32, 16, 32)
    y_row = ts.bloch_sheet_tensor(row, F, mesh_row, 1.05 * K0, 0.0, orders=(32, 16))
    y_one = ts.bloch_sheet_tensor(_cell(one), F, MeshSettings(16, 16, 16), 1.05 * K0, 0.0)
    assert y_row[0, 0].imag == pytest.approx(y_one[0, 0].imag, rel=0.03)  # same cells, doubled period


def test_slotted_square_keeps_square_outline_when_slot_turns():
    a = SlottedSquarePatch(5.0e-3, 3.0e-3, 0.8e-3, 0.0)
    b = SlottedSquarePatch(5.0e-3, 3.0e-3, 0.8e-3, 0.9)
    assert np.allclose(a.corners(), b.corners())
    assert a.feature_area == pytest.approx(25e-6 - 3.0e-3 * 0.8e-3)


def test_cross_slotted_circle_is_isotropic_and_quarter_turn_invariant():
    from src.dispersion.slotted_cells import CrossSlottedCircularPatch

    c = CrossSlottedCircularPatch(6.0e-3, 4.4e-3, 1.6e-3, 0.0)
    assert c.feature_area == pytest.approx(np.pi * 3.0e-3**2 - (2 * 4.4e-3 * 1.6e-3 - 1.6e-3**2), rel=1e-12)
    x = np.linspace(-3.2e-3, 3.2e-3, 129)
    xx, yy = np.meshgrid(x, x)
    m = c.in_feature(xx, yy)
    assert np.mean(m != m.T) < 0.005  # symmetric under x <-> y (grid-edge rounding only)
    assert np.mean(m != m[::-1, ::-1]) < 0.005
    cell = _cell(c)
    y = ts.bloch_sheet_tensor(cell, F, MESH, 0.0, 0.0)
    assert y[0, 0] == pytest.approx(y[1, 1], rel=1e-9)
    assert abs(y[0, 1]) < 1e-9 * abs(y[0, 0])
