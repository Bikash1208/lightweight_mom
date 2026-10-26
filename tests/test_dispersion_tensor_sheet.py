"""Phase 26 B: tensor sheet-admittance extraction and tensor transverse resonance."""

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.dispersion import tensor_sheet as ts
from src.dispersion.sheet_model import (
    find_beta_from_sheet_admittance,
    sheet_admittance_from_gamma,
    sheet_admittance_from_normal_incidence,
)
from src.errors import LwmomError
from src.frequency import free_space_wavenumber
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

F = 10.54e9
P = 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
MESH = MeshSettings(32, 32, 32)


def _cell(length, width, angle_deg=0.0):
    geo = Geometry("rectangle", length, width, rotation=np.radians(angle_deg))
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")


@pytest.fixture(scope="module")
def tensors():
    return {
        "rect": ts.sheet_admittance_tensor_for_cell(_cell(4.0e-3, 2.4e-3), F, MESH),
        "rect90": ts.sheet_admittance_tensor_for_cell(_cell(4.0e-3, 2.4e-3, 90.0), F, MESH),
        "rect_swapped": ts.sheet_admittance_tensor_for_cell(_cell(2.4e-3, 4.0e-3), F, MESH),
        "rect45": ts.sheet_admittance_tensor_for_cell(_cell(4.0e-3, 2.4e-3, 45.0), F, MESH),
        "square": ts.sheet_admittance_tensor_for_cell(_cell(3.0e-3, 3.0e-3), F, MESH),
        "square45": ts.sheet_admittance_tensor_for_cell(_cell(3.0e-3, 3.0e-3, 45.0), F, MESH),
    }


def test_square_patch_is_isotropic_and_matches_the_scalar_extraction(tensors):
    y = tensors["square"]
    assert abs(y[0, 0] - y[1, 1]) < 1e-9 * abs(y[0, 0])
    assert abs(y[0, 1]) < 1e-9 * abs(y[0, 0])
    scalar = sheet_admittance_from_normal_incidence(_cell(3.0e-3, 3.0e-3), F, MESH)
    assert y[0, 0] == pytest.approx(scalar, rel=1e-9)


def test_diagonal_square_is_still_isotropic(tensors):
    y = tensors["square45"]
    assert abs(y[0, 0] - y[1, 1]) < 1e-6 * abs(y[0, 0])
    assert abs(y[0, 1]) < 1e-6 * abs(y[0, 0])


def test_reciprocity_symmetric_tensor(tensors):
    for name in ("rect", "rect45", "rect90"):
        assert ts.reciprocity_error(tensors[name]) < 1e-9


def test_quarter_turn_equals_exchanging_length_and_width(tensors):
    # the 90-degree-rotated mask IS the swapped-axes mask on this mesh: bit-for-bit up to rounding
    assert np.allclose(tensors["rect90"], tensors["rect_swapped"], rtol=1e-9, atol=1e-15)
    assert ts.rotation_equivariance_error(tensors["rect90"], tensors["rect"], np.pi / 2) < 1e-9


def test_diagonal_cell_has_principal_axes_on_the_diagonal_with_the_rotated_anisotropy(tensors):
    yi, delta, axis, asym = ts.decompose(tensors["rect45"])
    assert axis == pytest.approx(np.pi / 4, abs=1e-6)
    _yi0, delta0, axis0, _ = ts.decompose(tensors["rect"])
    assert min(axis0, np.pi - axis0) == pytest.approx(0.0, abs=1e-9)  # axis is defined mod pi
    assert delta == pytest.approx(delta0, rel=0.15)  # lattice effects, quantified by the sweep
    assert asym < 1e-9 * abs(yi)


def test_scalar_limit_reproduces_sheet_model_exactly(tensors):
    y = tensors["square"]
    scalar = sheet_admittance_from_normal_incidence(_cell(3.0e-3, 3.0e-3), F, MESH)
    ref = find_beta_from_sheet_admittance(scalar, SUB, F).beta_rad_per_m
    beta, character = ts.find_beta_tensor(scalar * np.eye(2), 0.0, SUB, F)
    assert beta == ref  # bit-identical
    assert character == 0.0
    for phi in (0.3, 1.1):  # rotated frame: identical to rounding of cos^2 + sin^2
        assert ts.find_beta_tensor(scalar * np.eye(2), phi, SUB, F)[0] == pytest.approx(ref, rel=1e-12)
    assert ts.find_beta_tensor(y, 0.7, SUB, F)[0] == pytest.approx(find_beta_from_sheet_admittance(y[0, 0], SUB, F).beta_rad_per_m, rel=1e-9)


def test_anisotropic_sheet_has_direction_dependent_beta_between_the_principal_values(tensors):
    y = tensors["rect"]
    b0 = ts.find_beta_tensor(y, 0.0, SUB, F)[0]
    b90 = ts.find_beta_tensor(y, np.pi / 2, SUB, F)[0]
    b45 = ts.find_beta_tensor(y, np.pi / 4, SUB, F)[0]
    assert b0 != b90
    assert min(b0, b90) < b45 < max(b0, b90)
    # principal directions decouple TM and TE: no hybrid content; at 45 deg the mode is a hybrid
    assert ts.find_beta_tensor(y, 0.0, SUB, F)[1] < 1e-12
    assert ts.find_beta_tensor(y, np.pi / 4, SUB, F)[1] > 1e-3
    # principal-axis betas are the scalar sheet_model betas of Y_xx and Y_yy
    assert b0 == pytest.approx(find_beta_from_sheet_admittance(y[0, 0], SUB, F).beta_rad_per_m, rel=1e-9)
    assert b90 == pytest.approx(find_beta_from_sheet_admittance(y[1, 1], SUB, F).beta_rad_per_m, rel=1e-9)


def test_dispersion_is_symmetric_under_a_half_turn_and_reflection(tensors):
    y = tensors["rect"]
    b = [ts.find_beta_tensor(y, p, SUB, F)[0] for p in (0.4, 0.4 + np.pi, -0.4, np.pi - 0.4)]
    assert b[1] == pytest.approx(b[0], rel=1e-9) and b[2] == pytest.approx(b[0], rel=1e-9)
    assert b[3] == pytest.approx(b[0], rel=1e-9)


def test_no_metal_tensor_is_zero_and_gives_bare_tm0():
    gamma = np.exp(1j * np.radians(138.0)) * np.eye(2)
    y = ts.sheet_admittance_tensor_from_gamma(gamma, SUB, F)
    ref = sheet_admittance_from_gamma(np.exp(1j * np.radians(138.0)), SUB, F)
    assert y[0, 0] == pytest.approx(ref) and y[1, 1] == pytest.approx(ref) and y[0, 1] == 0
    k0 = free_space_wavenumber(F)
    beta = ts.find_beta_tensor(0.0 * np.eye(2), 0.0, SUB, F)[0]
    assert k0 < beta < np.sqrt(2.2) * k0


def test_no_root_is_an_error_not_a_number():
    with pytest.raises(LwmomError):
        ts.find_beta_tensor(-1.0j * np.eye(2), 0.0, SUB, F)  # a very strongly inductive sheet: no bound TM0 in the interval
