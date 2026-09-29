"""Automatic mesh / Floquet validation: a priori checks and the refinement study."""

import warnings

import numpy as np
import pytest
from sweep_helpers import config_dict

from src.basis_functions import MeshSettings
from src.config import config_from_dict
from src.constants import C0
from src.discretization import (
    cells_per_wavelength,
    check_cell,
    check_frequency,
    grating_onset_frequency,
    is_specular_only,
    max_abs_difference,
    refinement_study,
)
from src.sparameters import compute_s_parameters


def cell_and_planes(**overrides):
    cfg = config_from_dict(config_dict(**overrides))
    return cfg.unit_cell, cfg.reference_planes, cfg.mesh


# --------------------------------------------------------------------------- cell checks


def test_aligned_rectangle_has_zero_area_error_and_the_counted_number_of_unknowns():
    cell, _, mesh = cell_and_planes(**{"geometry.length": 3e-3, "geometry.width": 2e-3})
    chk = check_cell(cell, mesh, 0.05)
    nx, ny = 6, 4  # 3 mm x 2 mm on a 0.5 mm mesh
    assert chk["metal_area_error"] < 1e-14 and chk["area_ok"] is True
    assert chk["n_unknowns"] == (nx - 1) * ny + nx * (ny - 1)  # interior x-edges + interior y-edges


def test_misaligned_rectangle_area_error_equals_the_independent_staircase_count():
    cell, _, mesh = cell_and_planes(**{"geometry.length": 5.3e-3})
    chk = check_cell(cell, mesh, 0.05)
    # centroids at odd multiples of 0.25 mm: |x| <= 2.65 keeps 10 columns (5.0 mm), width 2 mm exact
    assert chk["metal_area_error"] == pytest.approx((5.3e-3 - 5.0e-3) / 5.3e-3, rel=1e-12)
    assert chk["area_ok"] is False
    assert check_cell(cell, mesh, 0.06)["area_ok"] is True  # the threshold is the only difference


def test_slot_area_error_is_measured_on_the_metal_complement():
    cell, _, mesh = cell_and_planes(**{"geometry.type": "slot", "geometry.length": 5e-3, "geometry.width": 2e-3})
    chk = check_cell(cell, mesh, 0.05)
    assert chk["metal_area_error"] < 1e-14
    assert chk["n_unknowns"] > 0
    mis, _, mesh = cell_and_planes(**{"geometry.type": "slot", "geometry.length": 5.3e-3, "geometry.width": 2e-3})
    metal_true = 64e-6 - 5.3e-3 * 2e-3
    metal_mesh = 64e-6 - 5.0e-3 * 2e-3
    assert check_cell(mis, mesh, 1.0)["metal_area_error"] == pytest.approx(
        abs(metal_mesh - metal_true) / metal_true, rel=1e-9
    )


def test_mesh_too_coarse_for_the_pattern_is_an_error_not_a_bare_slab():
    cell, _, mesh = cell_and_planes(**{"geometry.length": 0.3e-3, "geometry.width": 0.3e-3})
    with pytest.raises(ValueError, match="no interior metal edges"):
        check_cell(cell, mesh, 0.05)


# ------------------------------------------------------------------- frequency checks


def test_cells_per_wavelength_matches_the_definition_and_scales_correctly():
    cell, _, mesh = cell_and_planes()
    f = 10e9
    n_re = np.sqrt(2.2 * (1 - 1j * 0.0009)).real
    expected = (C0 / f / n_re) / (8e-3 / 16)
    assert cells_per_wavelength(cell, mesh, f) == pytest.approx(expected, rel=1e-12)
    assert cells_per_wavelength(cell, mesh, 2 * f) == pytest.approx(expected / 2, rel=1e-12)
    finer = MeshSettings(32, 32, 32)
    assert cells_per_wavelength(cell, finer, f) == pytest.approx(2 * expected, rel=1e-12)


def test_cells_per_wavelength_uses_the_coarser_mesh_direction():
    cfg = config_from_dict(config_dict(**{"mesh.cells_x": 8, "mesh.cells_y": 16, "mesh.floquet_order": 16}))
    n_re = np.sqrt(2.2 * (1 - 1j * 0.0009)).real
    assert cells_per_wavelength(cfg.unit_cell, cfg.mesh, 10e9) == pytest.approx((C0 / 10e9 / n_re) / 1e-3, rel=1e-12)


def test_check_frequency_flags_wavelength_and_specular_condition():
    cell, _, mesh = cell_and_planes()
    ok = check_frequency(cell, mesh, 10e9, 10.0)
    assert ok["wavelength_ok"] and ok["specular_only"]
    assert not check_frequency(cell, mesh, 10e9, 100.0)["wavelength_ok"]
    assert not check_frequency(cell, mesh, 40e9, 10.0)["specular_only"]


@pytest.mark.parametrize("px, py", [(8e-3, 8e-3), (8e-3, 5e-3), (5e-3, 8e-3), (12e-3, 6e-3)])
def test_grating_onset_is_c_over_the_larger_period(px, py):
    cell, _, _ = cell_and_planes(period_x=px, period_y=py, **{"geometry.length": 3e-3, "geometry.width": 1e-3})
    f_c = grating_onset_frequency(cell)
    assert f_c == pytest.approx(C0 / max(px, py), rel=1e-14)
    assert is_specular_only(cell, 0.99 * f_c) is True
    assert is_specular_only(cell, 1.01 * f_c) is False


def test_specular_flag_agrees_with_the_solver_warning():
    """The flag must be true exactly when compute_s_parameters does not warn about grating orders."""
    cell, planes, _ = cell_and_planes()
    mesh = MeshSettings(8, 8, 8)
    f_c = grating_onset_frequency(cell)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compute_s_parameters(cell, 0.9 * f_c, mesh, planes)
    assert is_specular_only(cell, 0.9 * f_c)
    with pytest.warns(UserWarning, match="higher Floquet order"):
        compute_s_parameters(cell, 1.1 * f_c, mesh, planes)
    assert not is_specular_only(cell, 1.1 * f_c)


def test_warn_grating_false_silences_the_warning_without_changing_the_result():
    cell, planes, _ = cell_and_planes()
    mesh = MeshSettings(8, 8, 8)
    f = 1.1 * grating_onset_frequency(cell)
    with pytest.warns(UserWarning):
        loud = compute_s_parameters(cell, f, mesh, planes)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        quiet = compute_s_parameters(cell, f, mesh, planes, warn_grating=False)
    assert np.array_equal(loud.matrix(), quiet.matrix())


def test_above_the_grating_onset_the_specular_matrix_is_passive_but_not_unitary():
    """Power leaves in higher orders: the (0,0) S-matrix must lose norm even for a lossless cell."""
    cell, planes, _ = cell_and_planes(**{"substrate.tan_delta": 0.0})
    mesh = MeshSettings(16, 16, 32)
    f = 1.15 * grating_onset_frequency(cell)
    m = compute_s_parameters(cell, f, mesh, planes, warn_grating=False).matrix()
    eig = np.linalg.eigvalsh(np.eye(4) - m.conj().T @ m)
    assert eig.min() > -1e-9  # never creates power
    assert eig.max() > 1e-3  # but does lose some


# ---------------------------------------------------------------------- refinement study


def test_max_abs_difference_definition():
    a = np.array([[1.0, 2.0j], [0.0, 1.0]])
    b = a + np.array([[0.0, 0.0], [0.3 - 0.4j, 0.0]])
    assert max_abs_difference(a, b) == pytest.approx(0.5)
    assert max_abs_difference(a, a) == 0.0
    with pytest.raises(ValueError):
        max_abs_difference(a, np.zeros((3, 3)))


def test_refinement_study_is_the_difference_between_the_stated_discretizations():
    cell, planes, _ = cell_and_planes()
    mesh = MeshSettings(8, 8, 8)
    f = 10e9
    base = compute_s_parameters(cell, f, mesh, planes).matrix()
    a = compute_s_parameters(cell, f, MeshSettings(8, 8, 16), planes).matrix()
    b = compute_s_parameters(cell, f, MeshSettings(16, 16, 16), planes).matrix()
    floquet, mesh_err = refinement_study(cell, planes, mesh, f, base, 2)
    assert floquet == pytest.approx(np.max(np.abs(a - base)), rel=1e-12) and floquet > 1e-6
    assert mesh_err == pytest.approx(np.max(np.abs(b - a)), rel=1e-12) and mesh_err > 1e-6
    f3, m3 = refinement_study(cell, planes, mesh, f, base, 3)
    c = compute_s_parameters(cell, f, MeshSettings(8, 8, 24), planes).matrix()
    assert f3 == pytest.approx(np.max(np.abs(c - base)), rel=1e-12)


def test_full_pec_sheet_is_converged_exactly_at_any_discretization():
    """Analytic case: a cell that is all metal has S11 = -I and S21 = 0 for every mesh and order,
    so both refinement errors must vanish to rounding."""
    cell, planes, _ = cell_and_planes(**{"geometry.length": 8e-3, "geometry.width": 8e-3})
    for mesh in (MeshSettings(4, 4, 4), MeshSettings(8, 8, 8), MeshSettings(6, 10, 8)):
        base = compute_s_parameters(cell, 10e9, mesh, planes).matrix()
        floquet, mesh_err = refinement_study(cell, planes, mesh, 10e9, base, 2)
        assert floquet < 1e-10 and mesh_err < 1e-10
        assert np.allclose(base[:2, :2], -np.eye(2), atol=1e-9) and np.allclose(base[2:, :2], 0.0, atol=1e-9)


def test_floquet_error_decreases_with_the_order_at_a_fixed_mesh():
    cell, planes, _ = cell_and_planes()
    errors = []
    for order in (8, 16, 32):
        mesh = MeshSettings(16, 16, order)
        base = compute_s_parameters(cell, 10e9, mesh, planes).matrix()
        errors.append(refinement_study(cell, planes, mesh, 10e9, base, 2)[0])
    assert errors[0] > errors[1] > errors[2] > 0.0
    assert errors[0] / errors[1] > 2.0 and errors[1] / errors[2] > 2.0  # measured ~ order^-2


def test_mesh_error_decreases_with_the_mesh_at_a_fixed_floquet_to_mesh_ratio():
    cell, planes, _ = cell_and_planes()
    errors = []
    for n in (8, 16):
        mesh = MeshSettings(n, n, 2 * n)
        base = compute_s_parameters(cell, 10e9, mesh, planes).matrix()
        errors.append(refinement_study(cell, planes, mesh, 10e9, base, 2)[1])
    assert errors[0] > errors[1] > 0.0
    assert errors[0] / errors[1] > 1.3  # first order in the mesh size


@pytest.mark.parametrize("factor", [1, 0, -2, 2.0, 1.5, True, "2", None])
def test_refinement_factor_must_be_an_integer_of_at_least_two(factor):
    cell, planes, mesh = cell_and_planes()
    with pytest.raises(ValueError, match="integer >= 2"):
        refinement_study(cell, planes, mesh, 10e9, np.eye(4), factor)
