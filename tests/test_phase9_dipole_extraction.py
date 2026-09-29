"""Sec. 20: dipole-moment extraction reuses `basis_functions.RooftopBasis.evaluate` (unmodified);
cross-checked here against the basis's own closed-form Fourier-transform integral."""

import numpy as np
import pytest

from src.basis_functions import MeshSettings, RooftopBasis
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.phase9.dipole_extraction import local_dipole_moment

P = 4.8e-3


def _basis_with_currents(seed=0):
    lattice = RectangularLattice(P, P)
    geometry = Geometry(kind="rectangle", length=2.8e-3, width=2.8e-3)
    substrate = Substrate("test", 2.2, 0.0009, 1.575e-3)
    unit_cell = UnitCell(lattice, geometry, substrate, Metal())
    basis = RooftopBasis(unit_cell, 10, 10)
    rng = np.random.default_rng(seed)
    coeffs = rng.normal(size=basis.n_unknowns) + 1j * rng.normal(size=basis.n_unknowns)
    return basis, coeffs


def test_local_dipole_moment_over_the_whole_cell_matches_current_coefficients_at_dc():
    """Sec. 20: the whole-cell dipole moment must equal `current_coefficients(coeffs, 0, 0) *
    cell_area` -- the same DC Fourier-transform integral `basis_functions.py` already implements,
    computed here by a completely independent (real-space Riemann-sum) method."""
    basis, coeffs = _basis_with_currents()
    dipole = local_dipole_moment(basis, coeffs, 0.0, 0.0, P / 2, P / 2, n_samples=200)
    jx_dc, jy_dc = basis.current_coefficients(coeffs, 0.0, 0.0)
    expected_px = jx_dc * basis.lattice.cell_area
    expected_py = jy_dc * basis.lattice.cell_area
    assert dipole.px == pytest.approx(expected_px, rel=1e-2, abs=1e-12)
    assert dipole.py == pytest.approx(expected_py, rel=1e-2, abs=1e-12)


def test_local_dipole_moment_is_zero_for_zero_currents():
    basis, coeffs = _basis_with_currents()
    dipole = local_dipole_moment(basis, np.zeros_like(coeffs), 0.0, 0.0, P / 2, P / 2)
    assert dipole.px == 0.0
    assert dipole.py == 0.0


def test_local_dipole_moment_over_a_region_with_no_metal_is_zero():
    basis, coeffs = _basis_with_currents()
    # a small region right at the very corner of the cell, well outside the 2.8mm patch.
    dipole = local_dipole_moment(basis, coeffs, P / 2 - 0.1e-3, P / 2 - 0.1e-3, 0.05e-3, 0.05e-3, n_samples=8)
    assert abs(dipole.px) < 1e-15
    assert abs(dipole.py) < 1e-15
