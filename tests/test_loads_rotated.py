"""Phase 20 R3 item 2: loads whose conduction axis is not x or y (rotated rectangle region, sheet-impedance tensor)."""

import math
import warnings
from dataclasses import replace

import numpy as np
import pytest

from src.basis_functions import MeshSettings, RooftopBasis
from src.composite_geometry import AxisRect, CompositeGeometry
from src.excitation import ReferencePlanes
from src.geometry import UnitCell
from src.loads import LumpedLoad, _add_region, _add_rotated, absorbed_power, check_load_geometry, load_matrix
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import compute_s_parameters_oblique

SUB = Substrate("rt5880", 2.2, 0.0, 1.575e-3)
PLANES = ReferencePlanes(0.0, None, SUB.thickness)
W = 2 * math.pi * 10e9


def _cell(loads=(), half=4e-3):
    geo = CompositeGeometry((AxisRect(0.0, 0.0, 2 * half, 2 * half),))
    return UnitCell(RectangularLattice(12e-3, 12e-3), geo, SUB, Metal(), backing="pec", loads=tuple(loads))


def test_numeric_direction_canonicalisation_is_the_axis_aligned_load():
    a, b, c = (LumpedLoad(0, 0, 1e-3, 1e-3, d, R=3.0) for d in (0, 90, 180))
    assert a.direction == "x" and b.direction == "y" and c.direction == "x" and not a.is_rotated
    assert a == LumpedLoad(0, 0, 1e-3, 1e-3, "x", R=3.0)
    r = LumpedLoad(0, 0, 1e-3, 1e-3, 45, R=3.0)
    assert r.is_rotated and r.direction == 45.0 and r.angle == pytest.approx(math.pi / 4)
    assert LumpedLoad(0, 0, 1e-3, 1e-3, -45, R=3.0).direction == 135.0
    with pytest.raises(ValueError):
        LumpedLoad(0, 0, 1e-3, 1e-3, "z")
    basis = RooftopBasis(_cell(), 24, 24)
    assert np.array_equal(load_matrix(basis, [a], W), load_matrix(basis, [LumpedLoad(0, 0, 1e-3, 1e-3, "x", R=3.0)], W))


@pytest.mark.parametrize("angle", [30.0, 45.0, 135.0])
@pytest.mark.parametrize("n", [24, 48])
def test_uniform_current_along_the_axis_dissipates_one_half_R_I_total_squared(angle, n):
    """J = u (unit current density along the axis) inside the region: x-rooftop coefficients cos(a), y-rooftop sin(a) partition unity."""
    ld = LumpedLoad(0.5e-3, -0.25e-3, 1.5e-3, 0.8e-3, angle, R=7.0, transverse_conducts=False)
    basis = RooftopBasis(_cell(), n, n)
    i = np.zeros(basis.n_unknowns, complex)
    i[: basis.n_x] = math.cos(ld.angle)
    i[basis.n_x :] = math.sin(ld.angle)
    p = absorbed_power(basis, [ld], W, i)
    # the pointwise rooftop sum is constant, so the Gauss rule is exact; the 1e12-ohm transverse penalty adds ~1e-6 of rounding noise
    assert p == pytest.approx(0.5 * 7.0 * ld.width**2, rel=1e-4)
    ld_c = replace(ld, transverse_conducts=True)
    assert absorbed_power(basis, [ld_c], W, i) == pytest.approx(0.5 * 7.0 * ld.width**2, rel=1e-12)


def test_rotated_path_at_zero_angle_agrees_with_the_exact_axis_aligned_integrals():
    basis = RooftopBasis(_cell(), 48, 48)
    ld = LumpedLoad(0.5e-3, 0.0, 1.5e-3, 0.75e-3, "x", R=5.0, L=1e-9)
    zs = ld.impedance(W) * ld.sheet_factor
    exact = np.zeros((basis.n_unknowns,) * 2, complex)
    _add_region(exact, basis, "x", zs, ld.region())
    rot = np.zeros_like(exact)
    _add_rotated(rot, basis, replace(ld, direction=1e-9), zs)  # a (numerically) zero-angle rotated load
    nz = np.abs(exact).max()
    assert (
        np.abs(rot - exact).max() < 2e-3 * nz
    )  # composite Gauss with sub-cells h/4: not exact because the edges cut cells


def test_mirror_symmetry_of_the_physical_response():
    """Mirroring x -> -x maps a load at angle a at (x0, y0) to angle 180 - a at (-x0, y0): the cell is mirror symmetric, so the
    co-polar reflection is unchanged and the cross-polar one only changes sign."""

    def run(x0, ang):
        ld = LumpedLoad(x0, 1e-3, 1.5e-3, 0.8e-3, ang, R=40.0, L=0.3e-9, transverse_conducts=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return compute_s_parameters_oblique(_cell([ld]), 10e9, MeshSettings(48, 48, 25), PLANES, 0.0, 0.0).s11

    a, b = run(1.5e-3, 45), run(-1.5e-3, 135)
    assert abs(a[0, 0] - b[0, 0]) < 1e-9 and abs(a[1, 1] - b[1, 1]) < 1e-9
    assert abs(a[0, 1] + b[0, 1]) < 1e-9 and abs(a[1, 0] + b[1, 0]) < 1e-9
    assert abs(a[0, 1]) > 1e-4  # a genuinely rotated load couples the polarisations


def test_region_helpers_and_geometry_checks():
    ld = LumpedLoad(1e-3, 2e-3, 2e-3, 1e-3, 45.0, R=1.0)
    c = np.array(ld.corners())
    assert c.shape == (4, 2) and np.allclose(c.mean(axis=0), [1e-3, 2e-3])
    assert ld.contains(1e-3, 2e-3) and not ld.contains(1e-3 + 1.2e-3, 2e-3 - 1.2e-3 + 2e-3)
    x0, x1, y0, y1 = ld.region()
    assert x0 == pytest.approx(c[:, 0].min()) and y1 == pytest.approx(c[:, 1].max())
    cell = _cell(half=1e-3)  # patch too small for the region: the region is not metal
    with pytest.warns(UserWarning, match="not metal"):
        check_load_geometry([LumpedLoad(2.5e-3, 0.0, 2e-3, 1e-3, 45.0, R=1.0)], cell, 24, 24)


def test_blocking_transverse_conduction_on_a_staircased_diagonal_strip_kills_the_line():
    """Documented numerical limit: on a 45-degree strip of a few cells the Cartesian rooftops need J_x = J_y exactly to satisfy the block,
    so `transverse_conducts=False` also kills the longitudinal current. The CP-2 cell therefore uses conducting diode regions."""
    from dataclasses import replace as rep

    from src.ris.cells import ON, Cp2Cell, ideal_kit

    def co_pol_spread(tc):
        cell = rep(Cp2Cell(), transverse_conducts=tc)
        uc = cell.unit_cell((OFF_ := "OFF", ON, OFF_, OFF_), ideal_kit())  # the 45-degree stub ON
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s = compute_s_parameters_oblique(uc, 10e9, MeshSettings(54, 54, 28), PLANES, 0.0, 0.0).s11
        return abs(s[0, 1])  # x/y coupling of a 45-degree anisotropy axis

    assert co_pol_spread(True) > 20 * co_pol_spread(False)


# ---------------------------------------------------------------------------------------------------------------------------------------------
# Phase 28 (additive): blocked-transverse rotated loads lock the rooftop basis (documented R3 limit, now also a warning); conducting rotated loads are rotation covariant
def _strap_geo():
    """Two patches joined only by a 1.2 mm x 1.2 mm strap: a load on the strap is the only conduction path between them."""
    return CompositeGeometry((AxisRect(-2.1e-3, 0.0, 3.0e-3, 7.2e-3), AxisRect(2.1e-3, 0.0, 3.0e-3, 7.2e-3), AxisRect(0.0, 0.0, 1.2e-3, 1.2e-3)))


def _sheet(ld_list, mesh=45, geo=None):
    from src.hologram.coded_cell_eval import sheet_tensor

    geo = geo or CompositeGeometry((AxisRect(0.0, 0.0, 7.2e-3, 7.2e-3),))
    uc = UnitCell(RectangularLattice(9e-3, 9e-3), geo, Substrate("rt5880", 2.2, 0.0, 3.175e-3), Metal(), backing="pec", loads=tuple(ld_list))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sheet_tensor(uc, 10e9, mesh)


def test_blocked_rotated_load_locks_the_line_and_warns():
    """Minimal reproduction of the Phase 28 finding. Two patches joined only by a strap that carries a 1 milliohm load (a closed switch). Axis-aligned or conducting loads leave the
    link closed. A ROTATED load with the 1e12 ohm transverse penalty locks the Cartesian rooftops (no x/y combination has a zero transverse component), cuts the link and the cell tensor
    jumps towards that of the open gap although the load is a 1 milliohm closed switch (when the load spans the whole strip, as in the Phase 28 loop cell, ON == OFF exactly)."""
    geo = _strap_geo()
    closed = _sheet([], geo=geo)
    mk = lambda d, **kw: LumpedLoad(0.0, 0.0, 0.6e-3, 0.6e-3, d, R=1e-3, **kw)  # noqa: E731
    rel = lambda t: np.linalg.norm(t - closed) / np.linalg.norm(closed)  # noqa: E731
    axis = rel(_sheet([mk("x", transverse_conducts=False)], geo=geo))
    free = rel(_sheet([mk(20.0, transverse_conducts=True)], geo=geo))
    locked = rel(_sheet([mk(20.0, transverse_conducts=False)], geo=geo))
    assert axis < 1e-3 and free < 1e-3 and locked > 0.03 and locked > 100 * axis
    uc = UnitCell(RectangularLattice(9e-3, 9e-3), geo, SUB, Metal(), backing="pec", loads=(mk(20.0, transverse_conducts=False),))
    with pytest.warns(UserWarning, match="locks the rooftop basis"):
        check_load_geometry(uc.loads, uc, 45, 45)
    uc2 = replace(uc, loads=(mk(20.0, transverse_conducts=True),))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        check_load_geometry(uc2.loads, uc2, 45, 45)


@pytest.mark.parametrize("theta", [30.0, 45.0])
def test_a_conducting_rotated_load_changes_the_cell_tensor_like_the_rotated_unrotated_load(theta):
    """Physics: a load turned by theta in the middle of a solid patch changes the sheet tensor by R(theta) dY(0) R(theta)^T, to the rooftop staircase of the rotated region
    (2 % at 45 deg, 7 % at 30 deg at this mesh)."""
    base = _sheet([])
    d0 = _sheet([LumpedLoad(0.0, 0.0, 3.6e-3, 1.2e-3, "x", R=300.0)]) - base
    dt = _sheet([LumpedLoad(0.0, 0.0, 3.6e-3, 1.2e-3, theta, R=300.0, transverse_conducts=True)]) - base
    a = math.radians(theta)
    r = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    pred = r @ d0 @ r.T
    assert np.linalg.norm(d0) > 1e-4  # the load is visible
    assert np.linalg.norm(dt - pred) / np.linalg.norm(pred) < 0.12
