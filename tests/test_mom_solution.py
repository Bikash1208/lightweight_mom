"""Solved MoM systems: analytic cases, PEC boundary condition, mesh convergence, symmetries."""

import math

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.config import config_from_dict
from src.constants import ETA0
from src.excitation import PlaneWaveExcitation
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.solver import assemble_system, build_system, solve
from test_mom_excitation import PolarizedStub

F = 10e9
PX = PY = 8e-3
LAT = RectangularLattice(PX, PY)


def cell(kind, length, width, eps_r=2.2, tand=0.0009, h=1.575e-3, **kw):
    return UnitCell(LAT, Geometry(kind, length, width, **kw), Substrate("t", eps_r, tand, h), Metal())


def system(c, nx, order=None, exc=None):
    exc = exc or PlaneWaveExcitation(F)
    return assemble_system(c, exc, MeshSettings(nx, nx, 4 * nx if order is None else order))


def mean_current(sys_, exc=None):
    i = solve(sys_.z, sys_.v)
    return sys_.basis.current_coefficients(i, sys_.kx0, sys_.ky0), i


# ==================================================================== analytic cases

@pytest.mark.parametrize("eps_r, tand, h", [(1.0, 0.0, 1.575e-3), (2.2, 0.0009, 1.575e-3),
                                             (10.2, 0.05, 3e-3), (2.2, 0.0, 5.05e-3)])
@pytest.mark.parametrize("pol", ["x", "y"])
@pytest.mark.parametrize("nx, order", [(4, 2), (8, 4), (12, 17)])   # order >= nx/2: Z full rank
def test_full_pec_sheet_carries_2_over_eta0_for_any_substrate_mesh_and_order(eps_r, tand, h, pol, nx, order):
    """Full-cell PEC on the interface: total E vanishes on z = 0, so the air side sees a PEC mirror
    and J = 2 H_inc = 2 E0 / eta0, independent of the substrate behind the sheet. Uniform current is
    exactly representable, so this holds to rounding for every mesh and Floquet order."""
    c = cell("patch", 8e-3, 8e-3, eps_r, tand, h)
    s = system(c, nx, order, PlaneWaveExcitation(F, polarization=pol))
    j, i = mean_current(s)
    main, other = (0, 1) if pol == "x" else (1, 0)
    assert j[main] == pytest.approx(2 / ETA0, rel=1e-11)
    assert abs(j[other]) < 1e-12 / ETA0
    on_main = s.basis.direction == main
    assert i[on_main] == pytest.approx(2 / ETA0, rel=1e-11)                # every rooftop, not just the mean
    assert np.all(np.abs(i[~on_main]) < 1e-11 / ETA0)


@pytest.mark.parametrize("theta_deg", [20.0, 50.0])
@pytest.mark.parametrize("tand", [0.0, 0.0009])
def test_full_pec_sheet_oblique_te_current(theta_deg, tand):
    """TE oblique: J = 2 H_tan = 2 E0 cos(theta) / eta0. (The rooftop interpolation of a phase-varying
    current is approximate, so refine the mesh and require convergence.)"""
    th = math.radians(theta_deg)
    exc = PlaneWaveExcitation(F, theta=th, phi=math.pi / 2, polarization="x")
    errs = []
    for nx in (8, 16):
        s = system(cell("patch", 8e-3, 8e-3, tand=tand), nx, 4 * nx, exc)
        j, _ = mean_current(s)
        errs.append(abs(j[0] - 2 * math.cos(th) / ETA0) * ETA0 / (2 * math.cos(th)))
    assert errs[1] < errs[0] < 1e-3
    assert errs[1] < 2e-4


@pytest.mark.parametrize("theta_deg", [20.0, 50.0])
@pytest.mark.parametrize("phi_deg", [0.0, 90.0])
def test_full_pec_sheet_oblique_tm_current_is_2_over_eta0_independent_of_angle(theta_deg, phi_deg):
    """TM oblique: H_tan = E0/eta0 whatever the angle (the cos(theta) of E_tan cancels Z_TM = eta0 cos(theta)),
    so J = 2/eta0 along the plane of incidence. Exercises the TM impedance and TM excitation together."""
    th, ph = math.radians(theta_deg), math.radians(phi_deg)
    exc = PolarizedStub(F, th, ph, tm=True)
    s = system(cell("patch", 8e-3, 8e-3), 16, 64, exc)
    j, _ = mean_current(s)
    along = j[0] * math.cos(ph) + j[1] * math.sin(ph)
    across = -j[0] * math.sin(ph) + j[1] * math.cos(ph)
    assert along == pytest.approx(2 / ETA0, rel=2e-4)
    assert abs(across) < 1e-4 / ETA0


def test_no_metal_no_unknowns():
    from src.basis_functions import RooftopBasis
    with pytest.raises(ValueError):
        RooftopBasis(cell("patch", 1e-3, 1e-3), 8, 8)


# ================================================================== PEC boundary condition

def total_tangential_field(s, i, points, order):
    """E_inc(stack) + E_scat at z = 0 at ``points`` (normal incidence), Floquet sum to ``order``."""
    kx, ky, g = s.green.floquet_dyads(s.basis.lattice, s.kx0, s.ky0, order)
    j = s.basis.current_coefficients(i, kx, ky)
    es = np.einsum("mnij,mnj->mni", g, j)
    x, y = points
    ex = np.exp(-1j * np.multiply.outer(x, kx[:, 0]))
    ey = np.exp(-1j * np.multiply.outer(y, ky[0, :]))
    scat = np.stack([np.einsum("abm,abn,mn->ab", ex, ey, es[..., c]) for c in (0, 1)], axis=-1)
    return scat + s.e_stack


def metal_cell_centres(b):
    cx = -PX / 2 + (np.arange(b.nx) + 0.5) * b.dx
    cy = -PY / 2 + (np.arange(b.ny) + 0.5) * b.dy
    return np.meshgrid(cx, cy, indexing="ij")


def test_full_cell_metal_has_zero_total_tangential_field_everywhere():
    c = cell("patch", 8e-3, 8e-3)
    s = system(c, 8, 16)
    _, i = mean_current(s)
    tot = total_tangential_field(s, i, metal_cell_centres(s.basis), 32)
    assert np.abs(tot).max() < 1e-12 * np.abs(s.e_stack).max()


def test_pec_condition_holds_in_the_weak_sense_and_converges_with_the_mesh():
    """Galerkin enforces E_tan = 0 weakly. For a 6 x 2 mm strip: (a) the mean |E_tot|/|E_inc| away from the
    strip ends and (b) a smooth-window (cosine) average of E_tot both fall roughly like 1/N when the
    mesh is refined. Pointwise values at the extreme end cells do NOT vanish (edge singularity of the
    true current) and are deliberately excluded from (a)."""
    L, W = 6e-3, 2e-3
    interior, weak = [], []
    for nx in (16, 24, 32):
        s = system(cell("rectangle", L, W), nx)
        _, i = mean_current(s)
        X, Y = metal_cell_centres(s.basis)
        tot = total_tangential_field(s, i, (X, Y), 8 * nx)[..., 0] / abs(s.e_stack[0])
        m = s.basis.metal_mask
        ii, jj = np.nonzero(m)
        k = int(round(1e-3 / s.basis.dx))
        sel = np.zeros_like(m)
        sel[ii.min() + k : ii.max() - k + 1, jj.min() : jj.max() + 1] = True
        interior.append(np.mean(np.abs(tot[sel])))
        w = np.cos(np.pi * X / L) * np.cos(np.pi * Y / W) * m
        weak.append(abs(np.sum(w * tot)) / np.sum(w))
    for series in (interior, weak):
        assert series[0] > series[1] > series[2]
        assert series[2] < 0.6 * series[0]
    assert weak[2] < 0.1 and interior[2] < 0.08
    # first-order-like: doubling the mesh roughly halves the residual (16 -> 32 here)
    assert 1.6 < interior[0] / interior[2] < 2.6 and 1.6 < weak[0] / weak[2] < 2.6


@pytest.mark.parametrize("kind, length, width", [("rectangle", 4e-3, 2e-3), ("slot", 6e-3, 2e-3)])
def test_solution_has_no_normal_current_on_metal_boundary_and_no_current_in_gaps(kind, length, width):
    c = cell(kind, length, width)
    s = system(c, 16)
    _, i = mean_current(s)
    rng = np.random.default_rng(3)
    # points in non-metal cells: zero current; at cell centres of the boundary edges: normal J = 0
    x = rng.uniform(-PX / 2, PX / 2, 4000)
    y = rng.uniform(-PY / 2, PY / 2, 4000)
    j = s.basis.evaluate(x, y, i)
    gap = ~c.is_metal(x, y)
    far = gap & ~c.is_metal(x + 0.4 * s.basis.dx, y) & ~c.is_metal(x - 0.4 * s.basis.dx, y) \
        & ~c.is_metal(x, y + 0.4 * s.basis.dy) & ~c.is_metal(x, y - 0.4 * s.basis.dy)
    assert far.sum() > 100
    assert np.abs(j[far]).max() == 0.0


# ============================================================= mesh / Floquet convergence

def test_mean_current_converges_with_mesh_refinement():
    """Zeroth Floquet coefficient of the solved current (what the far field / S-parameters are built
    from) for a 6 x 2 mm dipole-like patch, on aligned meshes 8/16/32 with Floquet order 4N.
    Rooftops at PEC edges converge ~first order: successive differences must shrink by ~2."""
    vals = []
    for nx in (8, 16, 32):
        s = system(cell("rectangle", 6e-3, 2e-3), nx)
        j, _ = mean_current(s)
        vals.append(j[0])
        assert abs(j[1]) < 1e-12 * abs(j[0])                    # y-symmetric cell: no cross-pol at k0
    d = [abs(vals[1] - vals[0]), abs(vals[2] - vals[1])]
    assert d[1] < d[0] and 1.5 < d[0] / d[1] < 3.0
    # extrapolated (first order) limit lies beyond the finest mesh and the finest is within ~10%
    limit = vals[2] + (vals[2] - vals[1])
    assert abs(vals[2] - limit) < 0.1 * abs(limit)


def test_mean_current_converges_with_mesh_refinement_for_slot_complement():
    """Same for a slot (metal-complement currents on a much larger domain)."""
    vals = []
    for nx in (8, 16, 24):
        s = system(cell("slot", 6e-3, 2e-3), nx, 3 * nx)
        j, _ = mean_current(s)
        vals.append(j[0])
    d = [abs(vals[1] - vals[0]), abs(vals[2] - vals[1])]
    assert d[1] < d[0]
    assert abs(vals[2] - vals[1]) < 0.1 * abs(vals[2])


def test_mean_current_converges_with_floquet_order():
    vals = []
    for order in (8, 16, 32, 64):
        s = system(cell("rectangle", 6e-3, 2e-3), 16, order)
        j, _ = mean_current(s)
        vals.append(j[0])
    d = np.abs(np.diff(vals))
    assert np.all(np.diff(d) < 0.0)
    assert d[-1] < 1e-2 * abs(vals[-1])


# ========================================================================= symmetries

@pytest.mark.parametrize("kind, length, width", [("rectangle", 6e-3, 2e-3), ("slot", 6e-3, 2e-3),
                                                  ("dipole", 8e-3, 1e-3)])
def test_x_polarized_solution_has_mirror_parity(kind, length, width):
    """Symmetric cell + uniform E_x: J_x even in x and y, J_y odd in x and in y."""
    c = cell(kind, length, width)
    s = system(c, 16, 32)
    _, i = mean_current(s)
    rng = np.random.default_rng(9)
    x = rng.uniform(-PX / 2 + 1e-6, PX / 2 - 1e-6, 500)
    y = rng.uniform(-PY / 2 + 1e-6, PY / 2 - 1e-6, 500)
    j = s.basis.evaluate(x, y, i)
    jxm = s.basis.evaluate(-x, y, i)
    jym = s.basis.evaluate(x, -y, i)
    scale = np.abs(j).max()
    assert np.abs(jxm[:, 0] - j[:, 0]).max() < 1e-9 * scale
    assert np.abs(jxm[:, 1] + j[:, 1]).max() < 1e-9 * scale
    assert np.abs(jym[:, 0] - j[:, 0]).max() < 1e-9 * scale
    assert np.abs(jym[:, 1] + j[:, 1]).max() < 1e-9 * scale


@pytest.mark.parametrize("kind, length, width", [("rectangle", 6e-3, 2e-3), ("slot", 6e-3, 2e-3)])
def test_rotating_cell_and_polarization_by_90_degrees_rotates_the_current(kind, length, width):
    """Square lattice, isotropic stack: cell rotated 90 deg + y polarization gives the 90-deg-rotated
    current:  J_B(R r) = R J_A(r),  R(x, y) = (-y, x)."""
    sa = system(cell(kind, length, width), 16, 32, PlaneWaveExcitation(F, polarization="x"))
    sb = system(cell(kind, length, width, rotation=math.pi / 2), 16, 32, PlaneWaveExcitation(F, polarization="y"))
    _, ia = mean_current(sa)
    _, ib = mean_current(sb)
    rng = np.random.default_rng(10)
    x = rng.uniform(-PX / 2 + 1e-6, PX / 2 - 1e-6, 400)
    y = rng.uniform(-PY / 2 + 1e-6, PY / 2 - 1e-6, 400)
    ja = sa.basis.evaluate(x, y, ia)
    jb = sb.basis.evaluate(-y, x, ib)                                    # at R r
    rotated = np.stack([-ja[:, 1], ja[:, 0]], axis=-1)                    # R J_A
    assert np.abs(jb - rotated).max() < 1e-8 * np.abs(ja).max()


def test_x_and_y_polarization_agree_for_a_square_symmetric_patch():
    """A 4 x 4 mm square patch is symmetric under x<->y: currents for x and y polarization map onto each other."""
    sx = system(cell("patch", 4e-3, 4e-3), 16, 32, PlaneWaveExcitation(F, polarization="x"))
    sy = system(cell("patch", 4e-3, 4e-3), 16, 32, PlaneWaveExcitation(F, polarization="y"))
    jx, _ = mean_current(sx)
    jy, _ = mean_current(sy)
    assert jx[0] == pytest.approx(jy[1], rel=1e-10)
    assert abs(jx[1]) < 1e-12 * abs(jx[0])


# =================================================== residual / passivity / config route

def test_solution_satisfies_the_linear_system_and_the_power_balance_identity():
    s = system(cell("rectangle", 6e-3, 2e-3), 16, 32)
    i = solve(s.z, s.v)
    assert np.linalg.norm(s.z @ i - s.v) < 1e-10 * np.linalg.norm(s.v)
    # power delivered by the currents to the scattered field is positive (lossy substrate + radiation)
    assert np.real(np.conj(i) @ s.z @ i) > 0.0
    assert np.real(np.conj(i) @ s.v) == pytest.approx(np.real(np.conj(i) @ s.z @ i), rel=1e-9)


def test_config_route_builds_a_consistent_solvable_system():
    from test_mom_excitation import config_dict
    cfg = config_from_dict(config_dict(10e-3, 10e-3))
    s = build_system(cfg)
    i = solve(s.z, s.v)
    assert s.z.shape == (s.basis.n_unknowns,) * 2 and s.v.shape == (s.basis.n_unknowns,)
    assert np.all(np.isfinite(i))
    direct = assemble_system(cfg.unit_cell, cfg.excitation, cfg.mesh)
    assert np.array_equal(direct.z, s.z) and np.array_equal(direct.v, s.v)


def test_solve_validation():
    with pytest.raises(ValueError):
        solve(np.eye(3), np.ones(2))
    with pytest.raises(ValueError):
        solve(np.ones((2, 3)), np.ones(2))
    with pytest.raises(np.linalg.LinAlgError):
        solve(np.zeros((2, 2)), np.ones(2))
    assert solve(2 * np.eye(2), np.array([2.0, 4.0j])) == pytest.approx([1.0, 2.0j])


def test_solve_handles_several_right_hand_sides_with_one_factorization():
    s = system(cell("rectangle", 6e-3, 2e-3), 8, 16)
    rng = np.random.default_rng(12)
    v = np.stack([s.v, rng.normal(size=s.v.size) + 1j * rng.normal(size=s.v.size)], axis=1)
    both = solve(s.z, v)
    assert both.shape == (s.v.size, 2)
    for k in range(2):
        assert both[:, k] == pytest.approx(solve(s.z, v[:, k]), rel=1e-12)
    with pytest.raises(ValueError):
        solve(s.z, np.ones((s.v.size + 1, 2)))
    with pytest.raises(ValueError):
        solve(s.z, np.ones((s.v.size, 2, 2)))
