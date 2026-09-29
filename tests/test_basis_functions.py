"""Rooftop basis: mesh, support, continuity, periodic wrap, transform."""

import math

import numpy as np
import pytest

from src.basis_functions import MeshSettings, RooftopBasis
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

PX = PY = 8e-3
SUB = Substrate("t", 2.2, 0.0009, 1.575e-3)
LAT = RectangularLattice(PX, PY)


def cell(kind, length, width, **kw):
    return UnitCell(LAT, Geometry(kind, length, width, **kw), SUB, Metal())


def unit_vector(basis, index):
    e = np.zeros(basis.n_unknowns, complex)
    e[index] = 1.0
    return e


# ------------------------------------------------------------------ mesh / counts

def test_mesh_settings_validation():
    MeshSettings(2, 2, 1)
    MeshSettings(8, 8, 4)                       # 2*4 + 1 = 9 >= 8: minimum full-rank truncation
    for kw in (dict(cells_x=1), dict(cells_y=0), dict(floquet_order=-1),
               dict(cells_x=4.0), dict(cells_y=True), dict(floquet_order="3")):
        args = dict(cells_x=4, cells_y=4, floquet_order=2)
        args.update(kw)
        with pytest.raises(ValueError):
            MeshSettings(**args)


def test_floquet_order_must_resolve_the_mesh():
    """2*order + 1 >= max(cells): otherwise harmonics alias at the rooftop centres and Z is singular."""
    for cells, order in [(8, 3), (16, 7), (9, 3), (4, 0)]:
        with pytest.raises(ValueError, match="too small"):
            MeshSettings(cells, cells, order)
    MeshSettings(9, 4, 4)                       # anisotropic: uses the larger of the two
    with pytest.raises(ValueError, match="too small"):
        MeshSettings(4, 20, 8)


def test_unknown_counts_for_aligned_rectangle():
    """4 x 2 mm patch on a 1 mm mesh = 4 x 2 metal cells.
    x-edges between horizontally adjacent metal cells: 3 per row x 2 rows = 6;
    y-edges between vertically adjacent metal cells: 4 per column x 1 = 4."""
    b = RooftopBasis(cell("rectangle", 4e-3, 2e-3), 8, 8)
    assert b.metal_mask.sum() == 8
    assert (b.n_x, b.n_y, b.n_unknowns) == (6, 4, 10)
    assert b.dx == pytest.approx(1e-3, rel=1e-15) and b.dy == pytest.approx(1e-3, rel=1e-15)


def test_ordering_x_directed_first_and_centers_on_mesh_lines():
    b = RooftopBasis(cell("rectangle", 4e-3, 2e-3), 8, 8)
    assert np.all(b.direction[: b.n_x] == 0) and np.all(b.direction[b.n_x :] == 1)
    xm = (b.center_x[: b.n_x] + PX / 2) / b.dx           # x-directed: x on a mesh line (integer)
    ym = (b.center_y[: b.n_x] + PY / 2) / b.dy - 0.5     # ... y at a cell centre (half-integer)
    assert np.allclose(xm, np.round(xm)) and np.allclose(ym, np.round(ym))
    xm = (b.center_x[b.n_x :] + PX / 2) / b.dx - 0.5     # y-directed: swapped
    ym = (b.center_y[b.n_x :] + PY / 2) / b.dy
    assert np.allclose(xm, np.round(xm)) and np.allclose(ym, np.round(ym))


def test_every_unknown_sits_between_two_metal_cells_and_none_elsewhere():
    """PEC condition on J: no rooftop straddles a metal/non-metal boundary."""
    for kind, L, W in [("rectangle", 4e-3, 2e-3), ("slot", 6e-3, 2e-3), ("dipole", 8e-3, 1e-3)]:
        c = cell(kind, L, W)
        b = RooftopBasis(c, 8, 8)
        step = np.array([b.dx, b.dy])
        for k in range(b.n_unknowns):
            centre = np.array([b.center_x[k], b.center_y[k]])
            normal = np.eye(2)[b.direction[k]]
            for sign in (-1.0, 1.0):
                p = centre + sign * normal * step / 2
                assert c.is_metal(*p)                       # both neighbouring cell centroids
        # and conversely: count all metal-metal edges by brute force (with periodic wrap)
        m = b.metal_mask
        expect_x = sum(m[i, j] and m[(i - 1) % 8, j] for i in range(8) for j in range(8))
        expect_y = sum(m[i, j] and m[i, (j - 1) % 8] for i in range(8) for j in range(8))
        assert (b.n_x, b.n_y) == (expect_x, expect_y)


def test_slot_is_represented_by_metal_complement():
    r = RooftopBasis(cell("rectangle", 6e-3, 2e-3), 8, 8)
    s = RooftopBasis(cell("slot", 6e-3, 2e-3), 8, 8)
    assert np.array_equal(s.metal_mask, ~r.metal_mask)
    assert s.metal_mask.sum() == 64 - 12
    # no rooftop inside or on the slot
    slot = Geometry("slot", 6e-3, 2e-3)
    for k in range(s.n_unknowns):
        assert not slot.in_feature(s.center_x[k], s.center_y[k]) or (
            abs(s.center_x[k]) >= 3e-3 - 1e-12 or abs(s.center_y[k]) >= 1e-3 - 1e-12
        )
    # wrap-around: the complement is connected across the cell boundary, so it has edges
    # on the periodic seam (x = -Px/2 and y = -Py/2)
    assert np.any(np.isclose(s.center_x[: s.n_x], -PX / 2))
    assert np.any(np.isclose(s.center_y[s.n_x :], -PY / 2))


def test_isolated_patch_has_no_seam_edges_but_full_length_dipole_has():
    p = RooftopBasis(cell("patch", 4e-3, 4e-3), 8, 8)
    assert not np.any(np.isclose(p.center_x[: p.n_x], -PX / 2))
    d = RooftopBasis(cell("dipole", 8e-3, 1e-3), 8, 8)     # spans the full period along x
    assert np.any(np.isclose(d.center_x[: d.n_x], -PX / 2))
    assert d.n_x == 8 and d.n_y == 0                        # one full row of 8 wrapped edges


def test_metal_too_small_for_mesh_is_an_error():
    with pytest.raises(ValueError, match="no interior metal edges"):
        RooftopBasis(cell("patch", 1e-3, 1e-3), 8, 8)       # a single metal cell has no interior edge


def test_wrong_length_coefficients_rejected():
    b = RooftopBasis(cell("rectangle", 4e-3, 2e-3), 8, 8)
    with pytest.raises(ValueError):
        b.evaluate(0.0, 0.0, np.ones(b.n_unknowns + 1))
    with pytest.raises(ValueError):
        b.current_coefficients(np.ones(3), 0.0, 0.0)


# ------------------------------------------------------------------------- support

@pytest.mark.parametrize("which", ["x", "y"])
def test_single_rooftop_support_peak_and_linearity(which):
    b = RooftopBasis(cell("rectangle", 4e-3, 4e-3), 8, 8)
    k = 0 if which == "x" else b.n_x
    xc, yc = b.center_x[k], b.center_y[k]
    nrm = 0 if which == "x" else 1          # axis along which the rooftop varies linearly
    tgt = 1 - nrm                           # axis along which it is a pulse
    d_norm = b.dx if which == "x" else b.dy
    d_tan = b.dy if which == "x" else b.dx
    e = unit_vector(b, k)

    def J(dn, dt):
        pt = np.array([xc, yc])
        pt[nrm] += dn
        pt[tgt] += dt
        return b.evaluate(pt[0], pt[1], e)

    comp = 0 if which == "x" else 1
    other = 1 - comp
    assert J(0.0, 0.0)[comp] == pytest.approx(1.0, abs=1e-14)            # peak = 1 at the edge
    assert J(0.0, 0.0)[other] == 0.0
    for frac in (0.25, 0.5, 0.75):                                         # linear ramps both sides
        assert J(frac * d_norm, 0.0)[comp] == pytest.approx(1 - frac, abs=1e-14)
        assert J(-frac * d_norm, 0.0)[comp] == pytest.approx(1 - frac, abs=1e-14)
    # support: |normal| < d_norm and |tangential| < d_tan / 2, zero outside
    assert abs(J(1.001 * d_norm, 0.0)[comp]) == 0.0
    assert abs(J(0.0, 0.51 * d_tan)[comp]) == 0.0
    assert abs(J(0.0, -0.51 * d_tan)[comp]) == 0.0
    assert abs(J(0.0, 0.49 * d_tan)[comp]) > 0.0
    # the far edges of the support carry zero current (continuity to zero)
    assert abs(J(d_norm * (1 - 1e-9), 0.0)[comp]) < 1e-8


def test_support_is_exactly_two_mesh_cells():
    b = RooftopBasis(cell("rectangle", 4e-3, 4e-3), 16, 16)
    k = b.n_x // 2
    e = unit_vector(b, k)
    x = np.linspace(-PX / 2, PX / 2, 321)[:-1] + 1e-9
    X, Y = np.meshgrid(x, x, indexing="ij")
    J = b.evaluate(X, Y, e)
    nz = np.abs(J[..., 0]) > 0
    xs, ys = X[nz], Y[nz]
    assert xs.min() >= b.center_x[k] - b.dx and xs.max() <= b.center_x[k] + b.dx
    assert ys.min() >= b.center_y[k] - b.dy / 2 - 1e-9 and ys.max() <= b.center_y[k] + b.dy / 2
    assert np.all(J[..., 1] == 0.0)
    # integral of the rooftop = dx * dy (triangle area dx times pulse height 1 * dy)
    area = np.sum(J[..., 0].real) * (x[1] - x[0]) ** 2
    assert area == pytest.approx(b.dx * b.dy, rel=2e-2)


# ------------------------------------------------------------------ continuity

def test_normal_component_continuous_across_every_mesh_line():
    """Divergence-conforming: J_x is continuous across x = const lines (J_y across y = const),
    for arbitrary coefficients; the tangential component may jump."""
    b = RooftopBasis(cell("rectangle", 6e-3, 4e-3), 16, 16)
    rng = np.random.default_rng(2)
    coeffs = rng.normal(size=b.n_unknowns) + 1j * rng.normal(size=b.n_unknowns)
    eps = 1e-9
    lipschitz_tol = 10 * eps / b.dx * np.abs(coeffs).max()     # |dJ| <= |c| * 2*eps/dx across the gap
    y = rng.uniform(-PY / 2 + 0.1e-3, PY / 2 - 0.1e-3, 40)
    for line in range(1, 16):
        xl = -PX / 2 + line * b.dx
        left = b.evaluate(np.full(y.size, xl - eps), y, coeffs)[:, 0]
        right = b.evaluate(np.full(y.size, xl + eps), y, coeffs)[:, 0]
        assert np.allclose(left, right, rtol=0, atol=lipschitz_tol)
    x = rng.uniform(-PX / 2 + 0.1e-3, PX / 2 - 0.1e-3, 40)
    for line in range(1, 16):
        yl = -PY / 2 + line * b.dy
        lo = b.evaluate(x, np.full(x.size, yl - eps), coeffs)[:, 1]
        hi = b.evaluate(x, np.full(x.size, yl + eps), coeffs)[:, 1]
        assert np.allclose(lo, hi, rtol=0, atol=lipschitz_tol)


def test_tangential_component_jumps_between_cell_rows():
    """Rooftop tangential component is a pulse (discontinuous): documents the (standard) behaviour."""
    b = RooftopBasis(cell("rectangle", 4e-3, 4e-3), 8, 8)
    k = 0
    e = unit_vector(b, k)
    y_edge = b.center_y[k] + b.dy / 2
    inside = b.evaluate(b.center_x[k], y_edge - 1e-9, e)[0].real
    outside = b.evaluate(b.center_x[k], y_edge + 1e-9, e)[0].real
    assert inside == pytest.approx(1.0, abs=1e-6) and outside == 0.0


def test_rooftop_charge_is_piecewise_constant_and_neutral():
    """div f = -/+ 1/dx on the two halves of an x-directed rooftop; zero net charge."""
    b = RooftopBasis(cell("rectangle", 4e-3, 4e-3), 8, 8)
    e = unit_vector(b, 0)
    xc, yc, h = b.center_x[0], b.center_y[0], 1e-7

    def jx(x):
        return b.evaluate(x, yc, e)[0].real

    for x in (xc - 0.8 * b.dx, xc - 0.3 * b.dx):
        assert (jx(x + h) - jx(x - h)) / (2 * h) == pytest.approx(+1 / b.dx, rel=1e-5)
    for x in (xc + 0.3 * b.dx, xc + 0.8 * b.dx):
        assert (jx(x + h) - jx(x - h)) / (2 * h) == pytest.approx(-1 / b.dx, rel=1e-5)
    # net charge: integral of div over the support = J(right) - J(left) = 0 - 0
    assert jx(xc - b.dx) == 0.0 and jx(xc + b.dx) == 0.0


# ----------------------------------------------------------- periodic wrap / Bloch

@pytest.mark.parametrize("kx0, ky0", [(0.0, 0.0), (137.0, -95.0)])
def test_bloch_periodicity_of_current(kx0, ky0):
    """J(r + P) = J(r) exp(-j k0.P) everywhere, including wrapped (seam) unknowns."""
    b = RooftopBasis(cell("slot", 6e-3, 2e-3), 8, 8)       # metal touches every cell edge
    rng = np.random.default_rng(4)
    coeffs = rng.normal(size=b.n_unknowns) + 1j * rng.normal(size=b.n_unknowns)
    x = rng.uniform(-PX, PX, 300)
    y = rng.uniform(-PY, PY, 300)
    base = b.evaluate(x, y, coeffs, kx0, ky0)
    for i, j in [(1, 0), (0, -1), (2, 1), (-1, 3)]:
        shifted = b.evaluate(x + i * PX, y + j * PY, coeffs, kx0, ky0)
        assert np.allclose(shifted, base * np.exp(-1j * (kx0 * i * PX + ky0 * j * PY)), atol=1e-10)


def test_current_flows_continuously_across_the_cell_boundary():
    """A wire spanning the period: J_x is continuous across x = +-Px/2 up to the Bloch phase."""
    b = RooftopBasis(cell("dipole", 8e-3, 1e-3), 8, 8)
    kx0 = 210.0
    coeffs = np.arange(1, b.n_unknowns + 1).astype(complex)
    y = np.full(5, b.center_y[0])
    eps = 1e-9
    tol = 10 * eps / b.dx * np.abs(coeffs).max()
    inside = b.evaluate(np.full(5, PX / 2 - eps), y, coeffs, kx0, 0.0)[:, 0]
    across = b.evaluate(np.full(5, PX / 2 + eps), y, coeffs, kx0, 0.0)[:, 0]
    assert np.allclose(inside, across, rtol=0, atol=tol)
    assert np.abs(inside).max() > 0.9                       # non-trivial: real current at the seam
    # Bloch relation on the far side: J(Px/2 + eps) = J(-Px/2 + eps) exp(-j kx0 Px)
    wrapped = b.evaluate(np.full(5, -PX / 2 + eps), y, coeffs, kx0, 0.0)[:, 0] * np.exp(-1j * kx0 * PX)
    assert np.allclose(across, wrapped, rtol=0, atol=tol)


@pytest.mark.parametrize("kind", ["patch", "dipole"])
def test_partition_of_unity_reproduces_uniform_current(kind):
    """Full-cell metal: all-ones coefficients reproduce J_x = 1 (or J_y = 1) exactly everywhere."""
    b = RooftopBasis(cell("rectangle", 8e-3, 8e-3) if kind == "patch" else cell("patch", 8e-3, 8e-3), 8, 8)
    rng = np.random.default_rng(5)
    x = rng.uniform(-PX, PX, 500)
    y = rng.uniform(-PY, PY, 500)
    ones_x = np.concatenate([np.ones(b.n_x), np.zeros(b.n_y)])
    ones_y = np.concatenate([np.zeros(b.n_x), np.ones(b.n_y)])
    Jx = b.evaluate(x, y, ones_x)
    Jy = b.evaluate(x, y, ones_y)
    assert np.allclose(Jx[:, 0], 1.0, atol=1e-12) and np.allclose(Jx[:, 1], 0.0, atol=1e-12)
    assert np.allclose(Jy[:, 1], 1.0, atol=1e-12) and np.allclose(Jy[:, 0], 0.0, atol=1e-12)
    assert (b.n_x, b.n_y) == (64, 64)


def test_zero_normal_current_on_metal_boundary_for_any_coefficients():
    """PEC boundary: normal J vanishes on the edge of an isolated patch (structural)."""
    c = cell("rectangle", 4e-3, 2e-3)
    b = RooftopBasis(c, 8, 8)
    rng = np.random.default_rng(6)
    coeffs = rng.normal(size=b.n_unknowns) + 1j * rng.normal(size=b.n_unknowns)
    t = np.linspace(-1.999e-3, 1.999e-3, 30)
    # x = +-2 mm edges (normal = x): J_x = 0 there
    for xe in (-2e-3, 2e-3):
        assert np.allclose(b.evaluate(np.full(30, xe), t / 2, coeffs)[:, 0], 0.0, atol=1e-12)
    # y = +-1 mm edges (normal = y): J_y = 0 there
    for ye in (-1e-3, 1e-3):
        assert np.allclose(b.evaluate(t, np.full(30, ye), coeffs)[:, 1], 0.0, atol=1e-12)
    # and there is no current at all outside the metal
    out = b.evaluate(np.array([3e-3, -3e-3, 0.0]), np.array([0.0, 2e-3, 2.5e-3]), coeffs)
    assert np.all(out == 0.0)


# --------------------------------------------------------------------- transform

def numeric_transform(b, index, kx, ky, n=512):
    """Direct midpoint quadrature of integral f e^{+jk.r} dA over the cell (interior rooftops)."""
    x = -PX / 2 + (np.arange(n) + 0.5) * PX / n
    X, Y = np.meshgrid(x, x, indexing="ij")
    J = b.evaluate(X, Y, unit_vector(b, index))
    integrand = (J[..., 0] + J[..., 1]) * np.exp(1j * (kx * X + ky * Y))
    return integrand.sum() * (PX / n) ** 2


@pytest.mark.parametrize("which", ["x", "y"])
@pytest.mark.parametrize("kx, ky", [(0.0, 0.0), (400.0, 0.0), (0.0, 700.0), (1234.0, -2345.0), (5000.0, 4000.0)])
def test_fourier_transform_closed_form_matches_quadrature(which, kx, ky):
    b = RooftopBasis(cell("rectangle", 4e-3, 4e-3), 8, 8)
    k = 3 if which == "x" else b.n_x + 2
    exact = b.fourier_transform(kx, ky)[k]
    assert exact == pytest.approx(numeric_transform(b, k, kx, ky), rel=2e-4, abs=2e-4 * b.dx * b.dy)


def test_transform_at_k0_is_the_rooftop_area():
    b = RooftopBasis(cell("rectangle", 4e-3, 2e-3), 8, 8)
    phi = b.fourier_transform(0.0, 0.0)
    assert np.allclose(phi, b.dx * b.dy, rtol=1e-15)


def test_transform_has_sinc_zeros_at_mesh_aliases():
    """x-rooftop: double zero of the transform at kx = 2 pi m / dx (m != 0), i.e. Floquet order m = Nx*m'."""
    b = RooftopBasis(cell("rectangle", 4e-3, 2e-3), 8, 8)
    kx = 2 * np.pi * 8 * np.array([1, 2, 3]) / PX
    phi = b.fourier_transform(kx, np.zeros(3))[: b.n_x]
    assert np.allclose(np.abs(phi), 0.0, atol=1e-25)


@pytest.mark.parametrize("kx0, ky0", [(0.0, 0.0), (180.0, -90.0)])
def test_transform_of_bloch_extension_equals_unwrapped_transform(kx0, ky0):
    """The transform over ONE cell of the Bloch-periodic (wrapped) rooftop equals Phi(k_mn):
    validates fold + Bloch phases in `evaluate` against the closed form used in the assembly."""
    b = RooftopBasis(cell("slot", 6e-3, 2e-3), 8, 8)
    seam = int(np.flatnonzero(np.isclose(b.center_x[: b.n_x], -PX / 2))[0])        # wrapped x-rooftop
    seam_y = b.n_x + int(np.flatnonzero(np.isclose(b.center_y[b.n_x :], -PY / 2))[0])
    n = 512
    x = -PX / 2 + (np.arange(n) + 0.5) * PX / n
    X, Y = np.meshgrid(x, x, indexing="ij")
    for idx in (seam, seam_y):
        J = b.evaluate(X, Y, unit_vector(b, idx), kx0, ky0)
        for m, nn in [(0, 0), (1, 0), (-2, 3)]:
            kx, ky = kx0 + 2 * np.pi * m / PX, ky0 + 2 * np.pi * nn / PY
            num = ((J[..., 0] + J[..., 1]) * np.exp(1j * (kx * X + ky * Y))).sum() * (PX / n) ** 2
            assert b.fourier_transform(kx, ky)[idx] == pytest.approx(num, rel=3e-4, abs=3e-4 * b.dx * b.dy)


def test_current_coefficients_match_numeric_cell_average():
    b = RooftopBasis(cell("rectangle", 6e-3, 2e-3), 8, 8)
    rng = np.random.default_rng(8)
    coeffs = rng.normal(size=b.n_unknowns) + 1j * rng.normal(size=b.n_unknowns)
    n = 800
    x = -PX / 2 + (np.arange(n) + 0.5) * PX / n
    X, Y = np.meshgrid(x, x, indexing="ij")
    J = b.evaluate(X, Y, coeffs)
    for kx, ky in [(0.0, 0.0), (2 * np.pi / PX, 0.0), (0.0, -4 * np.pi / PY)]:
        num = np.stack([(J[..., c] * np.exp(1j * (kx * X + ky * Y))).mean() for c in (0, 1)])
        assert b.current_coefficients(coeffs, kx, ky) == pytest.approx(num, rel=2e-4, abs=1e-6)
    out = b.current_coefficients(coeffs, np.zeros((3, 2)), np.zeros((3, 2)))
    assert out.shape == (3, 2, 2)
