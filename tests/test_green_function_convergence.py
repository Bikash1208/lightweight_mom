"""Convergence of Floquet-harmonic sums built from the spectral dyad.

The point-source Floquet sum of G diverges at the source (physical singularity),
so convergence is tested on the quantity a Galerkin MoM will actually need, a
bilinear form with a smooth test current f:

    Z[f] = -PxPy * sum_{m,n} conj(J_mn) . G(kxm, kyn) . J_mn ,
    J_mn = (1/(PxPy)) * F(kxm, kyn),   F(kx, ky) = int f exp(+j k.r) dA,

for f = rooftop(x; half-width W) * pulse(y; length L) directed along x, whose
transform is known in closed form:

    F(kx, ky) = W sinc^2(kx W / 2pi) * L sinc(ky L / 2pi)     (np.sinc convention).

(This is a stand-in current to test summation behaviour; it is NOT a MoM matrix.)
"""

import numpy as np
import pytest

from src.green_function import LayeredGreenFunction
from src.materials import Substrate
from src.periodic_bc import RectangularLattice

PX = PY = 8e-3
LAT = RectangularLattice(PX, PY)
W, L = 2e-3, 3e-3
ORDERS = [4, 8, 16, 32, 64, 128, 256]


def transform(kx, ky):
    return W * np.sinc(kx * W / (2 * np.pi)) ** 2 * L * np.sinc(ky * L / (2 * np.pi))


def bilinear_form(g, N, kx0=0.0, ky0=0.0, direction="x"):
    kx, ky, G = g.floquet_dyads(LAT, kx0, ky0, N)
    J = np.zeros(kx.shape + (2,), dtype=complex)
    if direction == "x":
        J[..., 0] = transform(kx, ky) / (PX * PY)
    else:  # same shape rotated 90 degrees
        J[..., 1] = transform(ky, kx) / (PX * PY)
    return -PX * PY * np.einsum("mni,mnij,mnj->", J.conj(), G, J)


def series(g, **kw):
    return np.array([bilinear_form(g, N, **kw) for N in ORDERS])


def gf(tand):
    return LayeredGreenFunction(Substrate("s", 2.2, tand, 1.575e-3), 10e9)


@pytest.mark.parametrize("tand", [0.0, 0.0009])
def test_partial_sums_converge_monotonically_with_floquet_order(tand):
    q = series(gf(tand))
    diffs = np.abs(np.diff(q))
    assert np.all(np.diff(diffs) < 0.0)  # Cauchy differences strictly shrink
    assert diffs[-1] < 1e-3 * abs(q[-1])  # and are small in relative terms
    # measured tail ~ N^-2 (ratio -> 4); demand clearly better than 1/N
    ratios = diffs[:-1] / diffs[1:]
    assert np.all(ratios[1:] > 3.0)
    assert ratios[-1] == pytest.approx(4.0, rel=0.05)


def test_richardson_extrapolation_confirms_algebraic_order_two():
    q = series(gf(0.0009))
    ext = q[1:] + (q[1:] - q[:-1]) / 3.0  # eliminates the c/N^2 term
    # extrapolated values agree far better than the raw partial sums
    assert abs(ext[-1] - ext[-2]) < 0.02 * abs(q[-1] - q[-2])
    assert abs(q[-1] - ext[-1]) / abs(ext[-1]) < 2e-4


def test_lossless_real_part_is_exact_after_propagating_orders():
    """Lossless: only the (0,0) harmonic radiates (Px < lambda0/2), so Re Z[f] is
    independent of truncation and all convergence is in the reactive part."""
    q = series(gf(0.0))
    assert np.allclose(q.real, q.real[0], rtol=1e-12, atol=0.0)
    assert q.real[0] > 0.0
    assert abs(q[-1].imag - q[0].imag) > 0.0


def test_lossy_real_part_converges_and_stays_positive():
    q = series(gf(0.0009))
    assert np.all(q.real > 0.0)
    assert abs(q.real[-1] - q.real[-2]) < 1e-5 * q.real[-1]


def test_converges_for_oblique_floquet_phase():
    g = gf(0.0009)
    q = series(g, kx0=0.35 * LAT.reciprocal_x, ky0=-0.2 * LAT.reciprocal_y)
    diffs = np.abs(np.diff(q))
    assert np.all(np.diff(diffs) < 0.0) and diffs[-1] < 1e-3 * abs(q[-1])


def test_converges_for_orthogonal_current_direction_and_agrees_by_symmetry():
    """Square lattice, isotropic stack: a y-directed rotated copy of f gives the same Z[f]."""
    g = gf(0.0009)
    qx = series(g, direction="x")
    qy = series(g, direction="y")
    assert qy == pytest.approx(qx, rel=1e-12)


def test_kernel_tail_decay_rates():
    """Asymptotics that make the sums converge: Z_TE ~ 1/kt, Z_TM ~ kt (both finite prefactors)."""
    g = gf(0.0009)
    kt = np.array([1e3, 1e4, 1e5]) * g.k0
    z_te, z_tm = g.equivalent_impedances(kt, 0.0)
    assert np.allclose(np.abs(z_te) * kt, np.abs(z_te[0]) * kt[0], rtol=1e-3)
    assert np.allclose(np.abs(z_tm) / kt, np.abs(z_tm[0]) / kt[0], rtol=1e-3)
