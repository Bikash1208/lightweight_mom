"""Impedance-matrix assembly: exact identities, symmetry/reciprocity, positivity, integration."""

import math

import numpy as np
import pytest

from src import mom_matrix
from src.basis_functions import RooftopBasis
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.materials import Metal, Substrate
from src.mom_matrix import assemble_impedance_matrix
from src.periodic_bc import RectangularLattice

F = 10e9
PX = PY = 8e-3
LAT = RectangularLattice(PX, PY)


def make(kind, length, width, tand=0.0009, eps_r=2.2, **kw):
    sub = Substrate("t", eps_r, tand, 1.575e-3)
    c = UnitCell(LAT, Geometry(kind, length, width, **kw), sub, Metal())
    return c, LayeredGreenFunction(sub, F)


class ConstantKernel:
    """Mock 'Green function' with a k-independent dyad, for exactly solvable identities."""

    def __init__(self, dyad):
        self.dyad = np.asarray(dyad, dtype=complex)

    def floquet_dyads(self, lattice, kx0, ky0, order):
        kx, ky = lattice.transverse_wavenumbers(kx0, ky0, *lattice.mode_grid(order))
        g = np.broadcast_to(self.dyad, kx.shape + (2, 2)).copy()
        return kx, ky, g


# ================================================= exactly solvable kernels (Parseval)

def rooftop_gram(b):
    """Analytic <f_i, f_j> for rooftops (vector inner product): 2/3 self, 1/6 neighbours."""
    n = b.n_unknowns
    gram = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if b.direction[i] != b.direction[j]:
                continue
            along, across = (0, 1) if b.direction[i] == 0 else (1, 0)
            d = [b.dx, b.dy]
            dt = abs([b.center_x[i] - b.center_x[j], b.center_y[i] - b.center_y[j]][along])
            ds = abs([b.center_x[i] - b.center_x[j], b.center_y[i] - b.center_y[j]][across])
            if ds > 1e-12:
                continue
            if dt < 1e-12:
                gram[i, j] = 2 * d[along] / 3 * d[across]
            elif abs(dt - d[along]) < 1e-12:
                gram[i, j] = d[along] / 6 * d[across]
    return gram


def test_identity_kernel_reproduces_minus_the_analytic_gram_matrix():
    """G = I  =>  Z = -(1/A) sum conj(Phi_i) Phi_j = -<f_i, f_j>  (Parseval). Validates the
    transform closed forms, the 1/A normalization, the x/y block placement and the index
    bookkeeping with NO electromagnetics involved. The Parseval sum converges like 1/M."""
    c, _ = make("patch", 4e-3, 4e-3)
    b = RooftopBasis(c, 8, 8)
    gram = rooftop_gram(b)
    errs = {}
    for order in (100, 200, 400):
        z = assemble_impedance_matrix(b, ConstantKernel(np.eye(2)), 0.0, 0.0, order)
        errs[order] = np.linalg.norm(z + gram) / np.linalg.norm(gram)
    assert errs[400] < 3e-3
    assert 1.6 < errs[100] / errs[200] < 2.5 and 1.6 < errs[200] / errs[400] < 2.5   # ~ 1/M tail
    z = assemble_impedance_matrix(b, ConstantKernel(np.eye(2)), 0.0, 0.0, 400)
    assert np.allclose(z.imag, 0.0, atol=1e-18)
    assert abs(z[0, 0] + gram[0, 0]) < 1e-2 * gram[0, 0]


def test_off_diagonal_kernel_couples_x_and_y_rooftops_exactly():
    """G = [[0,1],[1,0]] =>  Z_xy[i,j] = -integral of (x-rooftop i, scalar)(y-rooftop j, scalar).
    Overlap of tri x pulse with pulse x tri is dx/2 * dy/2 when the centres are offset by
    (dx/2, dy/2), else 0. Checks the G_xy entry and block wiring."""
    c, _ = make("patch", 4e-3, 4e-3)
    b = RooftopBasis(c, 8, 8)
    n_x = b.n_x
    expected = np.zeros((n_x, b.n_y))
    for i in range(n_x):
        for j in range(b.n_y):
            ox = abs(b.center_x[i] - b.center_x[n_x + j])
            oy = abs(b.center_y[i] - b.center_y[n_x + j])
            if abs(ox - b.dx / 2) < 1e-12 and abs(oy - b.dy / 2) < 1e-12:
                expected[i, j] = b.dx * b.dy / 4
    z = assemble_impedance_matrix(b, ConstantKernel([[0, 1], [1, 0]]), 0.0, 0.0, 400)
    assert np.linalg.norm(z[:n_x, n_x:] + expected) < 3e-3 * np.linalg.norm(expected)
    assert np.allclose(z[n_x:, :n_x], z[:n_x, n_x:].T, rtol=0, atol=1e-12 * np.abs(z).max())
    assert np.abs(z[:n_x, :n_x]).max() < 1e-20 and np.abs(z[n_x:, n_x:]).max() < 1e-20
    assert np.count_nonzero(expected) > 0


def test_assembly_equals_the_defining_double_sum():
    """Brute-force loop over harmonics with per-mode transforms == vectorized assembly."""
    c, g = make("rectangle", 4e-3, 2e-3)
    b = RooftopBasis(c, 8, 8)
    kx0, ky0, order = 37.0, -21.0, 2
    z = assemble_impedance_matrix(b, g, kx0, ky0, order)
    ref = np.zeros_like(z)
    for m in range(-order, order + 1):
        for n in range(-order, order + 1):
            kx = kx0 + 2 * np.pi * m / PX
            ky = ky0 + 2 * np.pi * n / PY
            dyad = g.spectral_dyad(kx, ky)
            phi = b.fourier_transform(kx, ky)
            for i in range(b.n_unknowns):
                for j in range(b.n_unknowns):
                    ref[i, j] += np.conj(phi[i]) * dyad[b.direction[i], b.direction[j]] * phi[j]
    ref *= -1.0 / (PX * PY)
    assert z == pytest.approx(ref, rel=1e-12, abs=1e-25)


def test_chunked_assembly_is_identical(monkeypatch):
    c, g = make("slot", 6e-3, 2e-3)
    b = RooftopBasis(c, 8, 8)
    full = assemble_impedance_matrix(b, g, 0.0, 0.0, 6)
    monkeypatch.setattr(mom_matrix, "_CHUNK_ELEMENTS", 500)      # forces many tiny chunks
    chunked = assemble_impedance_matrix(b, g, 0.0, 0.0, 6)
    assert chunked == pytest.approx(full, rel=1e-12, abs=1e-25)


# ============================================================ uniform-current identity

@pytest.mark.parametrize("order", [0, 2, 9])
@pytest.mark.parametrize("tand", [0.0, 0.0009])
def test_uniform_current_sees_only_the_k0_harmonic(order, tand):
    """Full-metal cell, all-ones x-coefficients = uniform J_x = 1: all harmonics except (0,0)
    vanish exactly (sinc zeros), so Z * 1_x = -dx dy G_xx(0) on x rows, 0 on y rows —
    independent of mesh and Floquet order."""
    c, g = make("patch", 8e-3, 8e-3, tand=tand)
    b = RooftopBasis(c, 8, 8)
    z = assemble_impedance_matrix(b, g, 0.0, 0.0, order)
    ones_x = np.concatenate([np.ones(b.n_x), np.zeros(b.n_y)])
    ones_y = 1.0 - ones_x
    g0 = g.spectral_dyad(0.0, 0.0)
    rows = z @ ones_x
    assert rows[: b.n_x] == pytest.approx(-b.dx * b.dy * g0[0, 0], rel=1e-12)
    assert np.abs(rows[b.n_x :]).max() < 1e-12 * abs(rows[0])
    rows = z @ ones_y
    assert rows[b.n_x :] == pytest.approx(-b.dx * b.dy * g0[1, 1], rel=1e-12)
    assert np.abs(rows[: b.n_x]).max() < 1e-12 * abs(rows[-1])


def test_free_space_sheet_impedance_is_eta0_over_2():
    """eps_r = 1: uniform current sees Z = eta0/2 (E = -eta0 J/2), so 1^T Z 1 / A' = eta0/2 * area."""
    from src.constants import ETA0

    c, g = make("patch", 8e-3, 8e-3, tand=0.0, eps_r=1.0)
    b = RooftopBasis(c, 8, 8)
    z = assemble_impedance_matrix(b, g, 0.0, 0.0, 3)
    ones = np.concatenate([np.ones(b.n_x), np.zeros(b.n_y)])
    total = ones @ z @ ones                    # = A * (eta0/2), A = cell area (all rooftops overlap)
    assert total == pytest.approx(PX * PY * ETA0 / 2, rel=1e-12)


# =========================================================== symmetry / reciprocity

CASES = [("rectangle", 6e-3, 2e-3), ("dipole", 8e-3, 1e-3), ("slot", 6e-3, 2e-3),
         ("patch", 4e-3, 4e-3)]


@pytest.mark.parametrize("kind, length, width", CASES)
@pytest.mark.parametrize("tand", [0.0, 0.0009, 0.05])
def test_normal_incidence_matrix_is_complex_symmetric(kind, length, width, tand):
    c, g = make(kind, length, width, tand=tand)
    b = RooftopBasis(c, 8, 8)
    z = assemble_impedance_matrix(b, g, 0.0, 0.0, 16)
    assert np.linalg.norm(z - z.T) <= 1e-12 * np.linalg.norm(z)
    if tand > 0:
        assert np.linalg.norm(z - z.conj().T) > 1e-6 * np.linalg.norm(z)   # NOT Hermitian: lossy


@pytest.mark.parametrize("kind, length, width", CASES[:3])
def test_oblique_reciprocity_Z_of_k0_transposed_equals_Z_of_minus_k0(kind, length, width):
    """Lorentz reciprocity for a Bloch problem: Z(k0)^T = Z(-k0), while Z(k0) alone is not symmetric."""
    c, g = make(kind, length, width)
    b = RooftopBasis(c, 8, 8)
    kx0, ky0 = 0.21 * g.k0, 0.12 * g.k0
    za = assemble_impedance_matrix(b, g, kx0, ky0, 16)
    zb = assemble_impedance_matrix(b, g, -kx0, -ky0, 16)
    assert np.linalg.norm(za.T - zb) <= 1e-12 * np.linalg.norm(za)
    assert np.linalg.norm(za - za.T) > 1e-4 * np.linalg.norm(za)


@pytest.mark.parametrize("kind, length, width", CASES)
def test_dissipation_operator_is_positive_semidefinite(kind, length, width):
    """Re(Z) = Hermitian part of a complex-symmetric Z: I^H Z I real part = 2 * power delivered >= 0."""
    c, g = make(kind, length, width)
    b = RooftopBasis(c, 8, 8)
    z = assemble_impedance_matrix(b, g, 0.0, 0.0, 16)
    herm = 0.5 * (z + z.conj().T)
    eig = np.linalg.eigvalsh(herm)
    assert eig.min() >= -1e-12 * eig.max()
    assert eig.min() > 0.0                                    # lossy substrate: strictly positive
    rng = np.random.default_rng(1)
    for _ in range(20):
        i = rng.normal(size=b.n_unknowns) + 1j * rng.normal(size=b.n_unknowns)
        assert np.real(np.conj(i) @ z @ i) > 0.0


def test_lossless_dissipation_operator_has_rank_two():
    """Lossless, one propagating order: only (0,0) radiates, into two polarizations => rank(Re Z) = 2."""
    c, g = make("rectangle", 6e-3, 2e-3, tand=0.0)
    b = RooftopBasis(c, 16, 16)
    z = assemble_impedance_matrix(b, g, 0.0, 0.0, 32)
    eig = np.linalg.eigvalsh(0.5 * (z + z.conj().T))
    assert eig.min() >= -1e-12 * eig.max()
    assert np.sum(eig > 1e-9 * eig.max()) == 2
    # reactive part is not small
    assert np.linalg.norm(z.imag) > 1e-3 * np.linalg.norm(z.real)


# ================================================== Green-function integration

def spatial_domain_entries(b, g, rows, cols, n):
    """Independent evaluation (normal incidence): sample J on an n x n grid, obtain E by a direct
    DFT of the sampled current (Floquet sum up to Nyquist), and integrate f_i . E on the grid.
    Uses neither the closed-form transforms nor `assemble_impedance_matrix`."""
    px, py = b.lattice.period_x, b.lattice.period_y
    xs = -px / 2 + np.arange(n) * px / n
    ys = -py / 2 + np.arange(n) * py / n
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    m = np.arange(-n // 2, n // 2)
    kx, ky = 2 * np.pi * m / px, 2 * np.pi * m / py
    ex, ey = np.exp(1j * np.outer(kx, xs)), np.exp(1j * np.outer(ky, ys))
    kxx, kyy = np.meshgrid(kx, ky, indexing="ij")
    dyad = g.spectral_dyad(kxx, kyy)
    tests = {}
    for i in rows:
        e_i = np.zeros(b.n_unknowns, complex)
        e_i[i] = 1.0
        tests[i] = b.evaluate(X, Y, e_i)
    out = np.zeros((len(rows), len(cols)), complex)
    for jj, j in enumerate(cols):
        e_j = np.zeros(b.n_unknowns, complex)
        e_j[j] = 1.0
        jf = b.evaluate(X, Y, e_j)
        jt = np.stack([ex @ jf[..., c] @ ey.T / n**2 for c in (0, 1)], axis=-1)   # J~ = mean(J e^{+jk.r})
        et = np.einsum("mnij,mnj->mni", dyad, jt)
        field = np.stack([ex.conj().T @ et[..., c] @ ey.conj() for c in (0, 1)], axis=-1)
        for ii, i in enumerate(rows):
            out[ii, jj] = -np.sum(tests[i] * field) * px * py / n**2
    return out


@pytest.mark.parametrize("kind, length, width", [("rectangle", 4e-3, 2e-3), ("dipole", 8e-3, 1e-3), ("slot", 6e-3, 2e-3)])
def test_assembly_agrees_with_independent_spatial_domain_evaluation(kind, length, width):
    """The physical Green's function, including G_xy coupling, wrapped (seam) rooftops and the
    metal-complement slot, gives the same Z entries whether the Floquet sum is done with the
    closed-form rooftop transforms or by DFT of the sampled current. The residual difference is
    the sampling error of the *check* (rooftops sampled on the DFT grid): it must shrink with n."""
    c, g = make(kind, length, width)
    b = RooftopBasis(c, 8, 8)
    idx = sorted({0, b.n_x // 2, b.n_x - 1, min(b.n_x, b.n_unknowns - 1), b.n_unknowns - 1})
    diffs = []
    for n in (64, 128):
        zs = spatial_domain_entries(b, g, idx, idx, n)
        za = assemble_impedance_matrix(b, g, 0.0, 0.0, n // 2)[np.ix_(idx, idx)]
        diffs.append(np.linalg.norm(zs - za) / np.linalg.norm(za))
    assert diffs[1] < 0.6 * diffs[0]
    assert diffs[1] < 0.06


# ================================================================ Floquet truncation

def test_matrix_entries_converge_with_floquet_order():
    """Cauchy differences of Z entries shrink ~ 1/M^2 (as in the Green-function phase)."""
    c, g = make("rectangle", 4e-3, 2e-3)
    b = RooftopBasis(c, 8, 8)
    orders = [8, 16, 32, 64, 128]
    zs = [assemble_impedance_matrix(b, g, 0.0, 0.0, m) for m in orders]
    diffs = np.array([np.linalg.norm(zs[i + 1] - zs[i]) for i in range(len(orders) - 1)])
    assert np.all(np.diff(diffs) < 0.0)
    ratios = diffs[:-1] / diffs[1:]
    assert np.all(ratios > 2.5)
    assert diffs[-1] < 2e-2 * np.linalg.norm(zs[-1])


def test_reactive_and_resistive_parts_are_separately_meaningful_lossless():
    """Lossless: Re Z is independent of Floquet order once (0,0) is included; Im Z converges."""
    c, g = make("rectangle", 4e-3, 2e-3, tand=0.0)
    b = RooftopBasis(c, 8, 8)
    z4 = assemble_impedance_matrix(b, g, 0.0, 0.0, 4)
    z32 = assemble_impedance_matrix(b, g, 0.0, 0.0, 32)
    assert np.allclose(z4.real, z32.real, rtol=0, atol=1e-12 * np.abs(z32.real).max())
    assert np.linalg.norm(z32.imag - z4.imag) > 1e-3 * np.linalg.norm(z32.imag)


def test_input_validation():
    c, g = make("rectangle", 4e-3, 2e-3)
    b = RooftopBasis(c, 8, 8)
    with pytest.raises(ValueError):
        mom_matrix.assemble_excitation_vector(b, np.ones(3), 0.0, 0.0)
    with pytest.raises(ValueError):
        assemble_impedance_matrix(b, g, 0.0, 0.0, -1)
