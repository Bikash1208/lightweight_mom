"""Phase 30 S1 (MODEL_ONLY): DIRECT MoM leaky-mode solve of a rotating-axis supercell -- the locality test.

The tensor/Floquet leakage model (`tensor_floquet.py`, `leakage_model.py`) replaces every cell of a supercell by the
admittance tensor of the *periodic* (uniform-neighbour) cell. This module solves the whole N-cell supercell as one periodic
MoM problem instead (grounded slab, `backing="pec"`), and finds the complex mode where the impedance matrix goes singular:

* mode "k": real frequency, complex Bloch wavenumber kx0 = beta - j*alpha  (alpha [Np/m] is the leakage + dielectric attenuation);
* mode "f": real kx0, complex frequency f = f_r + j*f_i, Q = f_r / (2 f_i) (what an HFSS eigenmode solve reports).

Nothing in the frozen solver is modified. The production `RooftopBasis.fourier_transform` and `spectral_dyad` cast k to float, so the
analytic continuation is done here, reusing the production closed forms (`LayeredGreenFunction._impedances_from_kz`, `RooftopBasis`
edge bookkeeping, `EllipticalPatch`):

* conj(Phi_i(k)) -> Phi_i(-k) (identical for real k; analytic for complex k); the rooftop amplitudes are even in k;
* radiating Floquet harmonics (Re kz^2 > 0) use the OUTGOING branch Re(kz) > 0 -- the leaky pole lives on the improper sheet;
  evanescent harmonics keep the decaying branch Im(kz) <= 0 (substrate kz: Im <= 0, the closed forms are even in it);
* Z_ij depends only on the grid offset of the two rooftops (block Toeplitz with Bloch phase), so Z is filled by one 2-D FFT per
  block instead of an (N x N x modes) product: bit-for-bit the same sum, O(N^2) memory/time. Verified against
  `mom_matrix.assemble_impedance_matrix` in `tests/test_leaky_supercell_mom.py`.

The mode is the zero of g(p) = 1 / (w^T Z(p)^-1 v) (analytic, simple zero at a simple singularity; one LU per evaluation), found by a
secant iteration in the complex plane.
"""

import math
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from ..basis_functions import RooftopBasis
from ..constants import C0
from ..geometry import UnitCell
from ..green_function import LayeredGreenFunction
from ..materials import Metal
from ..periodic_bc import RectangularLattice
from .elliptical_cells import EllipticalPatch

_KT_ZERO_RTOL = 1e-12


# --------------------------------------------------------------------------------------------- geometry
@dataclass(frozen=True)
class RotatingSlotSupercell:
    """`n_cells` identical super-elliptic patches (fixed outline) in a row along x, the centre slot of cell i turned by `angles[i]` [rad].
    Duck-types `geometry.Geometry` (like `EllipticalPatch`), centred on the supercell."""

    n_cells: int
    pitch: float
    major: float
    exponent: float
    slot_length: float
    slot_width: float
    angles: tuple

    kind = "patch_rotating_slot_supercell"
    polarity = "metal"

    def __post_init__(self):
        if len(self.angles) != self.n_cells:
            raise ValueError("need one slot angle per cell")

    @property
    def cells(self):
        n, p = self.n_cells, self.pitch
        return [
            EllipticalPatch(
                self.major,
                self.major,
                0.0,
                self.slot_length,
                self.slot_width,
                float(a),
                False,
                self.exponent,
                -(n * p) / 2 + (i + 0.5) * p,
                0.0,
            )
            for i, a in enumerate(self.angles)
        ]

    @property
    def feature_area(self):
        return sum(c.feature_area for c in self.cells)

    def corners(self):
        hx, hy = self.n_cells * self.pitch / 2, self.pitch / 2
        return np.array([[hx, hy], [-hx, hy], [-hx, -hy], [hx, -hy]])

    def in_feature(self, x, y):
        out = np.zeros(np.broadcast(np.asarray(x), np.asarray(y)).shape, dtype=bool)
        for c in self.cells:
            out |= c.in_feature(x, y)
        return out


def s1_angles(n_cells=7):
    """Slot angles of the S1 supercell: 180*i/N degrees (the STEP / HFSS construction), in radians."""
    return tuple(math.radians(180.0 * i / n_cells) for i in range(n_cells))


def build_supercell(
    angles, mesh_per_cell, substrate, pitch=3.2e-3, major=3.05e-3, exponent=4.0, slot_length=1.89e-3, slot_width=0.58e-3
):
    """(UnitCell, RooftopBasis) of the grounded N-cell supercell on an (N*m) x m rooftop mesh."""
    n = len(angles)
    geo = RotatingSlotSupercell(n, pitch, major, exponent, slot_length, slot_width, tuple(angles))
    cell = UnitCell(RectangularLattice(n * pitch, pitch), geo, substrate, Metal(), backing="pec")
    return cell, RooftopBasis(cell, n * mesh_per_cell, mesh_per_cell)


# --------------------------------------------------------------------------------------------- green function
def _kz_air(z):
    """kz = sqrt(z), z = k^2 - kt^2: OUTGOING branch (Re kz >= 0) where Re z >= 0, decaying (Im kz <= 0) where Re z < 0."""
    z = np.asarray(z, dtype=complex)
    return np.where(z.real < 0.0, -1j * np.sqrt(-z), np.sqrt(z))


def _kz_proper(z):
    """kz = sqrt(z) with Im kz <= 0 (the substrate; the closed forms are even in it, this only keeps exponentials bounded)."""
    r = np.sqrt(np.asarray(z, dtype=complex))
    return np.where(r.imag > 0.0, -r, r)


def green_dyad_complex(substrate, f_complex, kx, ky):
    """Spectral dyad (gxx, gxy, gyy) [ohm] of air | pattern | substrate | PEC at complex frequency / wavenumber (analytic continuation)."""
    omega = 2.0 * math.pi * f_complex
    k = omega / C0
    eps = substrate.eps_complex
    kx = np.asarray(kx, dtype=complex)
    ky = np.asarray(ky, dtype=complex)
    kt2 = kx**2 + ky**2
    a = _kz_air(k * k - kt2)
    b = _kz_proper(eps * k * k - kt2)
    shim = SimpleNamespace(omega=omega, substrate=substrate, backing="pec")
    z_te, z_tm = LayeredGreenFunction._impedances_from_kz(shim, a, b)
    zero = np.abs(kt2) <= (_KT_ZERO_RTOL * abs(k)) ** 2
    s = np.where(zero, 1.0, kt2)
    avg = 0.5 * (z_te + z_tm)
    gxx = np.where(zero, -avg, -(z_tm * kx**2 + z_te * ky**2) / s)
    gyy = np.where(zero, -avg, -(z_tm * ky**2 + z_te * kx**2) / s)
    gxy = np.where(zero, 0.0, -(z_tm - z_te) * kx * ky / s)
    return gxx, gxy, gyy


def _sinc(z):
    """sin(z)/z for complex z."""
    z = np.asarray(z, dtype=complex)
    small = np.abs(z) < 1e-8
    zs = np.where(small, 1.0, z)
    return np.where(small, 1.0 - z * z / 6.0, np.sin(zs) / zs)


# --------------------------------------------------------------------------------------------- impedance matrix
def _positions(basis):
    """Per unknown: type (0 = x, 1 = y), integer grid index (A, B) and half-cell flags (hx, hy) of the rooftop centre."""
    ni = basis.edge_index
    d = basis.direction
    a, b = ni[:, 0].astype(np.int64), ni[:, 1].astype(np.int64)
    hx = np.where(d == 1, 0.5, 0.0)  # y-rooftops sit half a cell along x
    hy = np.where(d == 0, 0.5, 0.0)
    return d, a, b, hx, hy


def toeplitz_tables(basis, substrate, f_complex, kx0, ky0, order_x, order_y):
    """S[(di, dj)] (nx x ny complex): the block-Toeplitz generators of Z (units of Z, before the Bloch phases)."""
    lat = basis.lattice
    nx, ny, dx, dy = basis.nx, basis.ny, basis.dx, basis.dy
    m = np.arange(-order_x, order_x + 1)
    n = np.arange(-order_y, order_y + 1)
    mm, nn = np.meshgrid(m, n, indexing="ij")
    kx = kx0 + lat.reciprocal_x * mm
    ky = ky0 + lat.reciprocal_y * nn
    gxx, gxy, gyy = green_dyad_complex(substrate, f_complex, kx, ky)
    sx, sy = _sinc(kx * dx / 2.0), _sinc(ky * dy / 2.0)
    amp = ((dx * dy) * sx**2 * sy, (dx * dy) * sx * sy**2)  # x-, y-rooftop transforms (even in k)
    g = {(0, 0): gxx, (0, 1): gxy, (1, 0): gxy, (1, 1): gyy}
    q = (mm % nx) * ny + (nn % ny)
    scale = -(nx * ny) / lat.cell_area
    tables = {}
    for (di, dj), gg in g.items():
        shift_x = 0.5 * (dj == 1) - 0.5 * (di == 1)
        shift_y = 0.5 * (dj == 0) - 0.5 * (di == 0)
        w = amp[di] * amp[dj] * gg * np.exp(2j * np.pi * (mm * shift_x / nx + nn * shift_y / ny))
        fr = np.bincount(q.ravel(), weights=w.real.ravel(), minlength=nx * ny)
        fi = np.bincount(q.ravel(), weights=w.imag.ravel(), minlength=nx * ny)
        tables[(di, dj)] = scale * np.fft.ifft2((fr + 1j * fi).reshape(nx, ny))
    return tables


def assemble_z_complex(basis, substrate, f_complex, kx0, ky0, order_x, order_y, chunk=2048):
    """Z (n x n complex) at complex f and/or complex kx0 (ky0 real), Floquet orders |m| <= order_x, |n| <= order_y.
    Equals `mom_matrix.assemble_impedance_matrix` for real inputs and order_x == order_y."""
    lat = basis.lattice
    tables = toeplitz_tables(basis, substrate, f_complex, kx0, ky0, order_x, order_y)
    d, a, b, hx, hy = _positions(basis)
    psi = (a + hx) * basis.dx * kx0 + (b + hy) * basis.dy * ky0
    ph = np.exp(1j * psi)
    n = basis.n_unknowns
    z = np.empty((n, n), dtype=complex)
    nx, ny = basis.nx, basis.ny
    for di in (0, 1):
        ri = np.nonzero(d == di)[0]
        for dj in (0, 1):
            cj = np.nonzero(d == dj)[0]
            t = tables[(di, dj)]
            for s in range(0, ri.size, chunk):
                rows = ri[s : s + chunk]
                da = (a[cj][None, :] - a[rows][:, None]) % nx
                db = (b[cj][None, :] - b[rows][:, None]) % ny
                z[np.ix_(rows, cj)] = t[da, db] * (ph[cj][None, :] / ph[rows][:, None])
    del lat
    return z


# --------------------------------------------------------------------------------------------- mode search
def probe_function(z, v, w):
    """w^T Z^-1 v  (one LU)."""
    lu = lu_factor(z, overwrite_a=True, check_finite=False)
    return complex(w @ lu_solve(lu, v, check_finite=False))


def probe_vectors(n, seed=1):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    return v, v.copy()


def find_mode(zfun, p0, p1, n_unknowns, tol=1e-9, max_iter=25, seed=1, log=None, probe=None):
    """Secant search in the complex plane for the zero of h(p) = 1/(w^T Z(p)^-1 v); `zfun(p)` returns Z. Returns (p, history)."""
    v, w = probe if probe is not None else probe_vectors(n_unknowns, seed)

    def h(p):
        return 1.0 / probe_function(zfun(p), v, w)

    hist = []
    p0, p1 = complex(p0), complex(p1)
    h0, h1 = h(p0), h(p1)
    hist += [(p0, h0), (p1, h1)]
    for it in range(max_iter):
        if log:
            log(f"  it {it}: p = {p1!r}  |h| = {abs(h1):.3e}")
        if h1 == h0:
            break
        p2 = p1 - h1 * (p1 - p0) / (h1 - h0)
        step = abs(p2 - p1)
        p0, h0, p1 = p1, h1, p2
        h1 = h(p1)
        hist.append((p1, h1))
        if step <= tol * max(abs(p1), 1e-30):
            break
    return p1, hist


def harmonic_probe(basis, kx0, order=1):
    """Probe vector that overlaps the surface-wave current: x-directed rooftops carrying the space harmonic kx0 + 2*pi*order/L.
    Returns (v, w) with w = conj(v), for `find_mode(..., probe=...)`. A random probe has ~1/sqrt(N) overlap with the mode and hides
    a Q ~ 10^2-10^3 pole in a frequency scan; this one does not."""
    kh = kx0 + 2.0 * math.pi * order / basis.lattice.period_x
    v = (basis.direction == 0) * np.exp(1j * kh * basis.center_x)
    return v, np.conj(v)


def nearest_eigenvalue(z, x0=None, steps=4, seed=1):
    """Eigenvalue of Z nearest zero and its vector: inverse iteration (one LU), Rayleigh quotient x^T Z x / x^T x (exact once x is the eigenvector)
    (analytic in the parameters, second-order accurate in the vector). `x0` warm-starts it (default: the probe-vector solution)."""
    n = z.shape[0]
    if x0 is None:
        x0, _ = probe_vectors(n, seed)
    zc = z.copy()
    lu = lu_factor(zc, overwrite_a=True, check_finite=False)
    x = x0 / np.linalg.norm(x0)
    for _ in range(steps):
        x = lu_solve(lu, x, check_finite=False)
        x /= np.linalg.norm(x)
    return complex(x @ (z @ x) / (x @ x)), x


def find_mode_eig(zfun, p0, p1, tol=1e-10, max_iter=30, log=None):
    """Secant search for the zero of the eigenvalue of Z(p) nearest zero (better conditioned than the probe form when the start is not
    already close). Returns (p, eigenvalue-history, eigenvector)."""
    lam0, x = nearest_eigenvalue(zfun(p0))
    lam1, x = nearest_eigenvalue(zfun(p1), x)
    p0, p1 = complex(p0), complex(p1)
    hist = [(p0, lam0), (p1, lam1)]
    for it in range(max_iter):
        if log:
            log(f"  it {it}: p = {p1!r}  |lam| = {abs(lam1):.3e}")
        if lam1 == lam0:
            break
        p2 = p1 - lam1 * (p1 - p0) / (lam1 - lam0)
        step = abs(p2 - p1)
        p0, lam0, p1 = p1, lam1, p2
        lam1, x = nearest_eigenvalue(zfun(p1), x)
        hist.append((p1, lam1))
        if step <= tol * max(abs(p1), 1e-30):
            break
    return p1, hist, x


def orders_for(basis, cells_per_period, factor=2):
    """Floquet orders: `factor` x (mesh cells per physical cell) harmonics of the single-cell lattice, i.e. N_cells times as many in the
    supercell lattice along x (the repo's F = 2n rule, applied per cell)."""
    m = basis.ny
    return factor * m * cells_per_period, factor * m
