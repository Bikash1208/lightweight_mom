"""Floquet (harmonic-balance) solution of a periodically MODULATED TENSOR sheet -- the tensor twin of
`leakage_model.floquet_complex_beta` (Phase 26 B2.4 correction).

The Floquet constant of a modulated surface is not the mean of the local constants (a second-order shift in
the modulation depth, plus coupling into the radiating n = -1 harmonic). Here the sheet is a staircase of
cells, cell `i` carrying the 2 x 2 tensor `Y_i` (rows/cols = x, y components) over `[i P, (i+1) P)`; the
period is `N P` along x. Harmonic `n` has transverse wavevector `(beta cos g + n K, beta sin g)` and sees the
same up-looking-air + down-looking-grounded-slab dyad as the scalar balance (`leakage_model._harmonic_block`,
reused unmodified); the sheet couples harmonics by the Fourier coefficients of `Y(x)`, which are 2 x 2.
Approximation: the cell tensors are those of ONE transverse wavenumber (the database/Bloch tensor), applied
to every harmonic -- spatial dispersion between harmonics is not modelled.

Labels: SOLVER_INTERNAL (a model of the model; not an HFSS number).
"""

import numpy as np

from src.dispersion.leakage_model import _balance_matrix, _harmonic_block, bound_beta  # noqa: F401
from src.frequency import angular_frequency, free_space_wavenumber


def staircase_fourier_tensor(tensors, n_max):
    """`Y_p` (shape `(4 n_max + 1, 2, 2)`, `p = -2 n_max .. 2 n_max`) of a staircase of `N` equal-width cells
    carrying `tensors` (N, 2, 2): `Y(x) = sum_p Y_p exp(-j p K x)`, `K = 2 pi / (N P)`."""
    t = np.asarray(tensors, dtype=complex)
    n = t.shape[0]
    p = np.arange(-2 * n_max, 2 * n_max + 1)
    out = np.zeros((p.size, 2, 2), dtype=complex)
    for ip, k in enumerate(p):
        if k == 0:
            out[ip] = t.mean(axis=0)
            continue
        # (1/L) int_{iP}^{(i+1)P} exp(+j k K x) dx, K P = 2 pi / N
        a = 2.0 * np.pi * k / n
        w = (np.exp(1j * a * (np.arange(n) + 1)) - np.exp(1j * a * np.arange(n))) / (1j * a * n)
        out[ip] = np.tensordot(w, t, axes=(0, 0))
    return out


def _balance_matrix_tensor(beta, y_coeffs, n_max, k_mod, gamma_rad, substrate, frequency_hz):
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    size = 2 * n_max + 1
    m = np.zeros((2 * size, 2 * size), dtype=complex)
    for i, n in enumerate(range(-n_max, n_max + 1)):
        for j, q in enumerate(range(-n_max, n_max + 1)):
            m[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = y_coeffs[(n - q) + 2 * n_max]
        m[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] += _harmonic_block(
            beta, n, k_mod, gamma_rad, substrate.eps_complex, substrate.thickness, omega, k0
        )
    return m


def floquet_complex_beta_tensor(tensors, period_m, substrate, frequency_hz, beta_start, n_max=6, steps=8, tol=1e-10):
    """Complex `beta = beta' - j alpha` of the TM0-like mode travelling along +x on the staircase of
    `tensors` (period `period_m` = N P). `beta_start`: the real constant of the unmodulated MEAN tensor
    (modulation switched on in `steps` equal steps: the root connected to the bound mode is the one found).
    Returns `(beta, n_iterations_last)`."""
    k_mod = 2.0 * np.pi / period_m
    y = staircase_fourier_tensor(tensors, n_max)
    centre = 2 * n_max
    beta = complex(beta_start)
    its = 0
    for s in np.linspace(1.0 / steps, 1.0, steps):
        scaled = y * s
        scaled[centre] = y[centre]

        def det(b, scaled=scaled):
            m = _balance_matrix_tensor(b, scaled, n_max, k_mod, 0.0, substrate, frequency_hz)
            sign, logdet = np.linalg.slogdet(m)
            return sign * np.exp(logdet)

        b0, b1 = beta, beta * (1.0 + 1e-5) + 1e-3j
        f0, f1 = det(b0), det(b1)
        for its in range(1, 80):
            if f1 == f0:
                break
            b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
            b0, f0, b1 = b1, f1, b2
            f1 = det(b1)
            if abs(b1 - b0) < tol * abs(b1):
                break
        beta = b1
    return beta, its


def mean_tensor_beta(tensors, substrate, frequency_hz, n_scan=2000):
    """Real constant of the area-averaged tensor along x (effective-medium reference)."""
    from src.dispersion.tensor_sheet import find_beta_tensor

    ym = np.mean(np.asarray(tensors, dtype=complex), axis=0)
    return find_beta_tensor(ym, 0.0, substrate, frequency_hz, n_scan=n_scan)[0]


# ---------------------------------------------------------------------------------- harmonic-dependent tensors
def floquet_complex_beta_dispersive(y_by_m, period_m, substrate, frequency_hz, beta_start, n_max, steps=8, tol=1e-10):
    """As `floquet_complex_beta_tensor`, but the coupling of harmonics `n <- q` uses the tensor staircase
    evaluated at the MEAN transverse wavenumber `beta + (n + q) K / 2` (symmetric, reciprocal-preserving).
    `y_by_m[s]`, `s = n + q = -2 n_max .. 2 n_max`: Fourier coefficients (`staircase_fourier_tensor`, shape
    `(4 n_max + 1, 2, 2)`) of the cells' tensors at `kx = beta_ref + s K / 2`."""
    k_mod = 2.0 * np.pi / period_m
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    size = 2 * n_max + 1
    beta = complex(beta_start)
    its = 0
    for s_scale in np.linspace(1.0 / steps, 1.0, steps):

        def det(b, s_scale=s_scale):
            m = np.zeros((2 * size, 2 * size), dtype=complex)
            for i, n in enumerate(range(-n_max, n_max + 1)):
                for j, q in enumerate(range(-n_max, n_max + 1)):
                    c = np.array(y_by_m[n + q][(n - q) + 2 * n_max])
                    if n != q:
                        c = c * s_scale
                    m[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = c
                m[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] += _harmonic_block(
                    b, n, k_mod, 0.0, substrate.eps_complex, substrate.thickness, omega, k0
                )
            sign, logdet = np.linalg.slogdet(m)
            return sign * np.exp(logdet)

        b0, b1 = beta, beta * (1.0 + 1e-5) + 1e-3j
        f0, f1 = det(b0), det(b1)
        for its in range(1, 80):
            if f1 == f0:
                break
            b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
            b0, f0, b1 = b1, f1, b2
            f1 = det(b1)
            if abs(b1 - b0) < tol * abs(b1):
                break
        beta = b1
    return beta, its
