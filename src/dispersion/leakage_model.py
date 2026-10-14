"""Phase 26 A: how fast a TM surface wave on a grounded slab, loaded by a periodically modulated patch
sheet, leaks into the radiating n = -1 Floquet harmonic -- the leakage constant `alpha` [Np/m].

Two independent evaluations of the SAME physics, so each checks the other:

* `leakage_alpha_first_order`: the closed form (Oliner-Hessel / Minatti adiabatic Floquet, first order in
  the modulation). Only the patch SHEET is modulated and it is local; the grounded slab is not (its
  admittance depends on the transverse wavenumber), so the radiated harmonic is driven by a sheet current
  and responds with ITS OWN slab admittance, not the surface wave's:

      E_0      = (eta0 kappa / k0) H_0           tangential E of the surface wave at the patch plane,
                                                 kappa = sqrt(beta^2 - k0^2)
      J_{-1}   = c1 E_0                          c1: amplitude of the n = -1 harmonic of Y_sheet [S]
      Y_c      = Y_up,c + Y_slab,c(k_t) + Y_sheet,mean      c = TM, TE at k_t = k0 sin(theta)
      P_rad    = 1/2 |J_{-1}|^2 [ cos^2(zeta) Re(Y_up,TM)/|Y_TM|^2 + sin^2(zeta) Re(Y_up,TE)/|Y_TE|^2 ]
      P_sw     = |H_0|^2 w                       surface-wave power per unit width (air + substrate,
                                                 `surface_wave_power_per_h2`)
      alpha    = P_rad / (2 P_sw)                (dP_sw/dx = -P_rad = -2 alpha P_sw)

  `zeta` is the in-plane angle between the propagation direction and the radiated harmonic, so the field
  along the propagation direction splits into a TM part (cos^2) and a TE part (sin^2). This is a
  sheet-modulation definition of leakage (the grounded slab is part of the background, not of the
  modulated sheet) -- the penetrable-sheet picture of Martini et al., EPJ AM 2020, not an impenetrable
  impedance surface.

* `floquet_complex_beta`: the exact (truncated Floquet harmonic balance, TM + TE coupled, every harmonic with
  its own grounded-slab admittance) complex propagation constant `beta - j alpha` of the periodic sheet.
  It contains every effect the closed form drops (finite modulation depth, binary rather than sinusoidal
  profile, the other harmonics, the n = 0 wave's own detuning). Unlike the closed form it assumes the
  pattern is locally a straight-stripe grating (the usual adiabatic assumption), exp(+j omega t).

Everything here is a MODEL of local leakage: a surface wave's `alpha` along a ray assumes the modulation
varies slowly on the scale of a surface-wavelength and that the surface wave is a plane TM wave locally
(`rho^-1/2` spreading is accounted for outside, in `hologram.leaky_forward`). It is model-only until HFSS.
"""

import numpy as np

from ..constants import EPS0, ETA0, MU0
from ..frequency import angular_frequency, free_space_wavenumber
from ..labels import EXPERIMENTAL
from .sheet_model import _y_air, _y_slab, find_beta_from_sheet_admittance

LEAKAGE_MODEL_LABELS = (EXPERIMENTAL, "LEAKAGE_FIRST_ORDER_ADIABATIC", "MODEL_ONLY")


# ---------------------------------------------------------------------------------- the surface wave
def surface_wave_power_per_h2(beta, eps_r, thickness_m, frequency_hz):
    """`w` [ohm m]: the TM0 surface-wave power flow per unit width, per `|H_0|^2` with `H_0` the magnetic
    field amplitude just above the patch plane: `P_sw = w |H_0|^2`.

    Uses the transverse-resonance identity `P_sw = |dB/dbeta| |E_0|^2 / 4` of a lossless reactive
    structure, `B(beta) = Im[Y_air(beta) + Y_slab(beta)]` the total TM susceptance at the patch plane (the
    patch layer's own susceptance does not depend on `beta`), `E_0 = H_0 kappa/(omega eps0)` the tangential
    E there. It is the air-side power `eta0 (beta/k0)/(4 kappa)` plus the slab's, evaluated with the H
    field at the slab's TOP (which differs from `H_0` by the sheet current -- integrating the slab with
    `H_0` instead, the obvious shortcut, is a few percent wrong; the test compares both)."""
    k0 = free_space_wavenumber(frequency_hz)
    omega = angular_frequency(frequency_hz)
    kappa = np.sqrt(beta**2 - k0**2)
    kz = np.sqrt(eps_r * k0**2 - beta**2)
    u = kz * thickness_m
    db_air = -omega * EPS0 * beta / kappa**3
    db_slab = -omega * EPS0 * eps_r * (beta / kz**3) * (1.0 / np.tan(u) + u / np.sin(u) ** 2)
    e0_per_h = ETA0 * kappa / k0
    return float(abs(db_air + db_slab) * e0_per_h**2 / 4.0)


def state_surface_reactances(y_sheets, beta, substrate, frequency_hz):
    """`X_s` [ohm] of each state: `Im(1 / (Y_sheet,s + Y_slab_TM(beta)))`, the impenetrable-surface reactance
    a TM surface wave of constant `beta` sees at the patch plane (the slab's own TM admittance at `beta`
    in parallel with the patch layer's `Y_sheet`, step 2 of `sheet_model`)."""
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    y_slab = _y_slab(beta, omega, substrate.eps_complex, k0, substrate.thickness)
    return np.array([float(np.imag(1.0 / (complex(y) + y_slab))) for y in y_sheets])


def mean_surface_wave(y_sheets, weights, substrate, frequency_hz):
    """`(beta, X_bar)`: the TM0 constant and surface reactance of the area-averaged sheet admittance
    `sum_s w_s Y_s` (a homogenised board), `X_bar = eta0 sqrt((beta/k0)^2 - 1)`."""
    w = np.asarray(weights, dtype=float)
    y_mean = complex(np.sum(w / w.sum() * np.asarray(y_sheets, dtype=complex)))
    mode = find_beta_from_sheet_admittance(y_mean, substrate, frequency_hz)
    return mode.beta_rad_per_m, mode.x_s_ohm


# ---------------------------------------------------------------------------------- closed form
def radiated_harmonic_geometry(beta, k_mod, gamma_rad, frequency_hz):
    """`(theta, zeta, kt)` of the n = -1 harmonic of a grating of modulation wavenumber `k_mod` seen by a
    surface wave of constant `beta` whose direction makes the angle `gamma_rad` with the modulation vector:
    `k_t = beta (cos gamma, sin gamma) - k_mod (1, 0)`. `theta` is NaN when the harmonic does not radiate."""
    k0 = free_space_wavenumber(frequency_hz)
    kx = beta * np.cos(gamma_rad) - k_mod
    ky = beta * np.sin(gamma_rad)
    kt = float(np.hypot(kx, ky))
    if kt >= k0:
        return float("nan"), float("nan"), kt
    zeta = float(np.arctan2(ky, kx) - gamma_rad)
    return float(np.arcsin(kt / k0)), zeta, kt


def leakage_alpha_first_order(beta, y_sheet_mean, c1_siemens, theta_rad, zeta_rad, substrate, frequency_hz):
    """`alpha` [Np/m] of the closed form in the module docstring. `y_sheet_mean`: the area-averaged patch
    sheet admittance [S] the surface wave of constant `beta` lives on; `c1_siemens`: the amplitude of the
    n = -1 harmonic of the sheet admittance (for `Y = Y_bar + dY cos(Kx)`: `dY / 2`; for a duty-`d`
    two-level pattern `|Y_1 - Y_0| sin(pi d) / pi`, `two_level_c1`); `theta_rad`, `zeta_rad`: the
    radiated harmonic's polar angle and in-plane angle to the propagation direction."""
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    kappa = np.sqrt(beta**2 - k0**2)
    j_per_h = c1_siemens * ETA0 * kappa / k0
    kt = k0 * np.sin(theta_rad)
    kz_air = k0 * np.cos(theta_rad)
    kz_sub = np.sqrt(substrate.eps_complex * k0**2 - kt**2)
    cot = 1.0 / np.tan(kz_sub * substrate.thickness)
    y_up_tm = omega * EPS0 / kz_air
    y_up_te = kz_air / (omega * MU0)
    y_tm = y_up_tm - 1j * (omega * EPS0 * substrate.eps_complex / kz_sub) * cot + y_sheet_mean
    y_te = y_up_te - 1j * (kz_sub / (omega * MU0)) * cot + y_sheet_mean
    bracket = np.cos(zeta_rad) ** 2 * y_up_tm / abs(y_tm) ** 2 + np.sin(zeta_rad) ** 2 * y_up_te / abs(y_te) ** 2
    p_rad_per_h2 = 0.5 * abs(j_per_h) ** 2 * bracket
    w = surface_wave_power_per_h2(beta, substrate.eps_r, substrate.thickness, frequency_hz)
    return float(p_rad_per_h2 / (2.0 * w))


def leakage_alpha_first_order_tensor(beta, y_sheet_mean, delta_siemens, theta_rad, substrate, frequency_hz):
    """`alpha` [Np/m] for a ROTATING-AXIS tensor sheet `Y = Y_bar I + delta R(a) diag(1, -1) R(a)^T` whose axis
    angle advances with the modulation phase, `a = psi / 2` (Fong / Minatti tensor impedance surface). The
    n = -1 harmonic of the tensor is `(delta/2) [[1, j], [j, -1]]`, which acting on the surface wave's field
    `E_0 rho_hat` gives a circular current of magnitude `delta E_0 / sqrt(2)`: equal TM and TE parts in
    quadrature, whatever the propagation direction (so no `zeta` dependence), each radiating with its own
    `Re(Y_up)/|Y_total|^2` -- the circularly polarised leakage. The caller supplies the mean-sheet `beta`."""
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    kappa = np.sqrt(beta**2 - k0**2)
    e0_per_h = ETA0 * kappa / k0
    kt = k0 * np.sin(theta_rad)
    kz_air = k0 * np.cos(theta_rad)
    kz_sub = np.sqrt(substrate.eps_complex * k0**2 - kt**2)
    cot = 1.0 / np.tan(kz_sub * substrate.thickness)
    y_up_tm = omega * EPS0 / kz_air
    y_up_te = kz_air / (omega * MU0)
    y_tm = y_up_tm - 1j * (omega * EPS0 * substrate.eps_complex / kz_sub) * cot + y_sheet_mean
    y_te = y_up_te - 1j * (kz_sub / (omega * MU0)) * cot + y_sheet_mean
    j_part = 0.5 * delta_siemens * e0_per_h  # |J_TM| = |J_TE| = delta E_0 / 2
    p_rad_per_h2 = 0.5 * j_part**2 * (y_up_tm / abs(y_tm) ** 2 + y_up_te / abs(y_te) ** 2)
    w = surface_wave_power_per_h2(beta, substrate.eps_r, substrate.thickness, frequency_hz)
    return float(p_rad_per_h2 / (2.0 * w))


def two_level_c1(v0, v1, duty=0.5):
    """n = +-1 Fourier amplitude of a two-level (binary) profile of duty `duty` (fraction of the period in
    level 1): `|v1 - v0| sin(pi duty) / pi`."""
    return float(abs(v1 - v0) * np.sin(np.pi * duty) / np.pi)


# ---------------------------------------------------------------------------------- Floquet balance
def two_level_fourier(y0, y1, n_max, duty=0.5):
    """`y_p`, `p = -2 n_max .. 2 n_max` (array of length `4 n_max + 1`): Fourier coefficients of a two-level
    profile (level 1 on a fraction `duty` of the period), `Y(x) = sum_p y_p exp(-j p K x)`."""
    p = np.arange(-2 * n_max, 2 * n_max + 1)
    out = np.empty(p.size, dtype=complex)
    for i, k in enumerate(p):
        if k == 0:
            out[i] = (1.0 - duty) * y0 + duty * y1
        else:
            out[i] = (y1 - y0) * (1.0 - np.exp(-2j * np.pi * k * duty)) / (2j * np.pi * k) * np.exp(1j * np.pi * k * duty)
    return out


def sinusoid_fourier(y_mean, y_amp, n_max):
    """Coefficients of `Y(x) = y_mean + y_amp cos(Kx)`."""
    p = np.arange(-2 * n_max, 2 * n_max + 1)
    out = np.zeros(p.size, dtype=complex)
    out[p == 0] = y_mean
    out[np.abs(p) == 1] = 0.5 * y_amp
    return out


def _kz(k, kt2):
    """`sqrt(k^2 - kt2)` for COMPLEX `kt2` (a leaky `beta`): the propagating branch `Re > 0` while
    `Re(k^2 - kt2) > 0` (this keeps the continuation of the bound mode on the leaky, improper sheet),
    the decaying branch `Im < 0` otherwise -- `periodic_bc.normal_wavenumber` for real `kt`."""
    arg = k * k - kt2
    root = np.sqrt(arg + 0j)
    if arg.real > 0.0:
        return root if root.real >= 0.0 else -root
    return root if root.imag <= 0.0 else -root


def _harmonic_block(beta, n, k_mod, gamma_rad, eps_c, thickness_m, omega, k0):
    """The 2 x 2 total admittance dyad (up-looking air + down-looking grounded slab) of harmonic `n`, in
    (x, y) components, for transverse wavevector `(beta cos g + n k_mod, beta sin g)`."""
    kx = beta * np.cos(gamma_rad) + n * k_mod
    ky = beta * np.sin(gamma_rad)
    kt2 = kx * kx + ky * ky
    kz_air = _kz(k0, kt2)
    kz_sub = _kz(np.sqrt(eps_c) * k0, kt2)
    cot = 1.0 / np.tan(kz_sub * thickness_m)
    y_tm = omega * EPS0 / kz_air - 1j * (omega * EPS0 * eps_c / kz_sub) * cot
    y_te = kz_air / (omega * MU0) - 1j * (kz_sub / (omega * MU0)) * cot
    ptm = np.array([[kx * kx, kx * ky], [kx * ky, ky * ky]]) / kt2
    return y_tm * ptm + y_te * (np.eye(2) - ptm)


def _balance_matrix(beta, y_coeffs, n_max, k_mod, gamma_rad, substrate, frequency_hz):
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    size = 2 * n_max + 1
    m = np.zeros((2 * size, 2 * size), dtype=complex)
    for i, n in enumerate(range(-n_max, n_max + 1)):
        for j, q in enumerate(range(-n_max, n_max + 1)):
            c = y_coeffs[(n - q) + 2 * n_max]
            m[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = c * np.eye(2)
        m[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] += _harmonic_block(
            beta, n, k_mod, gamma_rad, substrate.eps_complex, substrate.thickness, omega, k0
        )
    return m


def floquet_complex_beta(
    y_coeffs,
    period_m,
    gamma_rad,
    substrate,
    frequency_hz,
    beta_start,
    n_max=6,
    steps=8,
    tol=1e-10,
):
    """The complex `beta` (`beta' - j alpha`) of the TM0 mode of the periodic sheet `y_coeffs`
    (`two_level_fourier` / `sinusoid_fourier`, period `period_m` along x) for a wave whose real-part
    direction makes `gamma_rad` with the modulation vector. `beta_start`: the real constant of the
    UNMODULATED mean sheet (`y_coeffs[0]`); the modulation is switched on in `steps` equal steps
    (continuation, so the root that connects to the bound mode is the one found). Returns `(beta,
    n_iterations_last)`."""
    k_mod = 2.0 * np.pi / period_m
    y_coeffs = np.asarray(y_coeffs, dtype=complex)
    centre = 2 * n_max
    beta = complex(beta_start)
    its = 0
    for s in np.linspace(1.0 / steps, 1.0, steps):
        scaled = y_coeffs * s
        scaled[centre] = y_coeffs[centre]

        def det(b, scaled=scaled):
            m = _balance_matrix(b, scaled, n_max, k_mod, gamma_rad, substrate, frequency_hz)
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


def bound_beta(y_mean, substrate, frequency_hz):
    """The real TM0 constant of a uniform sheet `y_mean` (the continuation's starting point)."""
    return find_beta_from_sheet_admittance(complex(y_mean), substrate, frequency_hz).beta_rad_per_m


def _y_slab_air_check(beta, substrate, frequency_hz):
    """`Y_air(beta) + Y_slab(beta)` of the bare slab (zero at the bare TM0 root) -- used by the tests."""
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    return _y_air(beta, omega, k0) + _y_slab(beta, omega, substrate.eps_complex, k0, substrate.thickness)
