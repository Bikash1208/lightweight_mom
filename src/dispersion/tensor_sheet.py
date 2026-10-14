"""Phase 26 B: the 2 x 2 TENSOR sheet admittance of an anisotropic patch layer and the surface-wave dispersion
along any in-plane direction -- the tensor generalisation of `sheet_model.py` (additive; for a scalar
sheet it reproduces `find_beta_from_sheet_admittance` bit for bit).

Extraction (normal incidence, patch-plane reference, grounded slab; same assumptions as `sheet_model`: the
patch layer is a thin shunt sheet in parallel with the slab, and its admittance does not depend on the
incidence angle -- valid for a period much smaller than a wavelength):

    Gamma  = S11 (2 x 2, E_x/E_y components, `sparameters.compute_s_parameters`)
    Y_in   = Y0 (I - Gamma)(I + Gamma)^-1                  total admittance seen from air at z = 0
    Y_s    = Y_in - Y_slab(kt = 0) I                       the patch layer's own tensor sheet admittance

Y_s is symmetric for a lossless reciprocal cell (`reciprocity_error`) and has the exact tensor rotation
property in the homogeneous limit: a patch rotated by a gives `R(a) Y_s(0) R(a)^T`
(`rotation_equivariance_error` quantifies how far a real square-lattice cell is from it).

Dispersion. In the basis (k_hat, s_hat = z x k_hat) of an in-plane propagation direction, a TM-like and a TE-like
field are coupled only by the sheet's off-diagonal element. With the background (air + slab) admittances
`A(beta)` (TM) and `B(beta)` (TE) the bound mode satisfies

    (Y_kk + A)(Y_ss + B) - Y_ks^2 = 0   <=>   Y_eff(beta) + A(beta) = 0,  Y_eff = Y_kk - Y_ks^2 / (Y_ss + B)

(`Y_eff = Y_kk` exactly when `Y_ks = 0`, and for an isotropic sheet the equation IS `sheet_model`'s). The
TM0-like root is returned together with its hybrid character (`|E_s/E_k|`: the TE content of the mode).
"""

import numpy as np
from scipy.optimize import brentq

from ..constants import MU0
from ..errors import dispersion_sheet_model_no_root
from ..excitation import ReferencePlanes
from ..frequency import angular_frequency, free_space_wavenumber
from ..green_function import LayeredGreenFunction
from ..labels import EXPERIMENTAL
from ..sparameters import compute_s_parameters, surface_wave_interval
from .sheet_model import Y0, _y_air, _y_slab

TENSOR_SHEET_LABELS = (EXPERIMENTAL, "TENSOR_SHEET_FROM_NORMAL_INCIDENCE")


# ---------------------------------------------------------------------------------- extraction
def reflection_matrix(unit_cell, frequency_hz, mesh):
    """2 x 2 patch-plane reflection matrix of a grounded cell (E_x/E_y components)."""
    if getattr(unit_cell, "backing", "air") != "pec":
        raise ValueError("the tensor sheet extraction needs a grounded cell (backing='pec')")
    planes = ReferencePlanes(
        top_air_thickness=0.0, substrate_thickness=unit_cell.substrate.thickness, bottom_air_thickness=None
    )
    sp = compute_s_parameters(unit_cell, frequency_hz, mesh, planes, warn_grating=False, warn_surface_wave=False)
    return np.asarray(sp.s11, dtype=complex)


def sheet_admittance_tensor_from_gamma(gamma, substrate, frequency_hz):
    """`Y_s` (2 x 2 complex, siemens) from a 2 x 2 patch-plane reflection matrix (module docstring)."""
    g = np.asarray(gamma, dtype=complex)
    eye = np.eye(2)
    y_in = Y0 * (eye - g) @ np.linalg.inv(eye + g)
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    return y_in - _y_slab(0.0, omega, substrate.eps_complex, k0, substrate.thickness) * eye


def sheet_admittance_tensor_for_cell(unit_cell, frequency_hz, mesh):
    return sheet_admittance_tensor_from_gamma(
        reflection_matrix(unit_cell, frequency_hz, mesh), unit_cell.substrate, frequency_hz
    )


# ---------------------------------------------------------------------------------- tensor algebra
def rotation(angle_rad):
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array([[c, -s], [s, c]])


def rotate_tensor(y, angle_rad):
    """`R(a) Y R(a)^T`: the tensor of the same cell rotated counter-clockwise by `a`."""
    r = rotation(angle_rad)
    return r @ np.asarray(y) @ r.T


def decompose(y):
    """`(y_iso, delta, axis_rad, y_offdiag_asym)` of a (real-part-reactive) symmetric tensor: the mean
    `trace/2`, the anisotropy magnitude `delta` (half the eigenvalue split, >= 0), the angle of the
    eigenvector of the LARGER eigenvalue (mod pi) and the antisymmetric part `|Y_xy - Y_yx| / 2`."""
    y = np.asarray(y)
    yi = 0.5 * (y[0, 0] + y[1, 1])
    a = 0.5 * (y[0, 0] - y[1, 1])
    b = 0.5 * (y[0, 1] + y[1, 0])
    delta = np.hypot(a.imag, b.imag) if np.iscomplexobj(y) else np.hypot(a, b)
    axis = 0.5 * np.arctan2(b.imag, a.imag) if np.iscomplexobj(y) else 0.5 * np.arctan2(b, a)
    return yi, float(delta), float(np.mod(axis, np.pi)), abs(y[0, 1] - y[1, 0]) / 2.0


def reciprocity_error(y):
    """`|Y_xy - Y_yx| / (|Y_xx| + |Y_yy|)`: 0 for a reciprocal cell."""
    y = np.asarray(y)
    return float(abs(y[0, 1] - y[1, 0]) / max(abs(y[0, 0]) + abs(y[1, 1]), 1e-300))


def rotation_equivariance_error(y_rotated, y_reference, angle_rad):
    """`|| Y(rotated cell) - R(a) Y(reference) R(a)^T || / || Y(reference) ||` (Frobenius)."""
    pred = rotate_tensor(y_reference, angle_rad)
    return float(np.linalg.norm(np.asarray(y_rotated) - pred) / max(np.linalg.norm(y_reference), 1e-300))


# ---------------------------------------------------------------------------------- dispersion
def _b_te(beta, omega, eps_c, k0, thickness_m):
    """Background TE admittance: air side + shorted slab (the TE twin of `_y_air + _y_slab`)."""
    kz0 = -1j * np.sqrt(beta**2 - k0**2 + 0j)
    kz1 = np.sqrt(eps_c * k0**2 - beta**2)
    return kz0 / (omega * MU0) - 1j * (kz1 / (omega * MU0)) / np.tan(kz1 * thickness_m)


def tensor_mode_character(y_sheet, phi_rad, beta, substrate, frequency_hz):
    """`|E_s / E_k|` of the bound mode at `beta` (the TE content of the TM-like mode): from the second
    row of `(Y_s + diag(A, B)) E = 0`."""
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    r = rotation(phi_rad)
    y = r.T @ np.asarray(y_sheet, dtype=complex) @ r
    b = _b_te(beta, omega, substrate.eps_complex, k0, substrate.thickness)
    return float(abs(y[1, 0] / (y[1, 1] + b)))


def find_beta_tensor(y_sheet, phi_rad, substrate, frequency_hz, n_scan=2000):
    """`(beta, character)` of the TM0-like bound mode propagating at azimuth `phi_rad` on the tensor sheet
    `y_sheet` (module docstring). Raises `LwmomError` (`DISPERSION_SHEET_MODEL_NO_ROOT`) when there is no
    sign change -- never a fabricated number. For an isotropic / diagonal-in-the-propagation-frame sheet the
    scan, bracket and `brentq` call are the ones of `sheet_model.find_beta_from_sheet_admittance`."""
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    eps_c = substrate.eps_complex
    h = substrate.thickness
    green = LayeredGreenFunction(substrate, frequency_hz, backing="pec")
    lo, hi = surface_wave_interval(green)
    r = rotation(phi_rad)
    y = r.T @ np.asarray(y_sheet, dtype=complex) @ r
    y_kk, y_ks, y_ss = y[0, 0], 0.5 * (y[0, 1] + y[1, 0]), y[1, 1]
    # round-off off-diagonals (1e-19 S) of a cell that is diagonal in this frame are not coupling
    coupled = abs(y_ks) > 1e-10 * (abs(y_kk) + abs(y_ss))

    def y_eff(beta):
        if not coupled:
            return y_kk
        return y_kk - y_ks**2 / (y_ss + _b_te(beta, omega, eps_c, k0, h))

    def f_imag(beta):
        return (y_eff(beta) + _y_air(beta, omega, k0) + _y_slab(beta, omega, eps_c, k0, h)).imag

    betas = np.linspace(lo * (1.0 + 1e-4), hi * (1.0 - 1e-4), n_scan)
    vals = np.array([f_imag(b) for b in betas])
    changes = np.where(np.diff(np.sign(vals)) != 0)[0]
    roots = []
    for i in changes:
        beta = float(brentq(f_imag, betas[i], betas[i + 1]))
        # a sign change at a pole of Y_eff (TE resonance) is not a root
        if abs(f_imag(beta)) < 1e-3 * max(abs(vals[i]), abs(vals[i + 1]), 1e-300):
            roots.append(beta)
    if not roots:
        raise dispersion_sheet_model_no_root(
            f"no TM0-like root of the tensor transverse-resonance equation over ({betas[0]:.6g}, {betas[-1]:.6g}) rad/m"
        )
    beta = roots[0]
    character = tensor_mode_character(y_sheet, phi_rad, beta, substrate, frequency_hz) if coupled else 0.0
    return beta, character


def beta_azimuth_table(y_sheet, substrate, frequency_hz, n_phi=36):
    """`(phi_deg, beta/k0, character)` over `phi = 0 .. 180 deg` (the tensor is symmetric under `phi + pi`)."""
    k0 = free_space_wavenumber(frequency_hz)
    out = []
    for phi in np.linspace(0.0, np.pi, n_phi, endpoint=False):
        b, c = find_beta_tensor(y_sheet, phi, substrate, frequency_hz)
        out.append((float(np.degrees(phi)), b / k0, c))
    return out


__all__ = [
    "TENSOR_SHEET_LABELS",
    "beta_azimuth_table",
    "decompose",
    "find_beta_tensor",
    "reciprocity_error",
    "reflection_matrix",
    "rotate_tensor",
    "rotation",
    "rotation_equivariance_error",
    "sheet_admittance_tensor_for_cell",
    "sheet_admittance_tensor_from_gamma",
    "tensor_mode_character",
]
