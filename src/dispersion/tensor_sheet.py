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


def assemble_impedance_matrix_orders(basis, green, kx0, ky0, order_x, order_y):
    """`mom_matrix.assemble_impedance_matrix` with different Floquet orders along x and y (a long thin
    supercell needs `order_x ~ cells_x`, `order_y ~ cells_y` only; the square `floquet_order` of
    `MeshSettings` would multiply the work by the aspect ratio). Same formula, same chunking."""
    lat = basis.lattice
    mx, my = np.meshgrid(np.arange(-order_x, order_x + 1), np.arange(-order_y, order_y + 1), indexing="ij")
    kx, ky = lat.transverse_wavenumbers(kx0, ky0, mx, my)
    g = green.spectral_dyad(kx, ky)
    kx, ky = kx.ravel(), ky.ravel()
    gxx, gxy, gyy = g[..., 0, 0].ravel(), g[..., 0, 1].ravel(), g[..., 1, 1].ravel()
    nb, nxb = basis.n_unknowns, basis.n_x
    z = np.zeros((nb, nb), dtype=complex)
    chunk = max(1, 2_000_000 // nb)
    for s in range(0, kx.size, chunk):
        sl = slice(s, s + chunk)
        phi = basis.fourier_transform(kx[sl], ky[sl])
        phi_x, phi_y = phi[:nxb], phi[nxb:]
        cx, cy = np.conj(phi_x), np.conj(phi_y)
        z[:nxb, :nxb] += (cx * gxx[sl]) @ phi_x.T
        z[:nxb, nxb:] += (cx * gxy[sl]) @ phi_y.T
        z[nxb:, :nxb] += (cy * gxy[sl]) @ phi_x.T
        z[nxb:, nxb:] += (cy * gyy[sl]) @ phi_y.T
    return -z / lat.cell_area


def bloch_sheet_tensor(unit_cell, frequency_hz, mesh, kx0, ky0, basis=None, orders=None):
    """`Y_s(kx0, ky0)` (2 x 2, siemens, x/y components): the patch layer's tensor sheet admittance for a Bloch
    wave of transverse wavenumber `(kx0, ky0)` -- the evanescent regime `kt > k0` of a surface wave included.

    The cell is driven by an impressed tangential field of unit amplitude in the (0,0) Floquet harmonic
    (`E_stack`, x and y); the solved current's cell mean `J_00` and the total (0,0) tangential field
    `E_tot = E_stack + G_00 J_00` (`G_00` the grounded-slab spectral dyad, which contains the slab pole) give
    `Y_s = [J_x J_y] [E_tot,x E_tot,y]^-1`. At `kt = 0` this is the normal-incidence extraction of
    `sheet_admittance_tensor_for_cell` (tested); for `kt > k0` it is the quantity the transverse-resonance
    dispersion needs, with the cell-scale near field (angle dependence of a near-resonant cell) included."""
    from ..basis_functions import RooftopBasis
    from ..mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
    from ..solver import solve

    if getattr(unit_cell, "backing", "air") != "pec":
        raise ValueError("bloch_sheet_tensor needs a grounded cell (backing='pec')")
    green = LayeredGreenFunction(unit_cell.substrate, frequency_hz, backing="pec")
    basis = basis or RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    if orders is None:
        z = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    else:
        z = assemble_impedance_matrix_orders(basis, green, kx0, ky0, *orders)
    v = np.stack([assemble_excitation_vector(basis, e, kx0, ky0) for e in np.eye(2)], axis=1)
    currents = solve(z, v)
    j = np.stack([basis.current_coefficients(currents[:, k], kx0, ky0) for k in range(2)], axis=1)  # (2 comp, 2 rhs)
    g00 = np.asarray(green.spectral_dyad(kx0, ky0), dtype=complex)
    e_tot = np.eye(2) + g00 @ j
    return j @ np.linalg.inv(e_tot)


def bloch_beta(unit_cell, frequency_hz, mesh, phi_rad, beta_guess_over_k0=None, n_grid=7, refine=2):
    """`(beta, y_s)` of the lowest TM0-like bound mode along azimuth `phi_rad` from the Bloch sheet tensor at
    the OPERATING `kt` (module docstring): the root of `det Im(Y_s(kt) + A pp^T + B qq^T) = 0`. `Y_s(kt)` is
    sampled on a grid, interpolated, the root of the interpolant found and then refined with real MoM
    evaluations at the root (secant, `refine` steps)."""
    from scipy.interpolate import CubicSpline

    from ..basis_functions import RooftopBasis

    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    sub = unit_cell.substrate
    h = sub.thickness
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    p = np.array([np.cos(phi_rad), np.sin(phi_rad)])
    q = np.array([-np.sin(phi_rad), np.cos(phi_rad)])

    def y_of(beta):
        return bloch_sheet_tensor(unit_cell, frequency_hz, mesh, beta * p[0], beta * p[1], basis=basis)

    def residual_from(y, beta):
        a = (_y_air(beta, omega, k0) + _y_slab(beta, omega, sub.eps_complex, k0, h)).imag
        b = _b_te(beta, omega, sub.eps_complex, k0, h).imag
        m = y.imag + a * np.outer(p, p) + b * np.outer(q, q)
        return float(np.linalg.det(m))

    lo, hi = 1.004 * k0, 1.34 * k0
    grid = np.linspace(lo, hi, n_grid)
    ys = np.array([y_of(b) for b in grid])
    splines = [[CubicSpline(grid, ys[:, i, j].imag) for j in range(2)] for i in range(2)]

    def interp(beta):
        return np.array([[1j * float(splines[i][j](beta)) for j in range(2)] for i in range(2)])

    fine = np.linspace(lo, hi, 600)
    vals = np.array([residual_from(interp(b), b) for b in fine])
    idx = np.where(np.diff(np.sign(vals)) != 0)[0]
    roots = []
    for i in idx:
        r = float(brentq(lambda b: residual_from(interp(b), b), fine[i], fine[i + 1]))
        if abs(residual_from(interp(r), r)) < 1e-3 * max(abs(vals[i]), abs(vals[i + 1])):
            roots.append(r)
    if not roots:
        raise dispersion_sheet_model_no_root(f"no Bloch-dispersion root for phi={np.degrees(phi_rad):.1f} deg")

    def tm_fraction(beta):
        """Share of the mode's tangential E along the propagation direction (1 = pure TM-like)."""
        a = (_y_air(beta, omega, k0) + _y_slab(beta, omega, sub.eps_complex, k0, h)).imag
        b = _b_te(beta, omega, sub.eps_complex, k0, h).imag
        m = interp(beta).imag + a * np.outer(p, p) + b * np.outer(q, q)
        w, v = np.linalg.eigh(0.5 * (m + m.T))
        vec = v[:, int(np.argmin(np.abs(w)))]
        return float(abs(vec @ p) ** 2)

    # TM0-like mode: with a strongly capacitive sheet a TE-like bound root can also appear; take the root with
    # the largest TM content (a single hybrid root is the TM0 branch pushed by its TE neighbour, whatever its
    # fraction), the lower beta on ties
    beta = sorted(roots, key=lambda r: (-round(tm_fraction(r), 3), r))[0]
    # refine with the MoM itself (the interpolant is only a bracket)
    y = y_of(beta)
    for _ in range(refine):
        step = 2e-3 * k0
        f0 = residual_from(y, beta)
        y1 = y_of(beta + step)
        f1 = residual_from(y1, beta + step)
        if f1 == f0:
            break
        new = beta - f0 * step / (f1 - f0)
        if not (lo < new < hi):
            break
        beta, y = new, y_of(new)
    return float(beta), y


def bloch_frequencies(unit_cell, mesh, kx0, ky0, f_lo=8.0e9, f_hi=12.5e9, n_grid=7):
    """All bound-mode eigenfrequencies `[(f_hz, tm_fraction), ...]` (ascending) of the periodic cell at the FIXED
    Bloch wavevector `(kx0, ky0)` -- what an HFSS eigenmode solve with master/slave phases `kx0 P`, `ky0 P`
    returns near the surface-wave band. Roots in frequency of `det Im(Y_s(f) + A pp^T + B qq^T)`; `Y_s(f)` is
    sampled on `n_grid` frequencies, interpolated, searched on a fine grid and every root refined with real
    MoM evaluations. A TM-like and a TE-like root can lie within ~0.1 GHz of each other on a strongly loaded
    cell (the determinant then has no sign change at either alone), hence the dense search."""
    from scipy.interpolate import CubicSpline

    from ..basis_functions import RooftopBasis

    sub = unit_cell.substrate
    h = sub.thickness
    kt = float(np.hypot(kx0, ky0))
    p = np.array([kx0, ky0]) / kt
    q = np.array([-p[1], p[0]])
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    f_hi = min(f_hi, 0.998 * 299792458.0 * kt / (2 * np.pi))  # bound: kt > k0(f)

    def bg(f):
        k0 = free_space_wavenumber(f)
        om = angular_frequency(f)
        a = (_y_air(kt, om, k0) + _y_slab(kt, om, sub.eps_complex, k0, h)).imag
        b = _b_te(kt, om, sub.eps_complex, k0, h).imag
        return a, b

    def resid_y(y, f):
        a, b = bg(f)
        return float(np.linalg.det(y.imag + a * np.outer(p, p) + b * np.outer(q, q)))

    def y_of(f):
        return bloch_sheet_tensor(unit_cell, f, mesh, kx0, ky0, basis=basis)

    grid = np.linspace(f_lo, f_hi, n_grid)
    ys = np.array([y_of(f) for f in grid])
    sp = [[CubicSpline(grid, ys[:, i, j].imag) for j in range(2)] for i in range(2)]

    def interp(f):
        return np.array([[1j * float(sp[i][j](f)) for j in range(2)] for i in range(2)])

    fine = np.linspace(f_lo, f_hi, 1500)
    vals = np.array([resid_y(interp(f), f) for f in fine])
    out = []
    for i in np.where(np.diff(np.sign(vals)) != 0)[0]:
        f0 = float(brentq(lambda f: resid_y(interp(f), f), fine[i], fine[i + 1]))
        # one secant refinement with the real MoM
        d = 2e7
        r0 = resid_y(y_of(f0), f0)
        r1 = resid_y(y_of(f0 + d), f0 + d)
        if r1 != r0:
            f1 = f0 - r0 * d / (r1 - r0)
            if abs(f1 - f0) < 0.15e9:
                f0 = f1
        a, b = bg(f0)
        m = y_of(f0).imag + a * np.outer(p, p) + b * np.outer(q, q)
        w, v = np.linalg.eigh(0.5 * (m + m.T))
        out.append((f0, float(abs(v[:, int(np.argmin(np.abs(w)))] @ p) ** 2)))
    if not out:
        raise dispersion_sheet_model_no_root(f"no bound-mode eigenfrequency for k=({kx0:.1f},{ky0:.1f}) rad/m")
    return out


def bloch_frequency(unit_cell, mesh, kx0, ky0, f_lo=8.0e9, f_hi=12.5e9):
    """The TM0-like eigenfrequency: the lowest root of `bloch_frequencies` with the largest TM content."""
    roots = bloch_frequencies(unit_cell, mesh, kx0, ky0, f_lo, f_hi)
    return float(sorted(roots, key=lambda r: (-round(r[1], 2), r[0]))[0][0])


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
    "assemble_impedance_matrix_orders",
    "beta_azimuth_table",
    "bloch_beta",
    "bloch_frequencies",
    "bloch_frequency",
    "bloch_sheet_tensor",
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
