"""Axisymmetric finite-difference launcher: z-pin (radius a) through a grounded slab, uniform sheet admittance j*ys on rho >= rho_in (loaded annulus, bare inside), absorber layers at the top and the outer edge.  MODEL_ONLY.
Unknown u = rho * H_phi:  d_z(eps^-1 d_z u) + rho d_rho(eps^-1 rho^-1 d_rho u) + k0^2 u = 0;  pin cells (rho < a, z <= t): u = I(z)/(2 pi) prescribed; axis: zero flux; ground: d_z u = 0; sheet: u jumps by (ys/(w eps0)) * flux.
Outputs: P_in (flux out of the pin), Z_in = -(1/I0^2) int E_z I dz (complex, thin-pin model), launch = P_SW / P_in with P_SW from the TM0 projection of H_phi at rho_m.
Validation: bare and uniform loading against the spectral model `sidequest_launcher`."""
import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
EPS0, C0 = 8.8541878128e-12, 299792458.0


def solve(f, er, tand, h, t, ys, rho_in, a=0.635e-3, rmax=0.10, ztop=0.030, drho=0.127e-3, nsl=10, ab_r=0.025, ab_z=0.012, sig=2.0, fit_zmax=0.010, rho_m=None, disc_r=0.0):
    w = 2 * math.pi * f; k0 = w / C0
    dz = h / nsl; naz = int(round(ztop / dz)); nz = nsl + naz
    zc = np.concatenate([-h + (np.arange(nsl) + .5) * dz, (np.arange(naz) + .5) * dz])
    nr = int(round(rmax / drho)); rc = (np.arange(nr) + .5) * drho
    npin = int(round(a / drho)); a_eff = npin * drho
    R, Z = np.meshgrid(rc, zc, indexing="ij")
    base = np.where(Z < 0, er * (1 - 1j * tand), 1.0 + 0j)
    ramp_r = np.clip((R - (rmax - ab_r)) / ab_r, 0, 1); ramp_z = np.clip((Z - (ztop - ab_z)) / ab_z, 0, 1)
    eps = base * (1 - 1j * sig * np.maximum(ramp_r, ramp_z) ** 3)
    idx = np.arange(nr * nz).reshape(nr, nz)
    rows, cols, vals = [], [], []
    rhs = np.zeros((nr, nz), complex)
    diag = np.full((nr, nz), k0 * k0, complex)
    ef = 2.0 / (1 / eps[:-1] + 1 / eps[1:]); rf = rc[:-1] + 0.5 * drho
    cE = np.zeros((nr, nz), complex); cW = np.zeros((nr, nz), complex)
    cE[:-1] = (rc[:-1, None] / (drho * drho)) / (ef * rf[:, None]); cW[1:] = (rc[1:, None] / (drho * drho)) / (ef * rf[:, None])
    diag -= cE + cW
    rows += [idx[:-1].ravel()]; cols += [idx[1:].ravel()]; vals += [cE[:-1].ravel()]
    rows += [idx[1:].ravel()]; cols += [idx[:-1].ravel()]; vals += [cW[1:].ravel()]
    eN = eps[-1]; diag[-1] -= (rc[-1] / (drho * 0.5 * drho)) / (eN * rmax)     # outer wall u = 0
    rcz = 0.5 * dz / eps[:, :-1] + 0.5 * dz / eps[:, 1:]
    rcz[:, nsl - 1] = rcz[:, nsl - 1] + (ys / (w * EPS0)) * (rc >= rho_in)
    gz = 1.0 / rcz
    gz[rc < disc_r, nsl - 1] = 0.0                                              # printed PEC top disc (radius disc_r) at z = 0: no flux through it
    up = np.zeros((nr, nz), complex); dn = np.zeros((nr, nz), complex)
    up[:, :-1] = gz / dz; dn[:, 1:] = gz / dz
    diag -= up + dn
    rows += [idx[:, :-1].ravel()]; cols += [idx[:, 1:].ravel()]; vals += [up[:, :-1].ravel()]
    rows += [idx[:, 1:].ravel()]; cols += [idx[:, :-1].ravel()]; vals += [dn[:, 1:].ravel()]
    diag[:, -1] -= 1.0 / (dz * 0.5 * dz / eps[:, -1])                           # top wall u = 0
    rows += [idx.ravel()]; cols += [idx.ravel()]; vals += [diag.ravel()]
    r_ = np.concatenate(rows); c_ = np.concatenate(cols); v_ = np.concatenate(vals)
    if t > 0:
        I = np.array([(math.sin(k0 * (t - z)) if z > 0 else math.sin(k0 * t)) if z <= t else 0.0 for z in zc]); i0 = math.sin(k0 * t)
    else:  # pin ends at the top surface: uniform current through the slab (feeds the printed disc)
        I = np.array([1.0 if z < 0 else 0.0 for z in zc]); i0 = 1.0
    pin = np.zeros((nr, nz), bool); pin[:npin, :] = (zc <= t)[None, :]
    pin_i = idx[pin]
    keep = ~np.isin(r_, pin_i)                                                  # rows of the pin cells are replaced by identity
    r_ = np.concatenate([r_[keep], pin_i]); c_ = np.concatenate([c_[keep], pin_i]); v_ = np.concatenate([v_[keep], np.ones(len(pin_i), complex)])
    rhs[pin] = np.broadcast_to(I / (2 * math.pi), (nr, nz))[pin]
    A = sp.csc_matrix((v_, (r_, c_)), shape=(nr * nz,) * 2)
    u = spsolve(A, rhs.ravel()).reshape(nr, nz)
    du = (u[npin] - u[npin - 1]) / drho                                          # at rho = a_eff
    zi = zc <= t
    ez = -1j * du / (w * EPS0 * a_eff * eps[npin])
    zin = -np.sum((ez * I)[zi] * dz) / i0 ** 2
    pin_w = -0.5 * np.real(np.sum((ez * I)[zi] * dz))
    # energy-based input power: everything is finally absorbed (absorbers + slab tan d): P = 1/2 w eps0 sum Im(-eps) (|E_rho|^2 + |E_z|^2) dV
    Hh = u / rc[:, None]
    dHdz = np.gradient(Hh, dz, axis=1); dudr = np.gradient(u, drho, axis=0)
    Er = 1j * dHdz / (w * EPS0 * eps); Ez = -1j * dudr / (rc[:, None] * w * EPS0 * eps)
    vol = 2 * math.pi * rc[:, None] * drho * dz
    absmask = np.maximum(ramp_r, ramp_z) > 0
    p_abs = float(np.sum((0.5 * w * EPS0 * (-eps.imag) * (abs(Er) ** 2 + abs(Ez) ** 2) * vol)[absmask]))   # absorber cells only (the slab tan d is negligible; cells at the sheet jump are excluded)
    rho_m = rho_m or (rmax - ab_r - 0.012)
    j = int(np.argmin(abs(rc - rho_m))); rm = rc[j]
    from sidequest_a0_direct import beta_of
    beta = beta_of(ys if ys > 0 else 1e-9, er, tand, h, f)
    kz = math.sqrt(er * k0 ** 2 - beta ** 2); al = math.sqrt(beta ** 2 - k0 ** 2)
    ha = kz * math.sin(kz * h) / (er * al)
    psi = np.where(zc < 0, np.cos(kz * (zc + h)), ha * np.exp(-al * np.maximum(zc, 0)))
    sel = zc <= fit_zmax
    Hphi = u[j] / rm
    A_ = np.vdot(psi[sel], Hphi[sel]) / np.vdot(psi[sel], psi[sel])
    epsr = np.where(zc < 0, er, 1.0)
    psw = 2 * math.pi * rm * 0.5 * beta / (w * EPS0) * abs(A_) ** 2 * np.sum(abs(psi) ** 2 / epsr * dz)
    return {"P_in": pin_w, "R_in": 2 * pin_w / i0 ** 2, "Z_in": zin, "P_sw": psw, "launch": psw / pin_w, "P_abs": p_abs, "launch_energy": psw / p_abs, "beta_k0": beta / k0, "fit_resid": float(np.linalg.norm(Hphi[sel] - A_ * psi[sel]) / np.linalg.norm(Hphi[sel]))}


if __name__ == "__main__":
    f = 10.5e9
    for ys in (0.0, 7.38e-3):
        r = solve(f, 2.2, 0.0009, 1.575e-3, 6e-3, ys, 0.0)
        print("RT5880 1.575 ys %.2f mS (spectral: bare launch 0.611 R_in 42.4 / loaded 0.846, 56.0):" % (ys * 1e3), {k: (complex(round(v.real, 2), round(v.imag, 2)) if isinstance(v, complex) else round(float(v), 3)) for k, v in r.items()}, flush=True)
