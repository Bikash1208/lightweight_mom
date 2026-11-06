"""SWHB: Born-approximation board model of a single-feed centre-probe surface-wave (TM) holographic board on a grounded slab.  MODEL_ONLY.
Physics (no fitted constants):
 * mean sheet admittance per state -> TM0 beta, alpha_loss (transverse resonance, lossy slab + Re(Y)),
 * radial SW field E_i = A(rho) exp(-j phi(rho)),  A from the SW power P(rho) = 2 pi rho s |E|^2 (s: TM0 power per |E|^2 and width),
 * each cell radiates a sheet current J_i = dY_i E_i rho_hat_i (dY_i = Y_i - mean) with the exact grounded-slab element factor
   (E_theta ~ J_tm / Y_tm,  E_phi ~ cos(theta) J_te / Y_te, Y = Y_up + Y_dn + Ybar_sheet),
 * far field = coherent sum over cells; radiated power from the hemisphere integral; depletion of the SW per ring from the INCOHERENT element power (diffuse leakage of a quantised pattern)
   scaled by the global coherent/incoherent ratio kappa_c (iterated), plus the absorbed power Re(Y)|E|^2 and the dielectric loss,
 * launch efficiency eta_L is a parameter (the probe is not modelled).
Units SI; Y in siemens per square."""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import hankel2

C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 4e-7 * math.pi
ETA0 = math.sqrt(MU0 / EPS0)


class Env:
    def __init__(self, er, tand, h, f=10.5e9):
        self.er, self.tand, self.h, self.f = er, tand, h, f
        self.w = 2 * math.pi * f
        self.k0 = self.w / C0
        self.ec = er * (1 - 1j * tand)

    def kz1(self, kr):
        return self.k0 * np.sqrt(self.ec - kr * kr + 0j)

    def yup_tm(self, kr):
        return self.w * EPS0 / (self.k0 * np.sqrt(1 - kr * kr + 0j))

    def yup_te(self, kr):
        return self.k0 * np.sqrt(1 - kr * kr + 0j) / (self.w * MU0)

    def ydn_tm(self, kr):
        kz = self.kz1(kr)
        return -1j * (self.w * EPS0 * self.ec / kz) / np.tan(kz * self.h)

    def ydn_te(self, kr):
        kz = self.kz1(kr)
        return -1j * (kz / (self.w * MU0)) / np.tan(kz * self.h)

    def tm_power(self, kr):
        """TM0 power per unit width per |E0|^2 (peak, tangential E at the sheet), lossless slab."""
        b = kr * self.k0
        aa = math.sqrt(b * b - self.k0 ** 2)
        kz = math.sqrt(self.er * self.k0 ** 2 - b * b)
        return b * self.w * EPS0 * (1 / (4 * aa ** 3) + self.er / (2 * kz ** 2 * math.sin(kz * self.h) ** 2) * (self.h / 2 + math.sin(2 * kz * self.h) / (4 * kz)))

    def beta_alpha(self, ysheet_of_kr):
        """beta/k0 and alpha_loss [Np/m] for a sheet admittance function Y(kr) (complex S)."""
        def F(kr):
            return self.yup_tm(kr) * 0 + 1j * self.w * EPS0 / (self.k0 * math.sqrt(kr * kr - 1.0)) + self.ydn_tm(kr) + ysheet_of_kr(kr)

        ks = np.linspace(1.0 + 1e-3, 1.8, 800)
        v = np.array([F(k).imag for k in ks])
        ix = np.where(np.diff(np.sign(v)) != 0)[0]
        if len(ix) == 0:
            return None
        kr = brentq(lambda k: F(k).imag, ks[ix[0]], ks[ix[0] + 1], xtol=1e-10)
        d = (F(kr * (1 + 1e-6)).imag - F(kr * (1 - 1e-6)).imag) / (2e-6 * kr * self.k0)
        return kr, abs(F(kr).real / d)


def element_q(env, ysh_tm, ysh_te, nth=181, nph=72):
    """Radiated power of ONE cell per unit |J~|^2 [W per (A m)^2], and the angular grid weights (incoherent element power)."""
    th = (np.arange(nth) + 0.5) * (math.pi / 2) / nth
    ph = (np.arange(nph) + 0.5) * 2 * math.pi / nph
    kr = np.sin(th)
    ytm = env.yup_tm(kr) + env.ydn_tm(kr) + ysh_tm
    yte = env.yup_te(kr) + env.ydn_te(kr) + ysh_te
    dth, dph = (math.pi / 2) / nth, 2 * math.pi / nph
    # J along x: TM component cos(phi), TE component -sin(phi)
    ang_tm = np.sum(np.cos(ph) ** 2) * dph
    ang_te = np.sum(np.sin(ph) ** 2) * dph
    integ = np.sum((ang_tm / np.abs(ytm) ** 2 + np.cos(th) ** 2 * ang_te / np.abs(yte) ** 2) * np.sin(th)) * dth
    return env.k0 ** 2 / (8 * math.pi ** 2 * ETA0) * integ


def farfield(env, xy, J, ysh_tm, ysh_te, th, ph, p_cell, chunk=2000):
    """E_theta, E_phi (times r e^{jk0 r}) for cell currents J (n,2) [A/m] at xy (n,2) on the angle grid th (nt,), ph (np,) [rad]. Returns (nt, np) arrays."""
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    kx = env.k0 * np.sin(TH) * np.cos(PH)
    ky = env.k0 * np.sin(TH) * np.sin(PH)
    kr = np.sin(TH)
    ytm = env.yup_tm(kr) + env.ydn_tm(kr) + ysh_tm
    yte = env.yup_te(kr) + env.ydn_te(kr) + ysh_te
    sx = np.sinc(kx * p_cell / (2 * math.pi)) * np.sinc(ky * p_cell / (2 * math.pi))
    Jx = np.zeros(TH.shape, complex)
    Jy = np.zeros(TH.shape, complex)
    kxf, kyf = kx.ravel(), ky.ravel()
    for a in range(0, len(xy), chunk):
        ph_ = np.exp(1j * (np.outer(xy[a:a + chunk, 0], kxf) + np.outer(xy[a:a + chunk, 1], kyf)))
        Jx += (J[a:a + chunk, 0] @ ph_).reshape(TH.shape)
        Jy += (J[a:a + chunk, 1] @ ph_).reshape(TH.shape)
    Jx, Jy = Jx * p_cell ** 2 * sx, Jy * p_cell ** 2 * sx
    jtm = Jx * np.cos(PH) + Jy * np.sin(PH)
    jte = -Jx * np.sin(PH) + Jy * np.cos(PH)
    eth = 1j * env.k0 / (2 * math.pi) * (-jtm / ytm)
    eph = 1j * env.k0 * np.cos(TH) / (2 * math.pi) * (-jte / yte)
    return eth, eph


def radiated_power(eth, eph, th, ph):
    TH, _ = np.meshgrid(th, ph, indexing="ij")
    u = (np.abs(eth) ** 2 + np.abs(eph) ** 2) / (2 * ETA0)
    dth = th[1] - th[0]
    dph = ph[1] - ph[0]
    return float(np.sum(u * np.sin(TH)) * dth * dph), u


def run_state(env, xy, ybr, ycom_sw, ycom_rad, p, eta_l=1.0, rho_keep=6.5e-3, n_iter=5, grid=(2.0, 5.0), want_field=True, taper=None):
    """One state of a board.  xy (n,2) [m]; ybr (n,) complex: the per-cell admittance that varies between cells (or states); ycom_sw(kr): common sheet admittance for the SW solve (complex S, may be 0);
    ycom_rad = (tm, te) common admittance for the radiation directions.  Returns dict with fields, powers (fractions of the LAUNCHED power), beta, profile."""
    n = len(xy)
    rho = np.hypot(xy[:, 0], xy[:, 1])
    ybar = ybr.mean()
    dY = ybr - ybar
    if taper is not None:
        dY = dY * taper
    ba = env.beta_alpha(lambda kr: ycom_sw(kr) + ybar)
    if ba is None:
        return {"ok": False}
    kr_sw, a_loss = ba
    beta = kr_sw * env.k0
    s = env.tm_power(kr_sw)
    ysh_tm, ysh_te = ycom_rad[0] + ybar, ycom_rad[1] + ybar
    Q = element_q(env, ysh_tm, ysh_te)
    q = Q * np.abs(dY) ** 2 * p ** 4  # radiated power per cell per |E|^2 (incoherent)
    order = np.argsort(rho)
    rho_s = rho[order]
    ring = np.floor((rho_s - rho_s[0]) / p).astype(int)
    nr = ring.max() + 1
    kappa = 1.0
    th = np.arange(1.0, 90.0, grid[0]) * math.pi / 180
    phg = np.arange(0.0, 360.0, grid[1]) * math.pi / 180
    for it in range(n_iter):
        P = np.zeros(nr + 1)
        P[0] = 1.0
        Ecell = np.zeros(n)
        for k in range(nr):
            idx = order[ring == k]
            rk = rho[idx].mean()
            amp2 = P[k] / (2 * math.pi * rk * s)  # |E|^2 in the ring (peak)
            Ecell[idx] = amp2
            drop = kappa * q[idx].sum() * amp2 / 1.0 + 2 * a_loss * p * P[k]
            P[k + 1] = max(P[k] - drop, 0.0)
        if not want_field:
            break
        # phase and field vector
        arg = hankel2(1, beta * np.maximum(rho, 1e-4))
        ph_c = arg / np.abs(arg)
        E = np.sqrt(Ecell) * ph_c
        rh = xy / np.maximum(rho, 1e-9)[:, None]
        J = (dY * E)[:, None] * rh
        eth, eph = farfield(env, xy, J, ysh_tm, ysh_te, th, phg, p)
        prad, u = radiated_power(eth, eph, th, phg)
        sinc_sum = float(np.sum(q * Ecell))
        newk = prad / max(sinc_sum, 1e-30)
        if abs(newk - kappa) < 1e-3 * kappa:
            kappa = newk
            break
        kappa = 0.5 * (kappa + newk)
    out = {"ok": True, "beta_k0": kr_sw, "alpha_loss": a_loss, "kappa": kappa, "P_edge": float(P[nr]), "rho_ring": [], "P_ring": P.tolist(), "s": s}
    if want_field:
        out.update({"eth": eth, "eph": eph, "th": th, "ph": phg, "u": u, "P_rad_unit": prad})
    out["Pinc_sum"] = float(np.sum(q * Ecell)) if want_field else None
    out["rho_cells"] = rho
    out["Ecell2"] = Ecell
    return out
