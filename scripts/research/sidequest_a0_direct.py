"""Direct leakage rate A0 of a modulated grounded-slab sheet, first order (the -1 space harmonic radiates at theta0 = 25 deg).  MODEL_ONLY.
Sheet admittance Y(x) = Ybar + h0 * square(x) (+-h0, 50/50); mode: TM0 on the slab with Ybar (transverse resonance, repo `sheet_model`).
-1 harmonic current J = (2 h0/pi) E0;  E_-1 = -J / (Ybar + Y_up + Y_down)(kx = k0 sin th0);  P_rad = 1/2 Re(Y_up) |E_-1|^2 per area;  SW power per width S(E0) from the TM0 field;  alpha = P_rad / (2 S).
No fitted constant inside: the formula is checked against the HFSS back-out of the current cell (A0 = 0.844 Np/m at beta/k0 1.0731, dY 3.13 mS)."""
import math
import numpy as np
from scipy.optimize import brentq
from src.dispersion import sheet_model as sm
from src.constants import ETA0
EPS0 = 8.8541878128e-12
C0 = 299792458.0


def slab_y(kx, w, er, h, k0):
    kz = np.sqrt(er * k0 ** 2 - kx ** 2 + 0j)
    return -1j * (w * EPS0 * er / kz) / np.tan(kz * h)


def beta_of(ybar_s, er, tand, h, f):
    """TM0 beta (rad/m) for a sheet susceptance ybar_s [S] (imag), via the repo's transverse resonance."""
    from src.materials import Substrate
    sub = Substrate("s", er, tand, h)
    m = sm.find_beta_from_sheet_admittance(1j * ybar_s, sub, f)
    return m.beta_rad_per_m


def ybar_for_beta(beta_over_k0, er, h, f):
    w = 2 * math.pi * f; k0 = w / C0; b = beta_over_k0 * k0
    ya = 1j * w * EPS0 / math.sqrt(b * b - k0 * k0)
    return -(ya + slab_y(b, w, er, h, k0)).imag  # S (the sheet susceptance that gives this beta)


def a0_direct(ybar_s, dy_pp, er, h, f, th0_deg=25.0, tand=0.0):
    """alpha [Np/m] for a +-dy_pp/2 square modulation (complex dy_pp allowed) about Ybar = j*ybar_s; E and beta from the TM0 mode."""
    w = 2 * math.pi * f; k0 = w / C0
    b = beta_of(ybar_s, er, tand, h, f)
    aa = math.sqrt(b * b - k0 * k0)
    kz = math.sqrt(er * k0 ** 2 - b * b)
    s = b * w * EPS0 * (1 / (4 * aa ** 3) + er / (2 * kz ** 2 * math.sin(kz * h) ** 2) * (h / 2 + math.sin(2 * kz * h) / (4 * kz)))   # per |E0|^2, W/m
    kx = k0 * math.sin(math.radians(th0_deg)); kzu = math.sqrt(k0 ** 2 - kx ** 2)
    yup = w * EPS0 / kzu
    ytot = 1j * ybar_s + yup + slab_y(kx, w, er, h, k0)
    h0 = abs(dy_pp) / 2
    prad = 0.5 * yup * (2 * h0 / math.pi) ** 2 / abs(ytot) ** 2
    return prad / (2 * s), b / k0


if __name__ == "__main__":
    f = 10.5e9
    yb = ybar_for_beta(1.0731, 2.2, 1.575e-3, f)
    print("Ybar for beta/k0 1.0731 on RT5880 1.575: %.2f mS" % (yb * 1e3))
    for dy in (2.63e-3, 3.13e-3):
        print("dY %.2f mS: A0 = %.3f Np/m  (HFSS back-out 0.844)" % (dy * 1e3, a0_direct(yb, dy, 2.2, 1.575e-3, f)[0]))
