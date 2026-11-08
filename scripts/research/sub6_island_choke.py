"""SUB6 island cell: analytic bias choke per via (MODEL_ONLY / NOT_HFSS_VALIDATED).

Choke = high-Z microstrip lambda/4 (meandered) from the bottom via pad, ending in a RADIAL STUB (RF short, shunt) at its far end.  Reference plane = the ground plane
(top of the FR4, z = -3.210); the bottom copper (35 um) is the trace layer.  Z_in at the via pad:

    Z_in = Z0 (Zs + Z0 tanh(g l)) / (Z0 + Zs tanh(g l)),   g = alpha + j beta,  Zs = Z_radial(r1, r2, a) at the stub apex (continuation beyond the node open).

Formulas (all stated again in CHOKE.md):
 * microstrip (Hammerstad-Jensen, thickness-corrected width we = w + (t/pi)(1 + ln(4 pi w / t))):  eps_eff = (er+1)/2 + (er-1)/2 (1 + 12 h/we)^-1/2  [- 0.04 (1 - w/h)^2 for w/h < 1];
   Z0 = 60/sqrt(eps_eff) ln(8h/we + we/(4h))  (we/h <= 1);  Z0 = 120 pi / (sqrt(eps_eff) (we/h + 1.393 + 0.667 ln(we/h + 1.444)))  (we/h > 1);
   alpha_d = pi er (ee - 1) tan d / (lambda0 sqrt(ee) (er - 1)),  alpha_c = Rs / (Z0 we),  Rs = sqrt(pi f mu0 / sigma);  l_g/4 = lambda0 / (4 sqrt(eps_eff)).
 * radial stub (open at r2, sector angle a, height h, wavenumber k = k0 sqrt(eps_r,rad (1 - j tan d)), eta = eta0 / sqrt(eps_r,rad)):
   Zs = +j eta h / (a r1) * [N1(k r2') J0(k r1) - J1(k r2') N0(k r1)] / [N1(k r2') J1(k r1) - J1(k r2') N1(k r1)],  r2' = r2 + 0.5 h  (end fringing),  r1 = w / (2 sin(a/2));
   eps_r,rad = (er+1)/2 + (er-1)/2 (1 + 12 h / w_rad)^-1/2 with w_rad = a (r1 + r2)/2 (wide line).  The stub is sized so that Im Zs(f0) = 0 (series resonance = RF short).
 * via inductance L = 0.2 l [ln(4 l / d) + 1] nH (l, d in mm; l = 3.175 + 0.035 + t_FR4 + 0.035 long barrel, d = 0.3).
Shunt loading on the cell:  Z_v = R_via + j w L_via + Z_in, Y_load = 1 / |Z_v| (worst over 5.725-5.875 GHz).  Isolation = |Z_in| min over the band.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import jv, yv

ROOT = Path(__file__).resolve().parents[2]
ETA0, C0, MU0 = 376.730313668, 299792458.0, 4e-7 * math.pi
F0, BAND = 5.8e9, (5.725e9, 5.875e9)
ER, TAND = 4.4, 0.02
T_CU, SIGMA = 0.035e-3, 5.8e7
W_MIN, G_MIN = 0.15, 0.15
H_LIST = (0.1, 0.2, 0.4, 0.5, 0.8, 1.0, 1.6)
L_STACK_TOP = 3.175 + 0.035  # RT5880 + ground copper; + t_FR4 + 0.035 bottom copper


def ms_params(w_mm, h_mm, f=F0, er=ER, tand=TAND):
    w, h, t = w_mm * 1e-3, h_mm * 1e-3, T_CU
    u = w / h
    dw = (t / math.pi) * (1 + math.log(4 * math.pi * w / t)) if u < 1 / (2 * math.pi) else (t / math.pi) * (1 + math.log(2 * h / t))
    we = w + dw
    ue = we / h
    ee = (er + 1) / 2 + (er - 1) / 2 / math.sqrt(1 + 12 / ue)
    if ue < 1:
        ee += (er - 1) / 2 * 0.04 * (1 - ue) ** 2 * -1.0
        z0 = 60 / math.sqrt(ee) * math.log(8 / ue + ue / 4)
    else:
        z0 = ETA0 / (math.sqrt(ee) * (ue + 1.393 + 0.667 * math.log(ue + 1.444)))
    lam0 = C0 / f
    ad = math.pi * er * (ee - 1) * tand / (lam0 * math.sqrt(ee) * (er - 1))
    rs = math.sqrt(math.pi * f * MU0 / SIGMA)
    ac = rs / (z0 * we)
    return dict(eps_eff=ee, z0=z0, alpha_np_m=ad + ac, alpha_d=ad, alpha_c=ac, beta0_rad_m=2 * math.pi * f * math.sqrt(ee) / C0, quarter_wave_mm=lam0 / (4 * math.sqrt(ee)) * 1e3, w_eff_mm=we * 1e3)


def radial_z(f, r1_mm, r2_mm, ang_rad, h_mm, w_mm, er=ER, tand=TAND):
    h = h_mm * 1e-3
    r1, r2 = r1_mm * 1e-3, (r2_mm + 0.5 * h_mm) * 1e-3
    wr = ang_rad * 0.5 * (r1_mm + r2_mm)
    ee = (er + 1) / 2 + (er - 1) / 2 / math.sqrt(1 + 12 * h_mm / max(wr, 1e-6))
    k = 2 * math.pi * f / C0 * np.sqrt(ee * (1 - 1j * tand))
    eta = ETA0 / np.sqrt(ee)
    a, b = k * r1, k * r2
    num = yv(1, b) * jv(0, a) - jv(1, b) * yv(0, a)
    den = yv(1, b) * jv(1, a) - jv(1, b) * yv(1, a)
    return 1j * eta * h / (ang_rad * r1) * num / den


def stub_radius(h_mm, ang_deg, w_mm=W_MIN, f=F0):
    """smallest r2 (mm) with Im Zs(f0) = 0, crossing from capacitive (-) to inductive (+) (series resonance, RF short)."""
    ang = math.radians(ang_deg)
    r1 = w_mm / (2 * math.sin(ang / 2))
    rs = np.linspace(max(1.5 * r1, 0.4), 14.0, 600)
    im = np.array([radial_z(f, r1, r, ang, h_mm, w_mm).imag for r in rs])
    for i in range(len(rs) - 1):
        if im[i] < 0 <= im[i + 1] and abs(im[i + 1] - im[i]) < 5 * max(abs(im[i]), abs(im[i + 1]), 1e-9) + 1e3:
            return brentq(lambda r: radial_z(f, r1, r, ang, h_mm, w_mm).imag, rs[i], rs[i + 1], xtol=1e-6), r1
    return float("nan"), r1


def z_in(f, h_mm, w_mm, l_mm, r2_mm, ang_deg, tand=TAND):
    p = ms_params(w_mm, h_mm, f, tand=tand)
    ang = math.radians(ang_deg)
    r1 = w_mm / (2 * math.sin(ang / 2))
    zs = radial_z(f, r1, r2_mm, ang, h_mm, w_mm, tand=tand)
    g = p["alpha_np_m"] + 1j * p["beta0_rad_m"]
    t = np.tanh(g * l_mm * 1e-3)
    return p["z0"] * (zs + p["z0"] * t) / (p["z0"] + zs * t), zs


def via_l_nh(t_fr4_mm, d_mm=0.3):
    l = L_STACK_TOP + t_fr4_mm + 0.035
    return 0.2 * l * (math.log(4 * l / d_mm) + 1), l


def design(h_mm, ang_deg, w_mm=W_MIN):
    """quarter-wave length l (mm) that maximises the in-band minimum |Z_in|, with the stub sized for the RF short."""
    r2, r1 = stub_radius(h_mm, ang_deg, w_mm)
    if not np.isfinite(r2):
        return None
    p = ms_params(w_mm, h_mm)
    fs = np.linspace(BAND[0], BAND[1], 7)
    best = None
    for l in np.linspace(0.6 * p["quarter_wave_mm"], 1.4 * p["quarter_wave_mm"], 400):
        zi = [abs(z_in(f, h_mm, w_mm, l, r2, ang_deg)[0]) for f in fs]
        m = min(zi)
        if best is None or m > best[0]:
            best = (m, l)
    l = best[1]
    zin_band = [z_in(f, h_mm, w_mm, l, r2, ang_deg)[0] for f in fs]
    lv, lvia = via_l_nh(h_mm)
    zv = [1j * 2 * math.pi * f * lv * 1e-9 + zi for f, zi in zip(fs, zin_band)]
    return dict(h_mm=h_mm, w_mm=w_mm, angle_deg=ang_deg, r1_mm=r1, r2_mm=r2, l_mm=l, ms=p, zin_min_ohm=float(min(abs(z) for z in zin_band)),
                zin_f0=complex(z_in(F0, h_mm, w_mm, l, r2, ang_deg)[0]), zs_f0=complex(z_in(F0, h_mm, w_mm, l, r2, ang_deg)[1]), via_len_mm=lvia, via_L_nH=lv,
                Yload_max_uS=float(max(1.0 / abs(z) for z in zv) * 1e6), Zv_min_ohm=float(min(abs(z) for z in zv)))


if __name__ == "__main__":
    # formula checks
    z = radial_z(1e8, 0.1, 5.0, math.pi / 2, 0.5, 0.15)
    c = 8.854e-12 * ((4.4 + 1) / 2) * 0.5 * math.pi / 2 * ((5.25e-3) ** 2 - (1e-4) ** 2) / (0.5e-3)
    print("low-f radial Z", z, "quasi-static -j/(wC) ~", -1 / (2 * math.pi * 1e8 * c))
    for h in H_LIST:
        for a in (90, 120):
            d = design(h, a)
            if d:
                print(h, a, "Z0 %.0f ee %.2f lg/4 %.2f | r2 %.2f l %.2f |Zin|min %.0f ohm  Yload %.0f uS via %.2f mm %.2f nH" % (d["ms"]["z0"], d["ms"]["eps_eff"], d["ms"]["quarter_wave_mm"], d["r2_mm"], d["l_mm"], d["zin_min_ohm"], d["Yload_max_uS"], d["via_len_mm"], d["via_L_nH"]))
            else:
                print(h, a, "no stub root")
