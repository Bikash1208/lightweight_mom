"""Method 2 (fallback): the 2-quadrant board with ONE PROBE PER QUADRANT and an ideal SPDT. 150 mm hexagon, halves x > 0 / x < 0, probe A at (+x0, 0) feeds the +x half
(beam (25, 0) RHCP), probe B at (-x0, 0) the -x half (beam (25, 180) LHCP); static ideal-cell holograms (as ARCH_STUDY.md), conversion capped at the available alpha (3.4 Np/m,
Phase 26 `alpha_available`), eta_d = 0.9 per ray inside the own half. The hologram of the other half is illuminated by the active probe with the WRONG phase front
(designed for its own probe): that is the quadrant isolation. Probe-probe coupling: SW mutual impedance sqrt(R_A R_B) H0^(2)(beta d) with the hologram attenuation exp(-alpha d).
MODEL_ONLY; ideal polarisation by construction; no edge diffraction, no probe scattering of the other probe (terminated 50 ohm: included in S11/S21 only)."""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
from scipy.special import hankel2

import sidequest_arch_study as st
import sidequest_launcher_model as L

ALPHA_CAP = 3.4
ETA_D = 0.9
BETA_H = 1.065 * st.K0


def half_edge(x_probe, side, psi, poly):
    """Distance from the probe along psi to the boundary of the own half-hexagon (hex edge or x = 0)."""
    from src.hologram.leaky_forward import polygon_ray_distance

    d_hex = polygon_ray_distance(poly, (0.0, 0.0), psi) if False else None
    # ray from the probe: use a shifted polygon
    shifted = poly - np.array([x_probe, 0.0])
    d = polygon_ray_distance(shifted, (0.0, 0.0), psi)
    c = np.cos(psi)
    d_x0 = np.where(side * c < 0, np.abs(x_probe) / np.maximum(np.abs(c), 1e-9), np.inf)
    return np.minimum(d, d_x0)


def alpha_of_psi(x_probe, side, psi, poly):
    ell = half_edge(x_probe, side, psi, poly)
    return np.minimum(-np.log(1 - ETA_D) / (2 * np.maximum(ell - st.RHO_KEEP, 1e-3)), ALPHA_CAP)


def state_currents(x, y, poly, x0, active, pols, beams):
    """Cell currents (jx, jy) when probe `active` (+1: A at +x0, -1: B at -x0) is driven. Returns also the converted fraction and own/other masks."""
    sides = {1: (x0, +1), -1: (-x0, -1)}
    xp, side = sides[active]
    xo, side_o = sides[-active]
    own = (x * side) > 0
    out_jx = np.zeros(x.shape, complex)
    out_jy = np.zeros(x.shape, complex)
    w_total = 0.0
    cells = np.flatnonzero(np.hypot(x - xp, y) >= st.RHO_KEEP)
    th_own, ph_own = beams[active]
    th_oth, ph_oth = beams[-active]
    j_own = st.jones_to_current(st.JONES[pols[active]], th_own, ph_own)
    j_oth = st.jones_to_current(st.JONES[pols[-active]], th_oth, ph_oth)
    for i in cells:
        dx, dy = x[i] - xp, y[i]
        rho = np.hypot(dx, dy)
        psi = np.arctan2(dy, dx)
        a_own = float(alpha_of_psi(xp, side, np.array([psi]), poly)[0])
        if own[i]:
            amp = np.sqrt(2 * a_own / rho) * np.exp(-a_own * (rho - st.RHO_KEEP))
            ph = -st.K0 * np.sin(th_own) * (x[i] * np.cos(ph_own) + y[i] * np.sin(ph_own))
            cur = amp * np.exp(1j * ph)
            out_jx[i], out_jy[i] = j_own[0] * cur, j_own[1] * cur
            w_total += 2 * a_own * np.exp(-2 * a_own * (rho - st.RHO_KEEP)) / rho * st.PITCH**2 / (2 * np.pi)
        else:  # other half: wrong phase front, crossing attenuation
            # crossing point on x = 0
            t = (0.0 - xp) / dx if abs(dx) > 1e-12 else 1.0
            ell_own, ell_cross = abs(t) * rho, (1 - abs(t)) * rho
            a_oth = float(alpha_of_psi(xo, side_o, np.array([np.arctan2(dy, x[i] - xo)]), poly)[0])  # that cell's own-design alpha
            amp_inc = np.exp(-a_own * ell_own - a_oth * ell_cross) / np.sqrt(rho)
            rho_o = np.hypot(x[i] - xo, y[i])
            amp = np.sqrt(2 * a_oth) * amp_inc
            ph = -st.K0 * np.sin(th_oth) * (x[i] * np.cos(ph_oth) + y[i] * np.sin(ph_oth)) - st.BETA * (rho - rho_o)
            cur = amp * np.exp(1j * ph)
            out_jx[i], out_jy[i] = j_oth[0] * cur, j_oth[1] * cur
            w_total += 2 * a_oth * amp_inc**2 * rho * st.PITCH**2 / (2 * np.pi) / rho  # radiated by the wrong-phase cell
    return out_jx, out_jy, w_total, own


if __name__ == "__main__":
    poly, x, y, _ = st.aperture("hex")
    x0 = float(np.mean(x[x > 0]))
    print("n cells", x.size, "x0 (half centroid) = %.1f mm" % (x0 * 1e3))
    beams = {1: (np.radians(25), 0.0), -1: (np.radians(25), np.pi)}
    pols = {1: "RHCP", -1: "LHCP"}
    res = {"x0_mm": x0 * 1e3, "states": {}}
    for act in (1, -1):
        jx, jy, w_total, own = state_currents(x, y, poly, x0, act, pols, beams)
        et, ep = st.far_field(x, y, jx, jy)
        th_own, ph_own = beams[act]
        th_oth, ph_oth = beams[-act]
        mt = st.metrics(et, ep, th_own, ph_own, pols[act])
        inten = np.abs(et) ** 2 + np.abs(ep) ** 2
        lk = st.leak_db(et, ep, th_oth, ph_oth, mt["_inten_peak"])
        # fraction of the aperture power radiated by the own half
        pw_own = float(np.sum(np.abs(jx[own]) ** 2 + np.abs(jy[own]) ** 2))
        pw_all = float(np.sum(np.abs(jx) ** 2 + np.abs(jy) ** 2))
        gain = mt["D_dbi"] + 10 * np.log10(L.pin(6.0)["launch"] * w_total)
        res["states"][str(act)] = {"beam_peak": [mt["peak_theta_deg"], mt["peak_phi_deg"]], "D_dbi": mt["D_dbi"], "SLL_db": mt["SLL_db"], "converted_fraction": w_total,
                                   "gain_dbi": float(gain), "other_beam_level_db": lk, "isolation_db": float(-lk), "own_current_power_share": pw_own / pw_all}
        print(act, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res["states"][str(act)].items()}, flush=True)
    # probe-probe coupling (SW mutual impedance), probes calibrated as in the launcher model, other probe terminated by 50 ohm
    d = 2 * x0
    pe = L.probe_element()
    rs, rrad, xx = pe["R_sw"], pe["R_rad"], pe["X"]
    z11 = rrad + 1j * xx
    z12 = rs * hankel2(0, BETA_H * d) * np.exp(-ALPHA_CAP * d)
    zm = np.array([[z11, z12], [z12, z11]])
    zc = 50.0
    s = (zm - zc * np.eye(2)) @ np.linalg.inv(zm + zc * np.eye(2))
    res["probe_coupling"] = {"d_mm": d * 1e3, "S11_db": float(20 * np.log10(abs(s[0, 0]))), "S21_db": float(20 * np.log10(abs(s[1, 0]))), "z12_ohm": complex(z12).__repr__(),
                             "note": "S11 here is the bare probe against a 50 ohm terminated neighbour (coupling shifts S11 by the neighbour's mutual impedance only)"}
    print(res["probe_coupling"])
    Path("results/sidequest_hybrid/twoprobe_model.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
