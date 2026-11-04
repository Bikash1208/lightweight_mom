"""Side quest: switched surface-wave launcher = probe + ring of parasitic pins (each with one MA4AGP907 at its foot). MODEL_ONLY.
Model: 2-D cylindrical-wave network. Every pin couples to the TM0 surface wave with the radiation resistance R_sw(top) of the repo's pin-launch model
(`src/dispersion/pin_launch.py`: R_sw = 2 P_sw / |I|^2, R_rad = total incl. direct space radiation, X_est = wire reactance). Normalised impedance matrix
    Z~_mn = H0^(2)(beta rho_mn)  (m != n),   Z~_mm = (R_rad + j X_est)/R_sw + Z_load/R_sw ,   i_m = I_m sqrt(R_sw,m)
(sources: probe voltage; parasitic feet loaded by the diode). Far-field SW power per radian  P(psi) = |sum_m i_m exp(j beta r_m.u(psi))|^2 / (4 pi)  (mutual R = sqrt(R_m R_n) J0).
Diode: measured MA4AGP907 at 10.54 GHz: ON 2.44 - 3.86j ohm, OFF -536j ohm. beta = 1.018 k0 (bare slab around the launcher; the hologram average 1.065 is NOT used here).
Not modelled: pin-to-pin space-wave coupling, finite pin radius beyond the PinLaunch reactance estimate, the hologram's scattering back to the ring, the board edge.
"""
import functools
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")
import numpy as np
from scipy.special import hankel2, j0

from src.dispersion.pin_launch import Launcher, PinLaunch
from src.frequency import free_space_wavenumber
from src.ris import diodes_measured as dm

F0 = 10.54e9
K0 = free_space_wavenumber(F0)
Z_ON = dm.z_of_f("ON_10mA")(F0)
Z_OFF = dm.z_of_f("OFF_0V")(F0)
Z0 = 50.0


@functools.lru_cache(maxsize=None)
def pin(top_mm):
    s = PinLaunch(F0, 2.2, 1.575e-3, Launcher(top_above_surface_m=top_mm * 1e-3, pin_radius_m=0.305e-3), tan_delta=0.0009).split()
    return {"R_sw": 2 * s["p_sw_w_per_a2"], "R_rad": s["r_rad_ohm"], "X": float(s["x_est_ohm"]), "beta": s["beta_over_k0"] * K0, "launch": s["launch_fraction"]}


PROBE_CAL = True  # feed calibrated to the HFSS S11 of the single hexagon (Phase 26: R x 1.30; X set so that the bare probe has S11 = -15 dB): HFSS_CALIBRATED, one board


def probe_element(probe_top=6.0):
    e = dict(pin(probe_top))
    if PROBE_CAL:
        e["R_sw"] *= 1.30
        e["R_rad"] *= 1.30
        e["X"] = 17.8
    return e


def network(pos, tops, loads, probe_top=6.0):
    """pos: (N,2) element positions [m]; element 0 is the probe at the origin. tops/loads: per parasitic element. Returns normalised matrix and R_sw vector."""
    el = [probe_element(probe_top)] + [pin(t) for t in tops]
    beta = el[0]["beta"]
    n = len(el)
    zt = np.zeros((n, n), complex)
    rs = np.array([e["R_sw"] for e in el])
    for m in range(n):
        for k in range(n):
            if m == k:
                zl = 0.0 if m == 0 else loads[m - 1]
                zt[m, m] = (el[m]["R_rad"] + 1j * el[m]["X"] + zl) / rs[m]
            else:
                zt[m, k] = hankel2(0, beta * np.hypot(*(pos[m] - pos[k])))
    return zt, rs, beta


def solve(ring_n, rho_r, top, on_mask, probe_top=6.0, z_on=None, z_off=None, offset_deg=0.0, rot_deg=0.0, tops=None):
    z_on = Z_ON if z_on is None else z_on
    z_off = Z_OFF if z_off is None else z_off
    ang = np.radians(offset_deg + (np.arange(ring_n) + 0.5) * 360.0 / ring_n)
    pos = np.vstack([[0, 0], np.c_[rho_r * np.cos(ang), rho_r * np.sin(ang)]])
    tops = [top] * ring_n if tops is None else tops
    loads = [z_on if on else z_off for on in on_mask]
    zt, rs, beta = network(pos, tops, loads, probe_top)
    v = np.zeros(len(rs), complex)
    v[0] = 1.0
    i = np.linalg.solve(zt, v)  # normalised currents (probe driven by normalised voltage 1)
    zin = (zt[0, 0] - 0j) / 1.0
    i_p = i[0]
    z_in_norm = 1.0 / i_p
    z_in = z_in_norm * rs[0]  # ohm (probe self-impedance in the presence of the ring, incl. the probe's own R_rad + jX)
    s11 = (z_in - Z0) / (z_in + Z0)
    psi = np.linspace(-np.pi, np.pi, 1440, endpoint=False)
    u = np.c_[np.cos(psi), np.sin(psi)]
    field = np.exp(1j * beta * (u @ pos.T)) @ i  # sum i_m exp(j beta r_m . u)
    p_psi = np.abs(field) ** 2 / (4 * np.pi)
    dpsi = psi[1] - psi[0]
    p_tot = float(np.sum(p_psi) * dpsi)
    p_in = 0.5 * np.real(np.conj(i_p) * 1.0) * rs[0]  # accepted power ~ 0.5 Re(V I*) in ohm-normalised units
    # accepted power in the same normalisation: P = 0.5 Re(i^H z~ i) * R_sw summed ... use direct: P_in = 0.5 Re(v^H i) scaled by R of the probe
    p_acc = 0.5 * np.real(np.conj(v[0]) * i[0])
    # losses: diode ohmic (Re z_load) and direct space radiation of every pin (R_rad - R_sw)
    r_space = np.array([(probe_element(probe_top)["R_rad"] - rs[0])] + [(pin(t)["R_rad"] - pin(t)["R_sw"]) for t in tops]) / rs
    p_space = 0.5 * float(np.sum(np.abs(i) ** 2 * r_space))
    p_diode = 0.5 * float(sum(abs(i[k + 1]) ** 2 * (loads[k].real / rs[k + 1]) for k in range(ring_n)))
    return {"s11_db": float(20 * np.log10(abs(s11))), "z_in": complex(z_in), "psi": psi, "p_psi": p_psi, "p_sw": p_tot, "p_acc": float(p_acc), "p_space": p_space, "p_diode": p_diode,
            "i": i, "beta": beta, "pos": pos}


def sectors(res, centre_deg, half_deg):
    psi = np.degrees(res["psi"])
    d = (psi - centre_deg + 180) % 360 - 180
    return float(np.sum(res["p_psi"][np.abs(d) <= half_deg]) * (res["psi"][1] - res["psi"][0]))


def fb_2q(res):
    pf, pb = sectors(res, 0, 90), sectors(res, 180, 90)
    return float(10 * np.log10(pf / max(pb, 1e-30))), pf / (pf + pb)


if __name__ == "__main__":
    print("probe", pin(6.0))
    for t in (0.0, 2, 3, 4, 5, 6, 7, 8, 9, 10):
        p = pin(t)
        print("top", t, {k: round(v, 3) for k, v in p.items()})
    print("Z_ON", Z_ON, "Z_OFF", Z_OFF)
    # sanity: probe alone
    r = solve(4, 10e-3, 6.0, [False] * 4, z_off=1e9)
    print("probe + 4 open pins: S11 %.1f dB, SW power / accepted = %.3f" % (r["s11_db"], r["p_sw"] / r["p_acc"]), "(PinLaunch launch fraction 0.627, probe alone)")
