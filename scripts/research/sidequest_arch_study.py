"""Side quest item 2: idealised architecture study, hybrid HMS + sparse RIS (quadrant switching). MODEL_ONLY.

Idealised cells (no MoM): every cell is a tangential-current point source J_i = a_i m_i exp(-j phase_i) (Jx, Jy); the polarisation
Jones vector of the state is imposed ideally (perfect tensor cell). Surface-wave feed = isotropic cylindrical wave rho^-1/2
exp(-alpha rho) (the Phase 26 `leaky_forward` model, uniform conversion eta_d = 0.9 per ray in the active sector).
A: static continuous (phase-only) hologram per quadrant + gate ring (arc of `rows` cell rows at radius r_g, 1 diode / cell).
B: 1-bit amplitude hologram per quadrant, radiating cells carry a diode (non-radiating = average impedance, static).
Convention: pol Jones in the local (theta_hat, phi_hat) frame as `copol_component`: LP-H = theta_hat, LP-V = phi_hat,
RHCP/LHCP = repo E_R/E_L. Gain = D * LAUNCH * f_az * T_pass^2 * eta_conv (referred to accepted power), HFSS-OFFSET not applied.
Not modelled: TM/TE split vs ray/beam angle, edge diffraction, mutual coupling, gate cell physics (only an isolation number).
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.hologram.leaky_forward import polygon_ray_distance  # noqa: E402
from src.hologram.tensor_cp_design import F, K0, LAUNCH, PITCH, R_HEX, hexagon  # noqa: E402

BETA = 1.065 * K0  # Phase 26 mean surface-wave constant
ETA_D = 0.9
RHO_KEEP = 6.5e-3
TH = np.radians(np.arange(0.0, 90.5, 1.0))
PH = np.radians(np.arange(0.0, 360.0, 2.0))
TG, PG = np.meshgrid(TH, PH, indexing="ij")
JONES = {"LP-H": (1, 0), "LP-V": (0, 1), "RHCP": (1, 1j), "LHCP": (1, -1j)}
JONES = {k: np.array(v, dtype=complex) / np.linalg.norm(v) for k, v in JONES.items()}


def aperture(shape):
    if shape == "hex":
        poly = hexagon()
        n = int(np.ceil(R_HEX / PITCH)) + 3
        k = np.arange(-n, n + 1)
        xg, yg = np.meshgrid(PITCH * k, PITCH * k)
        x, y = xg.ravel(), yg.ravel()
        ok = np.ones(x.shape, bool)
        for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
            ok &= (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) >= -1e-9
    else:  # equal-area circle
        r = R_HEX * np.sqrt(3 * np.sqrt(3) / (2 * np.pi))
        poly = None
        n = int(np.ceil(r / PITCH)) + 2
        k = np.arange(-n, n + 1)
        xg, yg = np.meshgrid(PITCH * k, PITCH * k)
        x, y = xg.ravel(), yg.ravel()
        ok = np.hypot(x, y) <= r
    ok &= np.hypot(x, y) >= RHO_KEEP
    return poly, x[ok], y[ok], (None if poly is not None else r)


def redge(poly, rc, psi):
    return polygon_ray_distance(poly, (0.0, 0.0), psi) if poly is not None else np.full(np.shape(psi), rc)


def sector_mask(psi, centre, half):
    d = (psi - centre + np.pi) % (2 * np.pi) - np.pi
    return (d >= -half - 1e-9) & (d < half - 1e-9)


def far_field(x, y, jx, jy):
    """E_theta, E_phi on the grid for point currents (jx, jy)."""
    et = np.zeros(TG.size, complex)
    ep = np.zeros(TG.size, complex)
    th, ph = TG.ravel(), PG.ravel()
    for lo in range(0, th.size, 3000):
        s = slice(lo, lo + 3000)
        arr = np.exp(1j * K0 * np.sin(th[s])[:, None] * (np.cos(ph[s])[:, None] * x[None] + np.sin(ph[s])[:, None] * y[None]))
        sx, sy = arr @ jx, arr @ jy
        et[s] = np.cos(th[s]) * (np.cos(ph[s]) * sx + np.sin(ph[s]) * sy)
        ep[s] = -np.sin(ph[s]) * sx + np.cos(ph[s]) * sy
    return et.reshape(TG.shape), ep.reshape(TG.shape)


def jones_to_current(t, th0, ph0):
    m = np.array([[np.cos(th0) * np.cos(ph0), np.cos(th0) * np.sin(ph0)], [-np.sin(ph0), np.cos(ph0)]])
    return np.linalg.solve(m, t)


def quad_state(poly, rc, x, y, th0, ph0, pol, centre, half, rho_in, kind):
    """Unit-launched-power-per-sector current distribution of one quadrant. Returns (jx, jy, eta_conv, mask)."""
    rho, psi = np.hypot(x, y), np.arctan2(y, x)
    m = sector_mask(psi, centre, half) & (rho >= rho_in)
    re = redge(poly, rc, psi)
    alpha = -np.log(1 - ETA_D) / (2 * np.maximum(re - rho_in, 1e-3))
    amp = np.sqrt(2 * alpha / rho) * np.exp(-alpha * (rho - rho_in))
    phase = -K0 * np.sin(th0) * (x * np.cos(ph0) + y * np.sin(ph0))
    j = jones_to_current(JONES[pol], th0, ph0)
    w = amp * m
    if kind == "A":
        cur = w * np.exp(1j * phase)
    else:  # B: 1-bit amplitude hologram on the incident phase exp(-j beta rho)
        on = np.cos(BETA * rho + phase) > 0
        cur = w * on * np.exp(-1j * BETA * rho)
        w = w * on
    return j[0] * cur, j[1] * cur, m, w


def metrics(et, ep, th0, ph0, pol, p_ref=None):
    inten = np.abs(et) ** 2 + np.abs(ep) ** 2
    dth, dph = np.radians(1.0), np.radians(2.0)
    wgt = np.sin(TG) * dth * dph
    prad = float(np.sum(inten * wgt))
    i = np.unravel_index(np.argmax(inten), inten.shape)
    d_dbi = 10 * np.log10(4 * np.pi * inten[i] / prad)
    t = JONES[pol]
    co = np.abs(np.conj(t[0]) * et[i] + np.conj(t[1]) * ep[i]) ** 2
    cr = np.abs(np.conj(t[1]) * et[i] - np.conj(t[0]) * ep[i]) ** 2  # orthogonal Jones (t -> (-t1*, t0*) up to phase)
    to = np.array([-np.conj(t[1]), np.conj(t[0])])
    cr = np.abs(np.conj(to[0]) * et[i] + np.conj(to[1]) * ep[i]) ** 2
    cosd = np.cos(TG) * np.cos(TG[i]) + np.sin(TG) * np.sin(TG[i]) * np.cos(PG - PG[i])
    dist = np.degrees(np.arccos(np.clip(cosd, -1, 1)))
    hp = np.degrees(np.sqrt(4 * (4 * np.pi / 10 ** (d_dbi / 10)) / np.pi))  # HPBW estimate from D
    sll = 10 * np.log10(np.max(np.where(dist > 1.2 * hp, inten, 0)) / inten[i] + 1e-30)
    e = np.array([et[i], ep[i]])
    # axial ratio of the polarisation ellipse
    a = np.abs(e[0]) ** 2 + np.abs(e[1]) ** 2
    s3 = 2 * np.imag(np.conj(e[0]) * e[1])
    s1 = np.abs(e[0]) ** 2 - np.abs(e[1]) ** 2
    s2 = 2 * np.real(np.conj(e[0]) * e[1])
    chi = 0.5 * np.arcsin(np.clip(s3 / a, -1, 1))
    ar = float(20 * np.log10(1 / max(abs(np.tan(chi)), 1e-6))) if abs(np.tan(chi)) > 1e-6 else 99.0
    return {
        "peak_theta_deg": float(np.degrees(TG[i])),
        "peak_phi_deg": float(np.degrees(PG[i])),
        "D_dbi": float(d_dbi),
        "SLL_db": float(sll),
        "XPD_db": float(min(10 * np.log10(co / max(cr, 1e-300)), 60.0)),  # >60 dB = ideal by construction
        "AR_db": min(ar, 99.0),
        "_i": i,
        "_prad": prad,
        "_inten_peak": float(inten[i]),
    }


def leak_db(et, ep, th0, ph0, ref_peak):
    """Peak intensity within 8 deg of (th0, ph0) relative to ref_peak [dB]."""
    inten = np.abs(et) ** 2 + np.abs(ep) ** 2
    cosd = np.cos(TG) * np.cos(th0) + np.sin(TG) * np.sin(th0) * np.cos(PG - ph0)
    near = np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= 8.0
    return float(10 * np.log10(np.max(inten[near]) / ref_peak + 1e-30))


def study(shape, nq, kind, pols, rows=3, rg=0.013, iso_db=None, sparse=1.0):
    poly, x, y, rc = aperture(shape)
    beams = [(np.radians(25.0), np.radians(a)) for a in ([0, 180] if nq == 2 else [0, 90, 180, 270])]
    half = np.pi / nq
    rho = np.hypot(x, y)
    rho_in = (rg + rows * PITCH) if kind == "A" else RHO_KEEP
    q = []
    for (th0, ph0), pol in zip(beams, pols, strict=True):
        jx, jy, m, w = quad_state(poly, rc, x, y, th0, ph0, pol, ph0, half, rho_in, kind)
        if kind == "B" and sparse < 1.0:  # keep the sparse fraction of the radiating cells with the largest incident amplitude
            idx = np.flatnonzero(w > 0)
            keep = idx[np.argsort(-w[idx])[: int(round(sparse * idx.size))]]
            msk = np.zeros(w.shape, bool)
            msk[keep] = True
            jx, jy, w = jx * msk, jy * msk, w * msk
        et, ep = far_field(x, y, jx, jy)
        q.append({"jx": jx, "jy": jy, "w": w, "m": m, "et": et, "ep": ep, "beam": (th0, ph0), "pol": pol})
    f_az = 1.0 / nq
    # quadrant power normalisation: each quadrant carries LAUNCH * f_az * eta_conv of the accepted power
    res = {"shape": shape, "arch": kind, "nq": nq, "rows": rows if kind == "A" else None, "pols": pols, "states": []}
    area_q = None
    nrad = [int(np.sum(qq["w"] > 0)) for qq in q]
    for s, qq in enumerate(q):
        th0, ph0 = qq["beam"]
        mt = metrics(qq["et"], qq["ep"], th0, ph0, qq["pol"])
        # B: visible-radiation fraction of the 1-bit hologram relative to an ideal phase-only one with the same amplitudes
        kappa = 1.0
        eta = ETA_D
        if kind == "B":
            pol_j = jones_to_current(JONES[qq["pol"]], th0, ph0)
            norm = np.hypot(abs(pol_j[0]), abs(pol_j[1]))
            ph_i = -K0 * np.sin(th0) * (x * np.cos(ph0) + y * np.sin(ph0))
            # ideal reference: all cells of the sector radiate with the same amplitude model, phase-only
            rho_s, psi_s = np.hypot(x, y), np.arctan2(y, x)
            msec = qq["m"]
            re = redge(poly, rc, psi_s)
            al = -np.log(1 - ETA_D) / (2 * np.maximum(re - RHO_KEEP, 1e-3))
            aref = np.sqrt(2 * al / rho_s) * np.exp(-al * (rho_s - RHO_KEEP)) * msec
            ex, ey = far_field(x, y, pol_j[0] * aref * np.exp(1j * ph_i), pol_j[1] * aref * np.exp(1j * ph_i))
            m_ref = metrics(ex, ey, th0, ph0, qq["pol"])
            # power into the beam ~ peak intensity * Omega; compare per unit source power
            p_src_ref = float(np.sum((np.abs(pol_j[0]) ** 2 + np.abs(pol_j[1]) ** 2) * aref**2))
            p_src_b = float(np.sum(norm**2 * qq["w"] ** 2 * 1.0))
            kappa = (mt["_prad"] / p_src_b) / (m_ref["_prad"] / p_src_ref)
            # conversion stays eta_d IF the radiating cells' coupling is raised to compensate duty/thinning (cell-design
            # question, item 3); the required multiplier of the local alpha is reported as `alpha_mult`
            eta = ETA_D
        t_pass = 1.0
        gain = mt["D_dbi"] + 10 * np.log10(LAUNCH * f_az * t_pass**2 * eta * kappa)
        area = float(np.sum(qq["m"])) * PITCH**2
        d_ideal = 10 * np.log10(4 * np.pi * area / (299792458.0 / F) ** 2)
        st = {
            "state": s + 1, "beam_theta_deg": 25.0, "beam_phi_deg": float(np.degrees(ph0)), "pol": qq["pol"],
            "peak_theta_deg": mt["peak_theta_deg"], "peak_phi_deg": mt["peak_phi_deg"], "D_dbi": mt["D_dbi"],
            "gain_dbi": float(gain), "SLL_db": mt["SLL_db"], "XPD_db": mt["XPD_db"], "AR_db": mt["AR_db"],
            "eta_conv": float(eta), "kappa_vis": float(kappa), "active_cells": int(np.sum(qq["m"])), "radiating_cells": nrad[s],
            "alpha_mult": float(np.sum(qq["m"]) / max(nrad[s], 1)),
            "aperture_eff": float(10 ** ((mt["D_dbi"] - d_ideal) / 10)),
        }
        res["states"].append(st)
    # diode count
    if kind == "A":
        ring = (rho >= rg) & (rho < rg + rows * PITCH)
        res["diodes"] = int(np.sum(ring))
        res["gate_cells_per_quadrant"] = res["diodes"] / nq
    else:
        union = np.zeros(x.shape, bool)
        for qq in q:
            union |= qq["w"] > 0
        res["diodes"] = int(np.sum(union))
        res["cells_total"] = int(x.size)
        res["sparse_fraction_pct"] = 100.0 * res["diodes"] / x.size
    # leakage vs isolation / contrast (active state 1): off quadrants scaled by t (amplitude)
    sw = []
    et0, ep0 = q[0]["et"], q[0]["ep"]
    mt0 = metrics(et0, ep0, *q[0]["beam"], q[0]["pol"])
    for iso in (10, 14, 15, 16, 17, 18, 19, 20, 22, 25, 30, 40):
        t = 10 ** (-iso / 20)
        et = et0 + t * sum(qq["et"] for qq in q[1:])
        ep = ep0 + t * sum(qq["ep"] for qq in q[1:])
        lk = max(leak_db(et, ep, *qq["beam"], mt0["_inten_peak"]) for qq in q[1:])
        mm = metrics(et, ep, *q[0]["beam"], q[0]["pol"])
        sw.append({"iso_db": iso, "offbeam_peak_db": lk, "SLL_db": mm["SLL_db"], "XPD_db": mm["XPD_db"]})
    res["leak_sweep_state1"] = sw
    res["iso_for_minus15"] = next((s["iso_db"] for s in sw if s["offbeam_peak_db"] <= -15.0), None)
    for st in res["states"]:
        st.pop("_i", None)
    # S11 bound for A (reflective block): blocked azimuth fraction returns to the feed
    if kind == "A":
        res["S11_extra_amp_bound"] = float(np.sqrt(LAUNCH * (1 - f_az)))
    return res


if __name__ == "__main__":
    out = []
    sets = {2: ["RHCP", "LHCP"], 4: ["RHCP", "LHCP", "LP-H", "LP-V"]}
    for shape in ("hex", "circle"):
        for kind in ("A", "B"):
            for nq in (2, 4):
                r = study(shape, nq, kind, sets[nq])
                out.append(r)
                s1 = r["states"][0]
                print(shape, kind, nq, "diodes", r["diodes"], "iso15", r["iso_for_minus15"],
                      "g", [round(s["gain_dbi"], 1) for s in r["states"]], "sll", [round(s["SLL_db"], 1) for s in r["states"]],
                      "xpd", [round(s["XPD_db"], 1) for s in r["states"]], flush=True)
    r = [study("hex", nq, "A", sets[nq], rows=n) for nq in (2, 4) for n in (2, 4)]
    out += r
    for kind_f in (0.75, 0.5, 0.35):
        for nq in (2, 4):
            rr = study("hex", nq, "B", sets[nq], sparse=kind_f)
            rr["sparse_keep"] = kind_f
            out.append(rr)
            print("B sparse", kind_f, nq, "diodes", rr["diodes"], "g", [round(s["gain_dbi"], 1) for s in rr["states"]],
                  "sll", [round(s["SLL_db"], 1) for s in rr["states"]], flush=True)
    Path(__file__).resolve().parents[2].joinpath("results/sidequest_hybrid/arch_study.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8", newline="\n")
