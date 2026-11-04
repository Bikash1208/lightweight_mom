"""Side quest step 1: coherent feed reflection for A-2q with a closed gate on the blocked half-ring. MODEL_ONLY.

Cylindrical-harmonic multiple-scattering model of the surface wave (beta = 1.065 k0, lossless slab wave):
* probe = point source at the origin, outgoing H0^(2); it is transparent to returning n = 0 waves except for a scattering
  coefficient s_p (|s_p|^2 = 1 - LAUNCH: the returned power the probe absorbs into its port is LAUNCH, reciprocity);
* gate = ring at radius r (front face): blocked half-plane (psi in [90, 270] deg) reflects with Gamma_c (T = 0), open half
  transmits (Gamma = 0, T = 1, ideal pass state, IL = 0);
* inside the ring: A_n H_n^(2) + B_n H_n^(1), regular at the origin => A_n = B_n (n != 0), A_0 = 1 + s_p B_0; interface:
  B(psi) = Gamma(psi) A(psi); exterior outgoing field = A(psi) T(psi) re-expanded (diffraction of the half-aperture included);
* S11 = LAUNCH * B_0 (+ Gamma_0 of the bare feed, taken 0 = ideal match; HFSS bare baseline is -15 dB);
* power into the open half beyond the gate, relative to the un-reflected 0.5; then the aperture amplitude |E_ext| (phase assumed
  compensated by the hologram) is fed to the 2q state-1 pattern (25 deg, 0 deg, RHCP) for the directivity.
Reflection phase of the gate is absorbed in the radius: Gamma_c = +|Gamma| at the plane r.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import hankel1, hankel2

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_arch_study as st  # noqa: E402

BETA, LAUNCH, N = st.BETA, st.LAUNCH, 60
PSI = np.linspace(-np.pi, np.pi, 4096, endpoint=False)
NS = np.arange(-N, N + 1)
OPEN = (np.abs(PSI) < np.pi / 2).astype(float)  # state 1: open half is psi in (-90, 90) deg
FOUR = lambda f: np.array([np.mean(f * np.exp(-1j * k * PSI)) for k in np.arange(-2 * N, 2 * N + 1)])  # noqa: E731
G_OPEN = FOUR(OPEN)
TOEP = lambda c: np.array([[c[(m - n) + 2 * N] for n in NS] for m in NS])  # noqa: E731


def solve(r, gam, sp_amp=0.0):
    """Returns S11 (complex), p_open (relative to 0.5 -> in launched-power fraction), coefficients for the exterior field."""
    sp = sp_amp
    h1, h2 = hankel1(NS, BETA * r), hankel2(NS, BETA * r)
    blocked = 1.0 - OPEN
    Gk = gam * FOUR(blocked)
    Tg = TOEP(Gk)
    # unknown B_n.  A = B + e0*(1 + (sp-1) B_0)  (A_0 = 1 + sp B_0 ; A_n = B_n)
    e0 = (NS == 0).astype(float)
    # A_n = B_n + e0 * (1 + (sp - 1) B_0)
    m = len(NS)
    Amat = np.eye(m) + np.outer(e0, e0) * (sp - 1.0)  # A = Amat B + e0
    lhs = np.diag(h1) - Tg @ np.diag(h2) @ Amat
    rhs = Tg @ (h2 * e0)
    b = np.linalg.solve(lhs, rhs)
    a = Amat @ b + e0
    s11 = LAUNCH * b[N]
    # exterior outgoing coefficients: C_m H_m(r) = sum_n T_{m-n} A_n H_n(r),  T = open sector
    Tt = TOEP(G_OPEN)
    c = np.linalg.solve(np.diag(h2), Tt @ (h2 * a))
    return s11, c, a, b


def ext_field(c, rho, psi):
    rho, psi = np.broadcast_arrays(np.asarray(rho, float), np.asarray(psi, float))
    out = np.zeros(rho.shape, complex)
    for n, cn in zip(NS, c, strict=True):
        if abs(cn) > 1e-9:
            out += cn * hankel2(n, BETA * rho) * np.exp(1j * n * psi)
    return out


def p_open(c):
    """Power in the open half (unit = power of the bare n = 0 wave), from the far-field angular pattern."""
    f = sum(cn * (1j) ** int(n) * np.exp(1j * n * PSI) for n, cn in zip(NS, c, strict=True))
    return float(np.mean(OPEN * np.abs(f) ** 2))


def pattern_gain(c, r_front, rows=3):
    poly, x, y, rc = st.aperture("hex")
    rho, psi = np.hypot(x, y), np.arctan2(y, x)
    rho_in = r_front + rows * st.PITCH
    m = (np.abs(psi) < np.pi / 2 - 1e-9) & (rho >= rho_in)
    e = np.abs(ext_field(c, rho, psi))
    re = st.redge(poly, rc, psi)
    alpha = -np.log(1 - st.ETA_D) / (2 * np.maximum(re - rho_in, 1e-3))
    amp = e * np.exp(-alpha * (rho - rho_in)) * m
    th0, ph0 = np.radians(25.0), 0.0
    j = st.jones_to_current(st.JONES["RHCP"], th0, ph0)
    cur = amp * np.exp(-1j * st.K0 * np.sin(th0) * x)
    et, ep = st.far_field(x, y, j[0] * cur, j[1] * cur)
    mt = st.metrics(et, ep, th0, ph0, "RHCP")
    return mt


if __name__ == "__main__":
    lam_g = 2 * np.pi / BETA
    r0 = 0.013
    out = {"lambda_g_mm": lam_g * 1e3, "r0_mm": r0 * 1e3, "rows": []}
    print("lambda_g", lam_g * 1e3, "mm; p_open nominal", p_open(solve(r0, 0.0)[1]))
    for gam in (1.0, 0.5, 0.3, 0.1):
        for sp_ph in (0.0,):
            for dr in np.linspace(-lam_g / 4, lam_g / 4, 17):
                r = r0 + dr
                if r < 0.0065:
                    continue
                s11, c, a, b = solve(r, gam)
                out["rows"].append({"gamma": gam, "sp_phase_deg": float(np.degrees(sp_ph)), "dr_mm": float(dr * 1e3), "r_mm": float(r * 1e3),
                                    "S11_db": float(20 * np.log10(abs(s11) + 1e-12)), "S11_phase_deg": float(np.degrees(np.angle(s11))),
                                    "p_open": p_open(c)})
    # gain for the radii of interest
    base = pattern_gain(solve(r0, 0.0)[1], r0)
    out["baseline"] = {"D_dbi": base["D_dbi"], "gain_dbi": base["D_dbi"] + 10 * np.log10(LAUNCH * p_open(solve(r0, 0.0)[1]) * st.ETA_D)}
    print("baseline", out["baseline"])
    for row in out["rows"]:
        if row["gamma"] in (1.0, 0.3) and abs(row["dr_mm"] * 1e-3 / (lam_g / 4) * 4 % 2) < 1e-6 and row["sp_phase_deg"] == 0.0 and True:
            s11, c, a, b = solve(row["r_mm"] * 1e-3, row["gamma"])
            mt = pattern_gain(c, row["r_mm"] * 1e-3)
            row["D_dbi"] = mt["D_dbi"]
            row["SLL_db"] = mt["SLL_db"]
            # gain: launched power into the open half (p_open, absolute, of the SW power) x conversion x launch
            row["gain_dbi"] = mt["D_dbi"] + 10 * np.log10(LAUNCH * row["p_open"] * st.ETA_D)
            row["gain_delta_db"] = row["gain_dbi"] - out["baseline"]["gain_dbi"]
    Path(__file__).resolve().parents[2].joinpath("results/sidequest_hybrid/gate_reflection.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8", newline="\n")
    for g in (1.0, 0.5, 0.3, 0.1):
        rr = [r for r in out["rows"] if r["gamma"] == g and r["sp_phase_deg"] == 0.0]
        print("gamma", g, [(round(r["dr_mm"], 1), round(r["S11_db"], 1), round(r["p_open"], 2)) for r in rr])
