"""Beam pointing over the tolerance cases (RO3006 1.27, design (i), bare launch 0.71, <= 514 diodes).  MODEL_ONLY.
The hologram (cell assignment per position) is DESIGNED once for the nominal case (phase constant beta_nom, nominal cell library); every other case is the SAME board with the cell types kept (statics keep their geometry, switched cells stay switched) but the physics of that case:
its own levels, A0, beta/k0, background loss.  Reported per state: peak direction (theta, phi), peak gain, and the gain at the design directions (25, 0) / (25, 180); gains calibrated with the baseline HFSS offset (-1.78 dB).
    python sidequest_pointing.py <nominal_case>   -> results/sidequest_hybrid/k2_pointing_<nominal>.json"""
import json
import math
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
F = 10.5e9
OFFSET = -1.78


def cases():
    a = json.loads(Path("results/sidequest_hybrid/k2_tolerance_r225_summary.json").read_text())
    b = json.loads(Path("results/sidequest_hybrid/k2_tolerance_summary.json").read_text())
    out = dict(a)
    for k, nm in (("centre", "r2.20_e6.15"), ("e-", "r2.20_e6.00"), ("e+", "r2.20_e6.30")):
        out[nm] = b[k]
    return out


def setup(s):
    os.environ["K4_A0"] = "%.4f" % s["A0_cal"]
    os.environ["K4_BETA_K0"] = "%.4f" % s["beta_k0"]
    import numpy as np
    import sidequest_coded_k4 as K
    import sidequest_k2_default_design as D
    K.A0 = s["A0_cal"]
    K.BETA_K0 = s["beta_k0"]
    K.LAUNCH = 0.71
    K.AA_BG = 0.15 * (s["stack_par"][1] / 0.0009) ** 0.5 if s["stack_par"][1] > 0.0009 else 0.15
    D.setup(F, 150e-3, D.MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    return np, K, D, bd


def names_values(s):
    d = {"OFF_open": 1.0, **{k: round(float(v), 4) for k, v in s["levels"].items()}}
    return d


def nominal(a):
    warnings.simplefilter("ignore")
    s = a
    np, K, D, bd = setup(s)
    lev = sorted({round(float(v), 4) for v in s["levels"].values()} | {1.0})
    opt = D.options(lev, [(1.0, -1.0, 1.0)], [(0, 1)])
    eps = [0.0, np.radians(180)]
    lo, hi, r = 0.02, 1.5, None
    for _ in range(7):
        mid = 0.5 * (lo + hi)
        rr = D.run_design(bd, opt, eps, mid, 2.0, mu=30.0)
        if int(rr["nd"].sum()) > 514:
            lo = mid
        else:
            hi = mid
            r = rr
    if r is None:
        r = D.run_design(bd, opt, eps, hi, 2.0, mu=30.0)
    return r["m"], float(r["lam"]), int(r["nd"].sum())


def evaluate(a):
    warnings.simplefilter("ignore")
    name, s, s_nom, m_nom = a
    np, K, D, bd = setup(s)
    st = D.st
    vn, vc = names_values(s_nom), names_values(s)
    nom_levels = {round(v, 4): k for k, v in vn.items() if k != "OFF_open"}
    m = np.zeros_like(m_nom)
    sw = np.abs(m_nom[0, 0] - m_nom[1, 0]) > 1e-6
    for i in range(bd.n):
        if sw[i]:
            m[0, :, i], m[1, :, i] = 1.0, -1.0
        else:
            v = round(float(m_nom[0, 0, i]), 4)
            if abs(v - 1.0) < 1e-6:
                nv = 1.0
            else:
                key = min(nom_levels, key=lambda x: abs(x - v))
                nv = vc[nom_levels[key]]
            m[:, :, i] = nv
    loss = np.zeros((2, bd.n))
    loss[1, sw] = 1.0
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    for _ in range(3):
        amp, conv, edge = K.amp_profile(bd, m, loss)
    j, _ = bd.unit(amp)
    a_id = np.sqrt(np.mean(K.vec2(bd, m), axis=0)) * amp
    states = []
    for k in range(2):
        cur = (m[k, 0] * j["x"][0] + m[k, 1] * j["y"][0], m[k, 0] * j["x"][1] + m[k, 1] * j["y"][1])
        et, ep = st.far_field(bd.x, bd.y, *cur)
        jj = st.jones_to_current(st.JONES[bd.pols[k]], K.THETA, bd.phis[k])
        idl = a_id * np.exp(-1j * st.K0 * np.sin(K.THETA) * (bd.x * np.cos(bd.phis[k]) + bd.y * np.sin(bd.phis[k])))
        ei, epi = st.far_field(bd.x, bd.y, jj[0] * idl, jj[1] * idl)
        mi = st.metrics(ei, epi, K.THETA, bd.phis[k], bd.pols[k])
        i0 = mi["_i"]
        inten = np.abs(et) ** 2 + np.abs(ep) ** 2
        inten_i = np.abs(ei) ** 2 + np.abs(epi) ** 2
        g_design = mi["D_dbi"] + 10 * np.log10(K.LAUNCH * conv) + 10 * np.log10(inten[i0] / inten_i[i0])
        scale = 10 ** (g_design / 10) / inten[i0]
        ipk = np.unravel_index(np.argmax(inten), inten.shape)
        states.append({"theta_pk": float(np.degrees(st.TG[ipk])), "phi_pk": float(np.degrees(st.PG[ipk])), "gain_peak": float(10 * np.log10(scale * inten[ipk])) + OFFSET, "gain_design": float(g_design) + OFFSET,
                       "design_dir": [float(np.degrees(st.TG[i0])), float(np.degrees(st.PG[i0]))]})
    return {"case": name, "r": s["r"], "er": s["er"], "beta_k0": s["beta_k0"], "A0_cal": s["A0_cal"], "conv": float(conv), "edge": float(edge), "states": states}


if __name__ == "__main__":
    nom = sys.argv[1] if len(sys.argv) > 1 else "r2.25_e6.15"
    C = cases()
    with ProcessPoolExecutor(1) as ex:
        m_nom, lam, nd = ex.submit(nominal, C[nom]).result()
    print("nominal", nom, "lam %.3f diodes %d beta/k0 %.3f" % (lam, nd, C[nom]["beta_k0"]), flush=True)
    jobs = [(k, v, C[nom], m_nom) for k, v in C.items()]
    with ProcessPoolExecutor(min(len(jobs), 30)) as ex:
        res = list(ex.map(evaluate, jobs))
    res.sort(key=lambda r: (r["r"], r["er"]))
    for r in res:
        s1, s2 = r["states"]
        sinth = r["beta_k0"] - (C[nom]["beta_k0"] - math.sin(math.radians(25)))
        print("%-12s beta %.3f | s1 peak (%.0f, %.0f) %.1f dBi, design dir %.1f | s2 peak (%.0f, %.0f) %.1f, design %.1f | sin-rule theta %.1f" % (
            r["case"], r["beta_k0"], s1["theta_pk"], s1["phi_pk"], s1["gain_peak"], s1["gain_design"], s2["theta_pk"], s2["phi_pk"], s2["gain_peak"], s2["gain_design"], math.degrees(math.asin(min(sinth, 1)))), flush=True)
    Path(f"results/sidequest_hybrid/k2_pointing_{nom}.json").write_text(json.dumps({"nominal": nom, "offset_db": OFFSET, "rows": res}, indent=1), encoding="utf-8", newline="\n")
