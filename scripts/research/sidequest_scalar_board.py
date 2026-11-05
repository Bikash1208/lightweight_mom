"""Fallback: the current rim-disc cell as a SCALAR 1-bit switch (both diodes on one control line): HFSS (2026-11-05, 90 deg, x run) C/C 10.670, R/R 11.006 GHz -> isotropic sheet switch with Y 10.07 -> 6.94 mS (dY -3.13 mS,
alpha_beam 0.844 Np/m, beta/k0 1.0731 from `hfss_alpha_backout.py`). One switched cell = 2 diodes on one line; element value m_x = m_y = m (J = m rho_hat E: linear polarisation along the radial direction), K = 2: same polarisation (LP-H, theta_hat)
in two directions (phi 0 / 180, theta 25). Static cells: analog levels (round-11 library range, normal-incidence shapes: flagged). Engine = sidequest_coded_k4.  MODEL_ONLY; the loss factor AA0 is the MoM value (HFSS Q pending).
Output results/sidequest_hybrid/scalar_board.json"""
import itertools
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("K4_A0", "0.844")
os.environ.setdefault("K4_BETA_K0", "1.0731")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_coded_k4 as K

STAT, PAIRS = K.library()
VS_AN = np.linspace(STAT.min(), STAT.max(), 9)


def options(Kn, vs, with_switch=True):
    mx, my, nd, ld, wx, wy = [], [], [], [], [], []
    z = np.zeros(Kn)
    for a in vs:
        mx.append(np.full(Kn, a)); my.append(np.full(Kn, a)); nd.append(0); ld.append((z, z)); wx.append(-1); wy.append(-1)
    W = K.words(Kn)
    if with_switch:
        for wi, w in enumerate(W):
            val = np.array([-1.0 if b else 1.0 for b in w])  # ON = -1, OFF = +1 (dY < 0 when ON)
            on = np.array(w, float)
            mx.append(val); my.append(val); nd.append(2); ld.append((on, on)); wx.append(wi); wy.append(wi)
    return np.array(mx), np.array(my), np.array(nd, float), np.array(ld), np.array(wx), np.array(wy), W


def run(pol, lam, vs, eps, g=2.0, switch=True):
    bd = K.KBoard([0, 180], [pol, pol])
    r, _ = K.run(bd, options(2, vs, switch), eps, lam=lam, g=g, iters=3)
    s = K.summary(r)
    s["lam"] = lam
    return s


if __name__ == "__main__":
    out = {"cell": {"A0": K.A0, "beta_k0": K.BETA_K0, "AA0": K.AA0, "source": "HFSS 2026-11-05 C/C 10.670 / R/R 11.006 GHz (90 deg, x run)"}, "rows": []}
    for pol in ("LP-H", "LP-V"):
        for eps in ([0.0, 0.0], [0.0, np.pi]):
            for vs_name, vs in (("bit", np.array([-1.0, 1.0])), ("analog", VS_AN)):
                for lam in (0.1, 0.2, 0.28, 0.35, 0.55, 0.8, 1.2, 2.0, 1e3):
                    s = run(pol, lam, vs, eps, switch=lam < 1e2)
                    s.update({"pol": pol, "eps_B_deg": round(float(np.degrees(eps[1]))), "statics": vs_name})
                    out["rows"].append(s)
                    print(pol, s["eps_B_deg"], vs_name, lam, {k: s[k] for k in ("gain", "SLL", "XPD", "diodes", "static_pct", "lines_x_y", "ap_eff_pct")}, flush=True)
    Path("results/sidequest_hybrid/scalar_board.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
