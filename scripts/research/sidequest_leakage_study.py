"""Leakage / edge / launcher study for the K2-default board after the HFSS state-1 result (9.1 dBi, 13 % in theta < 30, 32 % at 60-90, 35 % back, strongest sidelobe -0.8 dB).  MODEL_ONLY.
Engine = design (i) of K2_DEFAULT (150 mm, LP-V, eps_B 180, g 2, mu 30), lam re-tuned so that diodes <= 514 (no more diodes).  Variants: A0 (cell leakage, Np/m for the 50/50 pattern), absorber rim (extra amplitude loss in the outer rim), launch fraction.
Calibration to the ONE HFSS result: gain offset (HFSS - model), split of the non-beam power (ring 60-90 / back / mid) and the peaking factor of the non-beam pattern (from the -0.8 dB sidelobe). All flagged.
    python sidequest_leakage_study.py [workers]  ->  results/sidequest_hybrid/k2_leakage_study.json"""
import json, os, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); warnings.simplefilter("ignore")
os.environ.setdefault("K4_A0", "0.844"); os.environ.setdefault("K4_BETA_K0", "1.0731")
import numpy as np
import sidequest_coded_k4 as K
import sidequest_k2_default_design as D
HF = {"gain": 9.1, "beam30": 0.13, "ring": 0.32, "back": 0.35, "sll": -0.8, "lam": 0.28}
R, F = 150e-3, 10.5e9
RIM_W = 12.8e-3


def build():
    stat, _p, _n = D.library()
    levels = sorted({round(v[0], 4) for v in stat.values()})
    D.setup(F, R, D.MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    opt = D.options(levels, [(1.0, -1.0, 1.0)], [(0, 1)])
    return bd, opt


def job(v):
    warnings.simplefilter("ignore")
    K.A0 = v.get("A0", 0.844)
    K.LAUNCH = v.get("launch", 0.627)
    bd, opt = build()
    rim = np.zeros(bd.n)
    if v.get("rim_db"):  # absorber rim: power absorption of rim_db [dB] across the outer RIM_W (amplitude loss a: 2 a w = ln(10^(db/10)))
        a = np.log(10 ** (v["rim_db"] / 10)) / (2 * RIM_W)
        rim = np.where(bd.rho > bd.redge - RIM_W, a, 0.0)
    orig = K.amp_profile
    def ap(bd_, m, loss):
        alpha = K.A0 * np.mean(K.vec2(bd_, m), axis=0)
        aa = K.AA_BG + 2 * K.AA0 * np.mean(loss, axis=0) + rim
        integ, ie = bd_.profile(alpha, aa)
        amp = np.sqrt(2 * K.A0 / bd_.rho) * np.exp(-integ)
        w = 2 * alpha * np.exp(-2 * integ) / bd_.rho * D.st.PITCH ** 2 / (2 * np.pi)
        return amp, float(w.sum()), float(np.mean(np.exp(-2 * ie)))
    K.amp_profile = ap
    eps = [0.0, np.radians(180)]
    lo, hi = 0.02, 1.2
    def run(l): return D.run_design(bd, opt, eps, l, 2.0, mu=30.0)
    r = run(v.get("lam", 0.28)) if v.get("lam") else None
    if r is None:
        for _ in range(7):
            mid = 0.5 * (lo + hi); rr = run(mid)
            if int(rr["nd"].sum()) > 514: lo = mid
            else: hi = mid; r = rr
        if r is None: r = run(hi)
    nd = int(r["nd"].sum())
    iso = D.isolation(r["states"])
    return {**v, "lam_used": float(r["lam"]), "diodes": 2 * 0 + nd, "gain_model": [s["gain_dbi"] for s in r["states"]], "sll_model": [s["sll_db"] for s in r["states"]], "conv": r["conv"], "edge_surv": r["edge"], "iso": iso}


if __name__ == "__main__":
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    V = [{"tag": "baseline", "lam": 0.28},
         *[{"tag": f"A0_{a}", "A0": a} for a in (1.5, 2.5, 4.0, 6.0, 9.0)],
         *[{"tag": f"rim_{d}dB", "lam": 0.28, "rim_db": d} for d in (10, 20)],
         *[{"tag": f"launch_{l}", "lam": 0.28, "launch": l} for l in (0.8, 0.9)],
         {"tag": "A0_4.0+rim10+launch0.8", "A0": 4.0, "rim_db": 10, "launch": 0.8}]
    with ProcessPoolExecutor(min(nw, len(V))) as ex:
        res = list(ex.map(job, V))
    Path("results/sidequest_hybrid/k2_leakage_study.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        print(r["tag"], "lam %.3f diodes %d gain %s conv %.3f edge_surv %.3f" % (r["lam_used"], r["diodes"], [round(x, 2) for x in r["gain_model"]], r["conv"], r["edge_surv"]), flush=True)
