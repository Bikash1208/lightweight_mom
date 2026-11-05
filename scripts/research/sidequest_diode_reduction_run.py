"""Run the fewer-diodes study: baselines (lambda sweeps), D1, D2, D3, D5, D6 alone and combinations, K = 2 and K = 4 (CP set B). Output results/sidequest_hybrid/diode_reduction.json (list of rows)."""
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_k4 as K
import sidequest_diode_reduction as D

OUT = Path("results/sidequest_hybrid/diode_reduction.json")
_OPT = {}


def run_spec(spec):
    try:
        r = D.execute(spec)
    except Exception as e:  # noqa: BLE001
        r = {"spec": spec, "error": repr(e)}
    return r


def _cost(a):
    name, lam, dphi, dth, eps = a
    if name not in _OPT:
        cls = D.CLS3 if name == "K4b" else None
        _OPT[name] = K.options(len(K.SETS[name][0]), D.VS_AN, D.G2, cls_x=cls, cls_y=cls)
    bd = D.board(name, dphi, dth)
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    return D.design_g(bd, amp, _OPT[name], eps, lam, 0.05, 2.0)[4]


def d5_parallel(pool, name, lam, passes=2):
    Kn = len(K.SETS[name][0])
    dphi, dth, eps = np.zeros(Kn), np.zeros(Kn), np.array(D.EPS_B[name], float)
    cur = _cost((name, lam, tuple(dphi), tuple(dth), eps))
    for _ in range(passes):
        for k in range(Kn):
            cands = []
            for a, b, e in itertools.product((-2.0, 0.0, 2.0), (-2.0, 0.0, 2.0), np.radians(np.arange(0, 360, 60)) if k else [eps[k]]):
                dp, dt, ep = dphi.copy(), dth.copy(), eps.copy()
                dp[k], dt[k], ep[k] = a, b, e
                cands.append((name, lam, tuple(dp), tuple(dt), ep))
            costs = pool.map(_cost, cands)
            i = int(np.argmin(costs))
            if costs[i] < cur - 1e-9:
                cur = costs[i]
                _, _, dp, dt, ep = cands[i]
                dphi, dth, eps = np.array(dp), np.array(dt), np.array(ep)
        print("D5", name, "pass", "cost", round(cur, 3), dphi.tolist(), dth.tolist(), np.degrees(eps).round().tolist(), flush=True)
    return dphi.tolist(), dth.tolist(), eps.tolist(), cur


def singles(name):
    lam0 = D.LAM_BASE[name]
    lams = [0.3, 0.4, 0.555, 0.8, 1.1, 1.5, 2.0] if name == "K4b" else [0.1, 0.18, 0.266, 0.4, 0.6, 0.9, 1.3]
    sets = []
    sets += [{"K": name, "id": "base", "lam": l} for l in lams]
    sets += [{"K": name, "id": "D1", "D1": f, "lam1": lam0 * 0.4} for f in (0.8, 0.6, 0.45, 0.3, 0.2, 0.1)]
    sets += [{"K": name, "id": f"D2_{m}", "D2": m, "lam": l} for m in ("1x2x", "1x2y", "2x2") for l in (lam0, lam0 * 1.5, lam0 * 2.2, lam0 * 3.3)]
    sets += [{"K": name, "id": f"D3_{m}", "D3": m, "lam": l} for m in ("one", "strict_x") for l in (lam0 * 0.6, lam0, lam0 * 1.5, lam0 * 2.2)]
    sets += [{"K": name, "id": "D6", "D6": rc, "lam": l} for rc in (35.0, 50.0, 65.0, 80.0, 95.0, 115.0) for l in (lam0 * 0.55, lam0)]
    return sets


if __name__ == "__main__":
    t0 = time.time()
    rows = []
    with Pool(4) as pool:
        specs = singles("K2") + singles("K4b")
        for r in pool.imap_unordered(run_spec, specs):
            r["id"] = r["spec"].get("id")
            rows.append(r)
        print("singles done", len(rows), round(time.time() - t0), flush=True)
        OUT.write_text(json.dumps(rows, indent=1, default=float), encoding="utf-8", newline="\n")
        d5 = {}
        for name in ("K2", "K4b"):
            d5[name] = d5_parallel(pool, name, D.LAM_BASE[name])
        specs = []
        for name in ("K2", "K4b"):
            dp, dt, ep, _ = d5[name]
            lam0 = D.LAM_BASE[name]
            for l in (lam0 * 0.6, lam0, lam0 * 1.5, lam0 * 2.2, lam0 * 3.3):
                specs.append({"K": name, "id": "D5", "lam": l, "dphi": dp, "dtheta": dt, "eps": ep})
            # combinations (D5 parameters carried)
            for rc in (50.0, 65.0, 80.0, 95.0):
                for l in (lam0 * 0.6, lam0, lam0 * 1.5):
                    specs.append({"K": name, "id": "C_D6_D5", "D6": rc, "lam": l, "dphi": dp, "dtheta": dt, "eps": ep})
                    specs.append({"K": name, "id": "C_D6_D2", "D6": rc, "D2": "2x2", "lam": l})
                    specs.append({"K": name, "id": "C_D6_D2_D5", "D6": rc, "D2": "2x2", "lam": l, "dphi": dp, "dtheta": dt, "eps": ep})
                    specs.append({"K": name, "id": "C_D6_D3", "D6": rc, "D3": "one", "lam": l})
                    specs.append({"K": name, "id": "C_D6_D3_D2_D5", "D6": rc, "D3": "one", "D2": "2x2", "lam": l, "dphi": dp, "dtheta": dt, "eps": ep})
            for f in (0.5, 0.3, 0.2):
                specs.append({"K": name, "id": "C_D1_D2", "D1": f, "D2": "2x2", "lam1": lam0 * 0.4})
                specs.append({"K": name, "id": "C_D1_D5", "D1": f, "lam1": lam0 * 0.4, "dphi": dp, "dtheta": dt, "eps": ep})
            for m in ("1x2x", "1x2y", "2x2"):
                for l in (lam0 * 1.5, lam0 * 2.2):
                    specs.append({"K": name, "id": "C_D2_D3_D5", "D2": m, "D3": "one", "lam": l, "dphi": dp, "dtheta": dt, "eps": ep})
        for r in pool.imap_unordered(run_spec, specs):
            r["id"] = r["spec"].get("id")
            rows.append(r)
            OUT.write_text(json.dumps(rows, indent=1, default=float), encoding="utf-8", newline="\n")
        print("all done", len(rows), round(time.time() - t0), flush=True)
    OUT.write_text(json.dumps({"rows": rows, "d5": d5}, indent=1, default=float), encoding="utf-8", newline="\n")
