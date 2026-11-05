"""Minimum control-line count for K = 4 (CP set RHCP/LHCP/RHCP/LHCP, RT5880 chosen cell, I1+I2+I5 options, lambda 0.555, eps (0, 0, 0, 180)).
A line = one complement class of ON-patterns (with an inverter); the 14 non-constant ON-patterns form 7 classes. For each L the design is restricted to L allowed classes on both axes (the same set; lines shared by x and y diodes of
the same class: L total, or 2L if the axes are wired separately) and the best subset is searched: exhaustive for L <= 3, greedy backward elimination for L = 4..6.  Score = mean gain - 0.5 dB per dB of SLL above -12.7 dB.
Output results/sidequest_hybrid/coded_k4_lines.json."""
import itertools
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_coded_k4 as K

BD = None
STAT, PAIRS = K.library()
VS_AN = np.linspace(STAT.min(), STAT.max(), 9)
G2 = [PAIRS[(2.8, 0.8)], PAIRS[(2.6, 0.8)]]
EPS = [0.0, 0.0, 0.0, np.pi]
CLASSES = sorted({K.canon(w) for w in K.words(4)})  # 7


def init():
    global BD
    BD = K.KBoard(*K.SETS["K4b"])


def evaluate(cls):
    cls = set(cls)
    opt = K.options(4, VS_AN, G2, cls_x=cls, cls_y=cls)
    r, _ = K.run(BD, opt, EPS, lam=0.555, mu=0.05, g=2.0, iters=3)
    s = K.summary(r)
    gm = float(np.mean(s["gain"]))
    sll = max(s["SLL"])
    s["score"] = gm - 0.5 * max(0.0, sll + 12.7)
    s["classes"] = [list(c) for c in sorted(cls)]
    s["gain_mean"] = gm
    s["SLL_worst"] = sll
    return s


if __name__ == "__main__":
    out = {}
    with Pool(4, initializer=init) as pool:
        full = pool.map(evaluate, [CLASSES])[0]
        out["7"] = full
        print("L=7", full["gain_mean"], full["SLL_worst"], full["static_pct"], full["diodes"], flush=True)
        for L in (1, 2, 3):
            subs = list(itertools.combinations(CLASSES, L))
            res = pool.map(evaluate, subs)
            best = max(res, key=lambda r: r["score"])
            out[str(L)] = best
            out[f"{L}_all"] = [{"classes": r["classes"], "gain_mean": r["gain_mean"], "SLL_worst": r["SLL_worst"], "static_pct": r["static_pct"], "diodes": r["diodes"]} for r in res]
            print(f"L={L}", round(best["gain_mean"], 2), best["SLL_worst"], best["static_pct"], best["diodes"], best["classes"], flush=True)
        cur = list(CLASSES)
        for L in (6, 5, 4):
            cand = [tuple(c for c in cur if c != x) for x in cur]
            res = pool.map(evaluate, cand)
            best = max(res, key=lambda r: r["score"])
            out[str(L)] = best
            cur = [tuple(c) for c in best["classes"]]
            print(f"L={L}", round(best["gain_mean"], 2), best["SLL_worst"], best["static_pct"], best["diodes"], best["classes"], flush=True)
    Path("results/sidequest_hybrid/coded_k4_lines.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
