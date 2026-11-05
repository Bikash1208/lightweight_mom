"""L = 2 classes per axis with every epsilon in {0, pi}^3 (the class population depends on epsilon): exhaustive over class pairs and epsilon, same options/score as sidequest_k4_lines.py."""
import itertools
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_coded_k4 as K
import sidequest_k4_lines as L

BD = None


def init():
    global BD
    BD = K.KBoard(*K.SETS["K4b"])


def ev(a):
    cls, e = a
    opt = K.options(4, L.VS_AN, L.G2, cls_x=set(cls), cls_y=set(cls))
    r, _ = K.run(BD, opt, [0.0, *[np.pi * b for b in e]], lam=0.555, mu=0.05, g=2.0, iters=3)
    s = K.summary(r)
    gm = float(np.mean(s["gain"]))
    return {"classes": [list(c) for c in cls], "eps_bits": list(e), "gain_mean": gm, "SLL_worst": max(s["SLL"]), "static_pct": s["static_pct"], "diodes": s["diodes"]}


if __name__ == "__main__":
    jobs = [(c, e) for c in itertools.combinations(L.CLASSES, 2) for e in itertools.product((0, 1), repeat=3)]
    with Pool(4, initializer=init) as pool:
        res = pool.map(ev, jobs)
    ok = [r for r in res if r["SLL_worst"] <= -12.7]
    best = max(res, key=lambda r: r["gain_mean"] - 0.5 * max(0, r["SLL_worst"] + 12.7))
    print("feasible (SLL <= -12.7):", len(ok), "of", len(res), "best:", best, flush=True)
    Path("results/sidequest_hybrid/coded_k4_lines_eps.json").write_text(json.dumps({"all": res, "best": best}, indent=1, default=float), encoding="utf-8", newline="\n")
