"""K = 2 and K = 4 runs (RT5880 chosen cell): option stages, epsilon choice, static-vs-SLL trade-off, comparators (dedicated K = 1 boards -> 4-feed + SP4T estimate, full 1-bit board)."""
import itertools
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_coded_k4 as K

OUT = Path("results/sidequest_hybrid/coded_k4_" + "_".join(sys.argv[1:] or ["all"]) + ".json")
LOG = {}
STAT, PAIRS = K.library()
VS_BIT = np.array([-1.0, 1.0])
VS_AN = np.linspace(STAT.min(), STAT.max(), 9)
G1 = [PAIRS[(2.8, 0.8)]]
G2 = [PAIRS[(2.8, 0.8)], PAIRS[(2.6, 0.8)]]
LIMIT = -12.7


def save():
    OUT.write_text(json.dumps(LOG, indent=1, default=float), encoding="utf-8", newline="\n")


def go(name, label, bd, opt, eps, target=None, lam=None, mu=0.0, gs=(1.0, 2.0)):
    t = time.time()
    best = None
    for g in gs:
        r, aux = K.run(bd, opt, eps, lam=lam, target=target, mu=mu, g=g)
        gm = np.mean([x["gain"] for x in r["states"]])
        ok = max(x["SLL"] for x in r["states"]) <= LIMIT
        if best is None or (ok, gm) > (best[0], best[1]):
            best = (ok, gm, r)
    LOG.setdefault(name, {})[label] = K.summary(best[2])
    save()
    print(name, label, LOG[name][label], round(time.time() - t), flush=True)
    return best[2]


if __name__ == "__main__":
    which = sys.argv[1:] or ["K2", "K4a", "K4b"]
    for name in which:
        bd = K.KBoard(*K.SETS[name])
        Kn = bd.K
        opt_bit = K.options(Kn, VS_BIT, G1)
        # epsilon: {0, pi}^(K-1) on the bit baseline (diodes <= 766 for K = 2; K >= 4: lambda of the K = 2 bit design)
        lam0 = 0.555
        best = None
        for e in itertools.product((0.0, np.pi), repeat=Kn - 1):
            eps = [0.0, *e]
            r, _ = K.run(bd, opt_bit, eps, lam=lam0, g=1.0, iters=2)
            key = (r["static_pct"], np.mean([x["gain"] for x in r["states"]]))
            if best is None or key > best[0]:
                best = (key, eps)
        eps = best[1]
        LOG.setdefault(name, {})["eps_deg"] = [round(float(np.degrees(x))) for x in eps]
        print(name, "eps", LOG[name]["eps_deg"], flush=True)
        target = 766 if Kn == 2 else None
        lamk = None
        r = go(name, "bit", bd, opt_bit, eps, target=target, lam=None if target else lam0)
        lam_fixed = None if target else lam0
        go(name, "I1", bd, K.options(Kn, VS_AN, G1), eps, target=target, lam=lam_fixed)
        go(name, "I1I2", bd, K.options(Kn, VS_AN, G2), eps, target=target, lam=lam_fixed)
        rfin = go(name, "I1I2I5", bd, K.options(Kn, VS_AN, G2), eps, target=target, lam=lam_fixed, mu=0.05)
        if Kn == 4:
            for lam in (0.2, 0.4, 0.8, 1.2, 2.0):
                go(name, f"tradeoff_lam{lam}", bd, K.options(Kn, VS_AN, G2), eps, lam=lam, mu=0.05, gs=(2.0,))
        # full 1-bit board: every element switched
        go(name, "full_1bit", bd, K.options(Kn, VS_BIT, G1, force_switched=True, allw=True), eps, lam=0.0, gs=(1.0,))
    # dedicated single-state boards (no switching): K = 1 per state, gives the ceiling and the 4-feed segmented estimate
    for name in ("K4a", "K4b"):
        ph, pol = K.SETS[name]
        for k in range(4):
            bd1 = K.KBoard([ph[k]], [pol[k]])
            opt = K.options(1, VS_AN, []) if False else None
    Path("/dev/null")
