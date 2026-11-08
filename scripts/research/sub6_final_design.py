"""SUB6 FINAL BOARD: the chosen calibrated V3 design (200 mm, N 240, placement a, taper 0.5, ring exclusion on) re-generated with its rotations (MODEL_ONLY)."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_final_cal as FC
import sub6_board_island2 as X
import sub6_board_model4 as R4

CONFIG = dict(side=0.2, n=240, place="a", taper=0.5, excl=1)


def final_design(mode="cal", config=None):
    c = config or CONFIG
    FC.MODE = mode
    T = R4.tables()
    D = R4.design(c["side"], c["n"], c["place"], c["taper"], c["excl"], T)
    sa = R4.site_arrays(D)
    ph = R4.physics(D, sa)
    dirs = {w: X.directions_w(D, ph, w) for w in X.W_GRID}
    single = {}
    for w in X.W_GRID:
        fl, info = R4.optimise_flips(D, sa, ph, dirs[w])
        single[w] = (fl, X.metrics_state(D, sa, ph, dirs[w], fl, FC.LAUNCH), info)
    wstar = max(X.W_GRID, key=lambda w: single[w][1]["minGA"])
    best = (wstar,) + single[wstar]
    best = (best[0], best[1], best[2], dict(best[3], ring_weight=0.0))
    for w in sorted({0.0, wstar}):
        for rwt in (0.1, 0.3, 1.0):
            fl3, info3 = X.optimise_flips_multi(D, sa, ph, dirs[w], ring_weight=rwt)
            m3 = X.metrics_state(D, sa, ph, dirs[w], fl3, FC.LAUNCH)
            if m3["margin"] > best[2]["margin"]:
                best = (w, fl3, m3, dict(info3, ring_weight=rwt))
    w, fl, m, info = best
    S = D["S"]
    phiR = np.arctan2(S.y, S.x)
    direction = dirs[w]
    chi = np.arctan2(direction[:, 1], direction[:, 0])
    axis = np.mod(0.5 * (chi + phiR) + fl * math.pi / 2, math.pi)  # slot long-axis angle (board frame), mod pi
    return dict(D=D, sa=sa, ph=ph, T=T, w=w, flips=fl, metrics=m, info=info, axis=axis, direction=direction)


if __name__ == "__main__":
    for mode in ("cal", "base"):
        r = final_design(mode)
        D = r["D"]
        live, mask = D["live"], D["mask"]
        print(mode, "sites", D["S"].n, "live", int(live.sum()), "switched live", int((mask & live).sum()), "rotated", int(r["flips"].sum()), "rotated switched", int((r["flips"] & mask).sum()), "w", r["w"], {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r["metrics"].items()})
