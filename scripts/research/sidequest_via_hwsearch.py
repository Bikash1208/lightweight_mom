"""Search over the bias-hardware parasitics and the knobs (C_s, L_ser): which hardware allows the 12-20 mS step at the tensor mean with margin >= 1.5 GHz?  MODEL_ONLY.   python sidequest_via_hwsearch.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_via_hw as W
import sidequest_via_cell as V
import sidequest_cellB_model as M

RES = ROOT / "results/sidequest_hybrid"


def main():
    P = W.patch()
    yp = P.y_sw(1.062).imag
    cs_grid = np.geomspace(0.1e-12, 6e-12, 70)
    out = []
    n_pass1 = 0
    for C_A in (10e-15, 15e-15, 25e-15, 35e-15):
        for C_p in (8e-15, 15e-15, 25e-15):
            for R_b in (1000.0, 2200.0, 4700.0, 10000.0, 22000.0):
                for L_x in (0.0, 0.5e-9, 1.0e-9, 1.5e-9):
                    rb = dict(W.RB, R_b=R_b, C_p=C_p, C_A=C_A, L_via_extra=0.05e-9 + L_x)
                    for C in cs_grid:
                        yon, yoff = W.branch2("ON", V.F0, C, rb), W.branch2("OFF", V.F0, C, rb)
                        step = (yoff - yon).imag * 1e3
                        mean = (yp + 0.5 * (yon + yoff).imag) * 1e3
                        if not (12.0 <= abs(step) <= 20.0) or abs(mean - 7.55) > 0.3 or yoff.real * 1e3 > 1.5:
                            continue
                        n_pass1 += 1
                        ok = True
                        for st in ("ON", "OFF"):
                            if any(abs(x - 10.5) < 1.5 for x in W.extrema(C, st, rb)):
                                ok = False
                                break
                        if not ok:
                            continue
                        dm, ds = [], []
                        for kw in (dict(dC=0.8), dict(dC=1.2), dict(dL=0.92), dict(dL=1.08), dict(dCs=0.95), dict(dCs=1.05), dict(dCA=0.7), dict(dCA=1.3), dict(dCp=0.7), dict(dCp=1.3), dict(dRb=0.8), dict(dRb=1.2), dict(dCV=0.7), dict(dCV=1.3)):
                            a, b = W.branch2("ON", V.F0, C, rb, **kw), W.branch2("OFF", V.F0, C, rb, **kw)
                            dm.append(abs(0.5 * ((a + b).imag - (yon + yoff).imag)) * 1e3)
                            ds.append(abs(abs((b - a).imag) - abs(step / 1e3)) * 1e3)
                        rm, rs = math.sqrt(sum(x * x for x in dm)), math.sqrt(sum(x * x for x in ds))
                        out.append(dict(C_A_fF=C_A * 1e15, C_p_fF=C_p * 1e15, R_b=R_b, L_x_nH=L_x * 1e9, C_s_pF=C * 1e12, step=step, mean=mean, Y_on=[yon.real * 1e3, yon.imag * 1e3], Y_off=[yoff.real * 1e3, yoff.imag * 1e3], rss_mean=rm, rss_step=rs, merit=rm / 1.7 + rs / (0.12 * abs(step))))
    print("pass-1 candidates", n_pass1, "with margin", len(out))
    out.sort(key=lambda r: r["merit"])
    for r in out[:12]:
        print({k: (round(v, 3) if isinstance(v, float) else (np.round(v, 2).tolist() if isinstance(v, list) else v)) for k, v in r.items()})
    (RES / "via_hw_search.json").write_text(json.dumps(out[:200], indent=1, default=float), encoding="utf-8", newline="\n")
    # which parasitic combinations are feasible at all
    feas = {}
    for r in out:
        feas.setdefault((r["C_A_fF"], r["C_p_fF"]), 0)
        feas[(r["C_A_fF"], r["C_p_fF"])] += 1
    print("feasible settings by (C_A fF, C_p fF):", feas)


if __name__ == "__main__":
    main()
