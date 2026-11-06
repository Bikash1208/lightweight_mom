"""Final cell numbers of the via-switched GB cell (RT5880 1.575, p 6.4, patch gap 0.4, MA4AGP907 measured chip PIN, series C_s knob): convergence, isotropy, step, beta ON / OFF / mean, frequency margins, tolerance.  MODEL_ONLY.
    python sidequest_via_final.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import brentq

import sidequest_cellB_model as M
import sidequest_via_cell as V
import sidequest_via_design as D

RES = ROOT / "results/sidequest_hybrid"
DIODE = "MA4AGP907 (measured)"


def rows(name):
    p = RES / name
    return json.loads(p.read_text()) if p.exists() else []


def group(rs, key):
    d = {}
    for r in rs:
        d.setdefault(key(r), []).append(r)
    return d


def lerp(r0, r1, t):
    out = []
    for a, b in zip(sorted(r0, key=lambda r: r["kr"]), sorted(r1, key=lambda r: r["kr"])):
        y0, y1 = complex(*a["Yuu_S"]), complex(*b["Yuu_S"])
        y = y0 + t * (y1 - y0)
        out.append(dict(a, Yuu_S=[y.real, y.imag]))
    return M.Patch(out)


def levels(P, C_s, L_ser, **kw):
    return D.br("ON", DIODE, C_s, L_ser, **kw), D.br("OFF", DIODE, C_s, L_ser, **kw)


def betas(P, C_s, L_ser, er=V.ER, **kw):
    yon, yoff = levels(P, C_s, L_ser, **kw)
    return {"ON": V.root(P, yon, er=er), "OFF": V.root(P, yoff, er=er), "mean": V.root(P, 0.5 * (yon + yoff), er=er), "step_mS": (yoff - yon).imag * 1e3,
            "sheet_on_mS": (P.y_sw(1.062) + yon).imag * 1e3, "sheet_off_mS": (P.y_sw(1.062) + yoff).imag * 1e3}


def knob_for_mean(P, L_ser=0.0, beta_mean=1.062):
    def f(C_pF):
        yon, yoff = levels(P, C_pF * 1e-12, L_ser)
        r = V.root(P, 0.5 * (yon + yoff))
        return (1.0 if r is None else r[0]) - beta_mean

    return brentq(f, 0.2, 3.0, xtol=1e-6) * 1e-12


def fscan(P, C_s, L_ser, fl):
    out = {}
    for st in ("ON", "OFF", "mean"):
        res = []
        for f in fl:
            w = 2 * math.pi * f
            sc = f / V.F0
            yon, yoff = levels(P, C_s, L_ser, f=f)
            yb = {"ON": yon, "OFF": yoff, "mean": 0.5 * (yon + yoff)}[st]
            r = V.root(P, yb, f=f, scale=sc)
            res.append((round(f / 1e9, 2), None if r is None else round(r[0], 3)))
        out[st] = res
    return out


def main():
    A = rows("cellB_rt5880_a.json") + rows("cellB_rt5880_b.json")
    if not A:
        A = rows("cellB_rt5880_explore.json")
    G = group([r for r in A if r["lam"] == "RT5880"], lambda r: (r["mesh"], r["psi"], round(r["gap"], 3)))
    P = {k: M.Patch(v) for k, v in G.items()}
    key0 = (64, 0.0, 0.4) if (64, 0.0, 0.4) in P else sorted(P)[0]
    base = P[key0]
    C_s = knob_for_mean(base)
    L_ser = 0.0
    out = {"label": "MODEL_ONLY: MoM patch grid + lumped via/diode/series-C branch", "cell": {"laminate": "RT5880 eps_r 2.2 tan d 0.0009", "h_mm": 1.575, "p_mm": 6.4, "gap_mm": 0.4, "patch_mm": 6.0, "diode": DIODE, "C_s_pF": C_s * 1e12, "L_ser_nH": L_ser * 1e9,
           "L_via_nH": M.MU0 * V.H * 1e9, "C_pad_fF_assumed": 25.0}}
    out["nominal"] = betas(base, C_s, L_ser)
    print("C_s %.3f pF" % (C_s * 1e12), {k: (v if not isinstance(v, tuple) else tuple(round(x, 3) for x in v)) for k, v in out["nominal"].items()})
    conv = {}
    for k, p in sorted(P.items()):
        if abs(k[2] - 0.4) < 1e-6:
            b = betas(p, C_s, L_ser)
            conv[f"mesh{k[0]}_psi{int(k[1])}"] = {kk: (None if v is None else (float(v[0]) if isinstance(v, tuple) else float(v))) for kk, v in b.items()}
    out["convergence_isotropy"] = conv
    for k, v in conv.items():
        print(k, {a: (None if b is None else round(b, 4)) for a, b in v.items()})
    # tolerance (mesh 64 psi 0)
    tol = {}
    nom = betas(base, C_s, L_ser)

    def rec(name, b):
        tol[name] = {"beta_mean": None if b["mean"] is None else b["mean"][0], "beta_OFF": None if b["OFF"] is None else b["OFF"][0], "beta_ON": None if b["ON"] is None else b["ON"][0], "step_mS": b["step_mS"], "sheet_on_mS": b["sheet_on_mS"], "sheet_off_mS": b["sheet_off_mS"]}

    for nm, kw in (("diode C -20%", dict(dC=0.8)), ("diode C +20%", dict(dC=1.2)), ("diode R -20%", dict(dR=0.8)), ("diode R +20%", dict(dR=1.2)), ("drill: via L -8%", dict(dL=0.92)), ("drill: via L +8%", dict(dL=1.08)),
                   ("C_s -10%", dict(dCs=0.9)), ("C_s +10%", dict(dCs=1.1)), ("pad C -5 fF", dict(C_pad=20e-15)), ("pad C +5 fF", dict(C_pad=30e-15)), ("kappa 0.7 (via L model)", dict(kappa=0.7))):
        rec(nm, betas(base, C_s, L_ser, **kw))
    g_lo, g_hi = G.get((64, 0.0, 0.2)), G.get((64, 0.0, 0.6))
    if g_lo and g_hi:
        for nm, sg in (("etch +0.03/edge (gap +0.06)", +0.06), ("etch -0.03/edge (gap -0.06)", -0.06)):
            Pt = lerp(G[key0], g_hi if sg > 0 else g_lo, abs(sg) / 0.2)
            rec(nm, betas(Pt, C_s, L_ser))
    for lam, e in (("RT5880_lo", 2.18), ("RT5880_hi", 2.22)):
        rr = [r for r in A if r["lam"] == lam]
        if rr:
            rec(f"eps_r {e} (assumed +-0.02)", betas(M.Patch(rr), C_s, L_ser, er=e))
    out["tolerance"] = tol
    for k, v in tol.items():
        print("%-30s" % k, {a: (None if b is None else round(b, 3)) for a, b in v.items()})
    # frequency margins
    fl = np.arange(7.0e9, 14.01e9, 0.5e9)
    out["fscan"] = fscan(base, C_s, L_ser, fl)
    for st in ("ON", "OFF", "mean"):
        print(st, out["fscan"][st])
    fo, ff = D.resonances(DIODE, C_s, L_ser)
    out["branch_series_resonances_GHz"] = {"ON": fo, "OFF": ff}
    print("branch series resonances GHz: ON %.2f OFF %.2f" % (fo, ff))
    (RES / "via_final.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
