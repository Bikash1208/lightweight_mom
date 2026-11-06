"""Item 2 final cell analysis, refined cell D2 (RO4003C 1.524, p 6.4, gap 0.64).  MODEL_ONLY.   python sidequest_cellB_final.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import brentq

import sidequest_cellB_core as K
import sidequest_cellB_model as M

RES = ROOT / "results/sidequest_hybrid"
H = 1.524e-3
ER, TD = K.LAMINATES["RO4003C"]


def rows(name):
    p = RES / name
    return json.loads(p.read_text()) if p.exists() else []


def grp(rs, key):
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


def brief(r):
    return {k: float(r[k]) for k in ("beta_ON", "beta_OFF", "beta_mean", "A0_rad", "alpha_loss_ON", "alpha_loss_OFF", "dY_mS", "ytot_rad_mS") if k in r} if r.get("ok") else {"ok": False}


def fscan(patch, L_ser, f_list, state, C_ap=0.0):
    """TM0 beta/k0 vs frequency (patch Y scaled with f, capacitive) for one state."""
    out = []
    for f in f_list:
        w = 2 * math.pi * f
        k0 = w / K.C0
        sc = f / M.F0

        def F(kr):
            ec = ER * (1 - 1j * TD)
            kz = k0 * np.sqrt(ec - kr * kr + 0j)
            ys = -1j * (w * M.EPS0 * ec / kz) / np.tan(kz * H)
            ya = 1j * w * M.EPS0 / (k0 * math.sqrt(kr * kr - 1))
            R, C = M.DIODE["MA4AGP907_measured"]
            L = M.MU0 * H + L_ser
            yb = 1 / (R + 1j * w * L) if state == "ON" else 1 / (1j * w * L + 1 / (1j * w * (C + C_ap)))
            return ya + ys + patch.y_sw(kr) * sc + yb

        ks = np.linspace(1.002, 1.48, 500)
        v = np.array([F(k).imag for k in ks])
        ix = np.where(np.diff(np.sign(v)) != 0)[0]
        out.append((f / 1e9, None if len(ix) == 0 else float(brentq(lambda k: F(k).imag, ks[ix[0]], ks[ix[0] + 1]))))
    return out


def main():
    A = rows("cellB_ypatch_d2a.json") + rows("cellB_ypatch_d2b.json")
    G = grp([r for r in A if r["lam"] == "RO4003C"], lambda r: (r["mesh"], r["psi"], round(r["gap"], 4)))
    P = {k: M.Patch(v) for k, v in G.items()}
    out = {"label": "MODEL_ONLY (MoM Bloch patch grid + lumped via/diode branch; via not in the MoM)", "cell": {"laminate": "RO4003C", "er": ER, "tand": TD, "h_mm": 1.524, "p_mm": 6.4, "gap_mm": 0.64, "patch_mm": 5.76}}
    # knob scan on mesh 60, psi 0
    base = P[(60, 0.0, 0.64)]
    scan = []
    for Ls in np.arange(0.0, 3.01e-9, 0.25e-9):
        r = M.evaluate(base, "x", ER, TD, H, L_ser=float(Ls))
        scan.append({"L_ser_nH": Ls * 1e9, **brief(r)})
    out["Lser_scan_m60_psi0"] = scan
    for s in scan:
        print("Lser %.2f" % s["L_ser_nH"], {k: round(v, 3) for k, v in s.items() if k != "L_ser_nH"})
    # choose L_ser: the largest A0 with beta_ON >= 1.10 (else the one with beta_ON closest)
    ok = [s for s in scan if "beta_ON" in s and s["beta_ON"] >= 1.10 and s["A0_rad"] >= 12]
    LS = (max(ok, key=lambda s: s["A0_rad"])["L_ser_nH"] if ok else max([s for s in scan if "beta_ON" in s], key=lambda s: s["beta_ON"])["L_ser_nH"]) * 1e-9
    out["chosen_L_ser_nH"] = LS * 1e9
    ev = lambda p, **kw: brief(M.evaluate(p, "x", kw.pop("er", ER), TD, kw.pop("h", H), L_ser=kw.pop("L_ser", LS), **kw))
    # convergence + isotropy
    conv = {}
    for (mesh, psi, gap), p in sorted(P.items()):
        if abs(gap - 0.64) < 1e-6:
            conv[f"mesh{mesh}_psi{int(psi)}"] = ev(p)
    out["convergence_isotropy"] = conv
    for k, v in conv.items():
        print(k, {a: round(b, 4) for a, b in v.items()})
    out["cross_pol_Yuv_over_Yuu_max"] = {f"mesh{k[0]}_psi{int(k[1])}": p.yuv_max for k, p in P.items() if abs(k[2] - 0.64) < 1e-6}
    nom = ev(base)
    out["nominal"] = nom
    # tolerance
    tol = {}
    g0 = G[(60, 0.0, 0.64)]
    lo = G.get((60, 0.0, round(0.64 - 2 * 6.4 / 60, 4)))
    hi = G.get((60, 0.0, round(0.64 + 2 * 6.4 / 60, 4)))
    if lo and hi:
        for nm, dg in (("etch +0.03/edge (gap +0.06)", 0.06), ("etch -0.03/edge (gap -0.06)", -0.06), ("geom gap +0.05", 0.05), ("geom gap -0.05", -0.05)):
            t = abs(dg) / (2 * 6.4 / 60)
            tol[nm] = ev(lerp(g0, hi if dg > 0 else lo, t))
    er = grp([r for r in A if r["lam"] in ("RO4003C_lo", "RO4003C_hi")], lambda r: r["lam"])
    for lam, e in (("RO4003C_lo", 3.33), ("RO4003C_hi", 3.43)):
        if lam in er:
            tol[f"er {e}"] = ev(M.Patch(er[lam]), er=e)
    for nm, kw in (("diode R -20%", dict(dR=0.8)), ("diode R +20%", dict(dR=1.2)), ("diode C -20%", dict(dC=0.8)), ("diode C +20%", dict(dC=1.2)),
                   ("via L -8% (drill/plating bound)", dict(dL=0.92)), ("via L +8% (drill/plating bound)", dict(dL=1.08)),
                   ("L_ser -15% (trace etch)", dict(L_ser=LS * 0.85)), ("L_ser +15% (trace etch)", dict(L_ser=LS * 1.15)),
                   ("pad C_ap = +5 fF (20% of an assumed 25 fF)", dict(C_ap=5e-15)), ("pad C_ap = -5 fF equivalent (C_ap 0, diode C -5 fF)", dict(dC=1 - 5 / 28.16)),
                   ("kappa 0.7 (via L model uncertainty)", dict(kappa=0.7)), ("h -0.05 mm (patch Y held)", dict(h=H - 0.05e-3)), ("h +0.05 mm (patch Y held)", dict(h=H + 0.05e-3))):
        tol[nm] = ev(base, **kw)
    out["tolerance"] = tol
    dev = {}
    for k, v in tol.items():
        if "beta_ON" not in v:
            dev[k] = None
            print(k, "NO ROOT")
            continue
        dev[k] = {"dbON": v["beta_ON"] - nom["beta_ON"], "dbOFF": v["beta_OFF"] - nom["beta_OFF"], "dA0_pct": 100 * (v["A0_rad"] / nom["A0_rad"] - 1), "beta_ON": v["beta_ON"]}
        print("%-52s dbON %+.4f dbOFF %+.4f dA0 %+.1f%% bON %.3f" % (k, dev[k]["dbON"], dev[k]["dbOFF"], dev[k]["dA0_pct"], v["beta_ON"]))
    out["tolerance_dev"] = dev
    # combined worst case (linear superposition of the spec factors: etch 0.03, er, diode R, diode C, via L 8 %, L_ser 15 %)
    spec = ["etch +0.03/edge (gap +0.06)", "etch -0.03/edge (gap -0.06)", "er 3.33", "er 3.43", "diode R -20%", "diode R +20%", "diode C -20%", "diode C +20%", "via L -8% (drill/plating bound)", "via L +8% (drill/plating bound)", "L_ser -15% (trace etch)", "L_ser +15% (trace etch)"]
    pairs = [(spec[i], spec[i + 1]) for i in range(0, len(spec), 2)]
    comb = {"dbON_worst": 0.0, "dbOFF_worst": 0.0, "dbON_rss": 0.0, "dbOFF_rss": 0.0, "dA0_worst_pct": 0.0, "dA0_rss_pct": 0.0, "beta_ON_min": nom["beta_ON"]}
    for a, b in pairs:
        if dev.get(a) is None or dev.get(b) is None:
            continue
        mxON = max(abs(dev[a]["dbON"]), abs(dev[b]["dbON"]))
        mxOFF = max(abs(dev[a]["dbOFF"]), abs(dev[b]["dbOFF"]))
        mxA = max(abs(dev[a]["dA0_pct"]), abs(dev[b]["dA0_pct"]))
        comb["dbON_worst"] += mxON
        comb["dbOFF_worst"] += mxOFF
        comb["dA0_worst_pct"] += mxA
        comb["dbON_rss"] += mxON ** 2
        comb["dbOFF_rss"] += mxOFF ** 2
        comb["dA0_rss_pct"] += mxA ** 2
        comb["beta_ON_min"] -= min(0, dev[a]["dbON"], dev[b]["dbON"]) * -1
    comb["dbON_rss"] = math.sqrt(comb["dbON_rss"])
    comb["dbOFF_rss"] = math.sqrt(comb["dbOFF_rss"])
    comb["dA0_rss_pct"] = math.sqrt(comb["dA0_rss_pct"])
    out["tolerance_combined"] = comb
    print("combined", {k: round(v, 4) for k, v in comb.items()})
    # frequency scan: resonance distance (TM0 root existence and beta vs f), both states
    fl = np.arange(7.0e9, 14.01e9, 0.5e9)
    out["fscan"] = {st: fscan(base, LS, fl, st) for st in ("ON", "OFF")}
    for st in ("ON", "OFF"):
        print(st, [(f, None if b is None else round(b, 3)) for f, b in out["fscan"][st]])
    # normal-incidence cell resonance (series/parallel zero of the total admittance seen by a plane wave), ON and OFF
    def ni(f, state):
        w = 2 * math.pi * f
        k0 = w / K.C0
        ec = ER * (1 - 1j * TD)
        kz = k0 * np.sqrt(ec)
        ys = -1j * (w * M.EPS0 * ec / kz) / np.tan(kz * H)
        yu = w * M.EPS0 / k0
        R, C = M.DIODE["MA4AGP907_measured"]
        L = M.MU0 * H + LS
        yb = 1 / (R + 1j * w * L) if state == "ON" else 1 / (1j * w * L + 1 / (1j * w * C))
        return (yu + ys + base.y_rad * f / M.F0 + yb).imag

    res_ni = {}
    for st in ("ON", "OFF"):
        fs = np.linspace(4e9, 20e9, 1600)
        v = np.array([ni(f, st) for f in fs])
        ix = np.where(np.diff(np.sign(v)) != 0)[0]
        res_ni[st] = [float(brentq(lambda f: ni(f, st), fs[i], fs[i + 1])) / 1e9 for i in ix]
    out["normal_incidence_resonances_GHz"] = res_ni
    print("normal-incidence resonances (Im Y_total = 0), GHz:", res_ni)
    (RES / "cellB_item2_final.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
