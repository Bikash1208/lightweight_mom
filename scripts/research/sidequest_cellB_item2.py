"""Item 2 analysis of candidate B (RO4003C 1.524 mm, p 6.4, square patch gap 0.8 mm, via + diode on the back).  MODEL_ONLY.
    python sidequest_cellB_item2.py  -> results/sidequest_hybrid/cellB_item2.json (+ prints the tables)"""
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

H = 1.524e-3
ER, TD = K.LAMINATES["RO4003C"]
RES = ROOT / "results/sidequest_hybrid"
LSER = 1.6e-9  # chosen knob (series inductance of the back trace), H


def rows_of(name):
    p = RES / name
    return json.loads(p.read_text()) if p.exists() else []


def group(rows, key):
    d = {}
    for r in rows:
        d.setdefault(key(r), []).append(r)
    return d


def lerp_patch(rows0, rows1, t):
    out = []
    for a, b in zip(sorted(rows0, key=lambda r: r["kr"]), sorted(rows1, key=lambda r: r["kr"])):
        y0, y1 = complex(*a["Yuu_S"]), complex(*b["Yuu_S"])
        y = y0 + t * (y1 - y0)
        out.append(dict(a, Yuu_S=[y.real, y.imag]))
    return M.Patch(out)


def ev(p, er=ER, tand=TD, h=H, **kw):
    kw.setdefault("L_ser", LSER)
    return M.evaluate(p, "RO4003C", er, tand, h, **kw)


def brief(r):
    return {k: float(r[k]) for k in ("beta_ON", "beta_OFF", "beta_mean", "A0_rad", "alpha_loss_ON", "alpha_loss_OFF", "dY_mS") if k in r} if r.get("ok") else {"ok": False}


def main():
    rows = rows_of("cellB_ypatch_final_m32_m64.json") + rows_of("cellB_ypatch_final_m96.json")
    G = group(rows, lambda r: (r["mesh"], r["psi"], round(r["gap"], 4)))
    P = {k: M.Patch(v) for k, v in G.items()}
    out = {"label": "MODEL_ONLY: MoM Bloch patch grid + lumped via/diode branch; via not in the MoM; NOT HFSS-checked", "design": {"laminate": "RO4003C", "er": ER, "tand": TD, "h_mm": 1.524, "p_mm": 6.4, "gap_mm": 0.8, "patch_mm": 5.6, "L_ser_nH": LSER * 1e9, "L_via_nH": M.MU0 * H * 1e9, "diode": "MA4AGP907_measured"}}
    # 1) convergence + isotropy
    conv, iso = {}, {}
    for (mesh, psi, gap), p in sorted(P.items()):
        if abs(gap - 0.8) < 1e-6:
            conv[f"mesh{mesh}_psi{int(psi)}"] = brief(ev(p))
    out["convergence_isotropy"] = conv
    for k, v in conv.items():
        print(k, {a: round(b, 4) for a, b in v.items()})
    base = P[(64, 0.0, 0.8)]
    base45 = P[(64, 45.0, 0.8)]
    # 2) cross-pol (patch tensor)
    out["cross_pol_max_Yuv_over_Yuu"] = {f"mesh{k[0]}_psi{int(k[1])}": p.yuv_max for k, p in P.items() if abs(k[2] - 0.8) < 1e-6}
    # 3) tolerance (mesh 64, psi 0 and 45 where available)
    tol = {}
    nom = ev(base)
    out["nominal_m64_psi0"] = brief(nom)
    gaprows = rows_of("cellB_ypatch_gap_m64.json")
    gg = group(gaprows, lambda r: round(r["gap"], 4))
    r08 = G[(64, 0.0, 0.8)]
    if 0.6 in gg and 1.0 in gg:
        for name, dg in (("etch_-0.03_each_edge(gap+0.06)", +0.06), ("etch_+0.03_each_edge(gap-0.06)", -0.06), ("geom_+0.05_gap", +0.05), ("geom_-0.05_gap", -0.05)):
            t = abs(dg) / 0.2
            other = gg[1.0] if dg > 0 else gg[0.6]
            tol[name] = brief(ev(lerp_patch(r08, other, t)))
    er = group(rows_of("cellB_ypatch_er.json"), lambda r: r["lam"]) if (RES / "cellB_ypatch_er.json").exists() else {}
    for lam, e in (("RO4003C_lo", 3.33), ("RO4003C_hi", 3.43)):
        if lam in er:
            tol[f"er_{e}"] = brief(ev(M.Patch(er[lam]), er=e))
    for name, kw in (("diode_R_-20%", dict(dR=0.8)), ("diode_R_+20%", dict(dR=1.2)), ("diode_C_-20%", dict(dC=0.8)), ("diode_C_+20%", dict(dC=1.2)),
                     ("via_L_-8%(drill/plating bound)", dict(dL=0.92)), ("via_L_+8%(drill/plating bound)", dict(dL=1.08)),
                     ("antipad_C_ap_+20fF", dict(C_ap=20e-15)), ("kappa_0.7(via L model)", dict(kappa=0.7))):
        tol[name] = brief(ev(base, **kw))
    for dh in (-0.05e-3, 0.05e-3):
        tol[f"h_{dh * 1e3:+.2f}mm(patch Y held)"] = brief(ev(base, h=H + dh))
    out["tolerance"] = tol
    for k, v in tol.items():
        if v.get("ok", True) is False:
            print(k, "NO ROOT"); continue
        print("%-34s dbON %+.4f dbOFF %+.4f dA0 %+.1f%% aloss %.2f/%.2f" % (k, v["beta_ON"] - nom["beta_ON"], v["beta_OFF"] - nom["beta_OFF"], 100 * (v["A0_rad"] / nom["A0_rad"] - 1), v["alpha_loss_ON"], v["alpha_loss_OFF"]))
    # 4) other diode
    dio = {}
    for d in ("MA4AGP907_measured", "MA4AGFCP910_datasheet_typ", "MA4AGFCP910_datasheet_max"):
        dio[d] = brief(ev(base, diode=d))
    out["diodes"] = dio
    for k, v in dio.items():
        print(k, {a: round(b, 3) for a, b in v.items()})
    # 5) static cells: hard-wired via, series L chosen so that beta = beta_mean of the switched pair
    bm = nom["beta_mean"]
    def bstat(Lst):
        r = M.solve_beta(base, ER, TD, H, 1.0 / (1j * M.W0 * Lst))
        return None if r is None else r[0]
    st = {}
    for Lst in np.arange(1.5e-9, 12e-9, 0.5e-9):
        b = bstat(Lst)
        st[f"{Lst * 1e9:.1f}nH"] = b
    out["static_via_L_scan_beta"] = st
    try:
        Ls = brentq(lambda L: (bstat(L) or 0) - bm, 1.5e-9, 11e-9)
        out["static_cell"] = {"L_total_nH": Ls * 1e9, "L_ser_beyond_via_nH": Ls * 1e9 - M.MU0 * H * 1e9, "beta_target_mean": bm}
    except Exception as e:
        out["static_cell"] = {"error": str(e), "beta_target_mean": bm}
    print("static", out["static_cell"])
    (RES / "cellB_item2.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
    return out


if __name__ == "__main__":
    main()
