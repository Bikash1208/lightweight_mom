"""CP3 / CP4 cell package, RO3006 1.27 mm, disc r 2.30 (lattice 6.4): (1) switched cell (rim 0.2, slit half-length 2.10, sheets DX / DY 0.4 x 1.0) OFF / ON at 90 and 110 deg along x, (2) neutral static cell (open slit half-length 1.951) at 110 deg.
Predictions MODEL_ONLY: loaded Bloch tensor (MoM mesh 64) + transverse resonance, with the 30 mm PEC lid (CP3) and without the lid (radiation boundary top, CP4); Q from the susceptance slope and the loss conductance; alpha implied = omega / (2 Q v_g), v_g from the 90 -> 110 deg frequencies.
Calibration of the model: HFSS r 2.20 check at 90 deg (OFF 10.0193 GHz Q 422, ON 10.3617 GHz Q 222) is evaluated with the same code: offsets and a Q scale factor are stored.
    python sidequest_cell_r230_package.py geom | pred | docs"""
import json
import math
import shutil
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

import sidequest_cell_ro3006_package as P
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_cell_RO3006_r2p30_cp3"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\cell_RO3006_r2p30_cp3")
R = 2.30
LS_SW = 2.10
LS_ST = 1.951  # neutral static cell: 0.929 x (r - 0.2)
C0 = 299792458.0


def fgrid(ph):
    return np.arange(8.4, 13.01, 0.4) * 1e9 if ph > 100 else np.arange(8.0, 11.61, 0.4) * 1e9


def job(a):
    import warnings as w

    w.simplefilter("ignore")
    import sidequest_cell_search as C
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.materials import Substrate

    kind, ls, ph, f = a
    gc.SUB = Substrate("ro3006", P.ER, P.TAND, P.H * 1e-3)
    L.R = R
    beta = np.radians(ph) / (P.P * 1e-3)
    y = C.bloch_loaded(L.cell(ls, P.POS, kind), f, 64, beta, 0.0)
    return {"kind": kind, "ls": ls, "ph": ph, "f": f, "Y_re_im_mS": [float(y[0, 0].real * 1e3), float(y[0, 0].imag * 1e3)]}


def stack(beta, f, lid):
    """slab + air above (with the 30 mm PEC lid or open): complex admittance seen by the sheet [S]."""
    from src.dispersion import sheet_model as sm
    from src.materials import Substrate

    sub = Substrate("s", P.ER, P.TAND, P.H * 1e-3)
    om, k0 = 2 * np.pi * f, 2 * np.pi * f / C0
    ya = P.y_air_lid(beta, om, k0) if lid else sm._y_air(beta, om, k0)
    return ya + sm._y_slab(beta, om, sub.eps_complex, k0, sub.thickness)


def resonance(rows, ph, lid):
    """eigenfrequency, Q and G from one Bloch curve (rows: list of dicts sorted by f)."""
    beta = np.radians(ph) / (P.P * 1e-3)
    fs = np.array([r["f"] for r in rows])
    re = CubicSpline(fs, [r["Y_re_im_mS"][0] * 1e-3 for r in rows])
    im = CubicSpline(fs, [r["Y_re_im_mS"][1] * 1e-3 for r in rows])
    fll = C0 * beta / (2 * np.pi)

    def bt(f):
        return stack(beta, f, lid).imag + float(im(f))

    fine = np.linspace(fs[0], fs[-1], 3000)
    v = np.array([bt(f) for f in fine])
    ix = np.where(np.diff(np.sign(v)) != 0)[0]
    roots = [brentq(bt, fine[i], fine[i + 1], xtol=1e3) for i in ix if fine[i] < fll]
    if not roots:
        return None
    f0 = roots[0]
    h = 2e6
    dB = (bt(f0 + h) - bt(f0 - h)) / (2 * np.pi * 2 * h)
    G = float(re(f0)) + stack(beta, f0, lid).real
    Q = 2 * np.pi * f0 * dB / (2 * G) if G > 0 else None
    return {"f_GHz": f0 / 1e9, "Q": Q, "G_mS": G * 1e3, "G_cell_mS": float(re(f0)) * 1e3, "G_slab_mS": float(stack(beta, f0, lid).real) * 1e3, "roots_GHz": [r / 1e9 for r in roots]}


def predictions():
    cases = [("open", LS_SW, 90.0), ("ON", LS_SW, 90.0), ("open", LS_SW, 110.0), ("ON", LS_SW, 110.0), ("open", LS_ST, 90.0), ("open", LS_ST, 110.0)]
    jobs = [(k, ls, ph, f) for k, ls, ph in cases for f in fgrid(ph)]
    with ProcessPoolExecutor(40) as ex:
        res = list(ex.map(job, jobs))
    Path(ROOT / "results/sidequest_hybrid/k2_cell_r230_bloch_raw.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    out = {"cases": {}}
    for k, ls, ph in cases:
        rows = sorted([r for r in res if r["kind"] == k and abs(r["ls"] - ls) < 1e-9 and r["ph"] == ph], key=lambda r: r["f"])
        name = ("switched" if ls == LS_SW else "neutral") + "_" + ("OFF" if k == "open" and ls == LS_SW else "ON" if k == "ON" else "static") + "_%d" % int(ph)
        out["cases"][name] = {"lid": resonance(rows, ph, True), "radiation_top": resonance(rows, ph, False)}
    # r 2.20 HFSS reference with the same code (saved mesh-64 Bloch data)
    ref = json.loads((ROOT / "results/sidequest_hybrid/k2_offcal.json").read_text())["raw"]
    refc = {}
    for mode, nm, fh, qh in (("open", "OFF", 10.0193, 422.0), ("ON", "ON", 10.3617, 222.0)):
        rows = sorted([x for x in ref if x["mesh"] == 64 and x["mode"] == mode], key=lambda x: x["f"])
        rows = [{"f": x["f"] * 1e9, "Y_re_im_mS": x["Y_mS"]} for x in rows]
        rr = resonance(rows, 90.0, True)
        refc[nm] = {"model_f_GHz": rr["f_GHz"], "model_Q": rr["Q"], "HFSS_f_GHz": fh, "HFSS_Q": qh, "Q_scale_HFSS_over_model": None if not rr["Q"] else qh / rr["Q"], "df_HFSS_minus_model_GHz": fh - rr["f_GHz"], "G_mS": rr["G_mS"]}
    out["r2.20_reference_90deg_lid"] = refc
    # group velocity and alpha implied (switched cell, both states; neutral static)
    for st in ("switched_OFF", "switched_ON", "neutral_static"):
        a, b = out["cases"][st + "_90"], out["cases"][st + "_110"]
        for top in ("lid", "radiation_top"):
            if a[top] and b[top]:
                dw = 2 * np.pi * (b[top]["f_GHz"] - a[top]["f_GHz"]) * 1e9
                db = np.radians(110.0 - 90.0) / (P.P * 1e-3)
                vg = dw / db
                for ph, d in (("90", a[top]), ("110", b[top])):
                    if d["Q"]:
                        d["alpha_implied_Np_per_m"] = float(2 * np.pi * d["f_GHz"] * 1e9 / (2 * d["Q"] * vg))
                        d["vg_over_c"] = float(vg / C0)
    return out


def geom():
    OUT.mkdir(parents=True, exist_ok=True)
    allres = {}
    for tag, ls, sheets in (("switched", LS_SW, True), ("neutral", LS_ST, False)):
        P.R = R
        P.LS = ls
        parts = P.build()
        if not sheets:
            parts = {k: v for k, v in parts.items() if k not in ("DX_SHEET", "DY_SHEET")}
        asm = cq.Assembly()
        cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "DX_SHEET": (0.9, 0.1, 0.1), "DY_SHEET": (0.1, 0.1, 0.9)}
        for nm, shp in parts.items():
            asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
        path = OUT / f"cell_r2p30_{tag}.step"
        asm.export(str(path), exportType="STEP")
        raw = path.read_text(encoding="utf-8")
        it = P.intended()
        snapped, stats = snap_step_text(raw, it, description=f"sidequest RO3006 r2.30 {tag} cell, mm, 0.001 mm grid")
        path.write_text(snapped, encoding="utf-8", newline="\n")
        bad = inexact_literals(path.read_text(encoding="utf-8"), it)
        vs = verify_snap(raw, snapped)
        if sheets:
            res, ok = P.checks(path)
        else:
            from sidequest_board_k2_lp_package import read_names
            names = read_names(path)
            res, ok = {}, True
            solids = names["METAL_PEC"].Solids()
            patch = [s for s in solids if s.BoundingBox().zmin > -1e-6]
            gnd = [s for s in solids if s.BoundingBox().zmax < -1.2]
            exp = (np.pi * R**2 - (2 * (2 * ls) * P.G - P.G**2)) * P.T
            c1 = {"ok": set(names) == {"METAL_PEC", "SUBSTRATE", "AIR"} and len(patch) == 1 and len(gnd) == 1 and abs(patch[0].Volume() - exp) < 2e-4, "names": sorted(names), "patch_volume": patch[0].Volume(), "expected": exp, "slit_half_mm": ls, "rim_mm": R - ls}
            sb = names["SUBSTRATE"].BoundingBox()
            ab = names["AIR"].BoundingBox()
            c2 = {"ok": abs(sb.zmin + P.H) < 1e-9 and abs(sb.zmax) < 1e-9 and len(names["AIR"].Faces()) == 6 and abs(ab.zmin + P.H + P.GT) < 1e-9 and abs(ab.zmax - P.LID) < 1e-9}
            res = {"neutral_cell": c1, "substrate_and_air_box": c2}
            ok = c1["ok"] and c2["ok"]
        exact = {"ok": len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001}
        res["EXACT_DIMS"] = exact
        ok = bool(ok and exact["ok"])
        allres[tag] = {"all_pass": ok, "checks": res}
        print(tag, {k: v["ok"] for k, v in res.items()}, "ALL", ok, flush=True)
        if not ok:
            sys.exit("checks failed")
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEPs", "cells": allres}, indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "geom"
    if what == "geom":
        geom()
    elif what == "pred":
        pr = predictions()
        (OUT / "predictions_raw.json").write_text(json.dumps(pr, indent=1, default=float), encoding="utf-8", newline="\n")
        for k, v in pr["cases"].items():
            print(k, {t: (None if v[t] is None else {x: round(y, 3) for x, y in v[t].items() if isinstance(y, float)}) for t in ("lid", "radiation_top")}, flush=True)
        print(json.dumps(pr["r2.20_reference_90deg_lid"], indent=0)[:900])
