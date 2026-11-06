"""A/B board (hfss_exports/sidequest_board_k2_lp: RT5880 1.575, 150 mm hexagon, 1430 cells, GA 180 + GB 180 switched cells (2 diodes each), 526 static open, 544 static strap):
4-state calibration of SWHB against the user's HFSS gains.  MODEL_ONLY.
HFSS (10.5 GHz, LP): (GA,GB) = (ON,OFF): 13.5 dBi @ (25,0), 1.2 @ (25,180); (ON,ON): 8.9 / 8.9; (OFF,ON): 1.4 @ (25,0), 13.3 @ (25,180); (OFF,OFF): 10.4 / 8.5.   SLL about -5 dB; about 12 % of the power within theta < 30, about 70 % to the horizon / back.
    python sidequest_ab_calib.py raw | grid"""
import csv
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_sll_def as SLD
import sidequest_swhb as S

RES = ROOT / "results/sidequest_hybrid"
PKG = ROOT / "hfss_exports/sidequest_board_k2_lp"
P = 6.4e-3
ENV = S.Env(2.2, 0.0009, 1.575e-3, 10.5e9)
HFSS = {("ON", "OFF"): (13.5, 1.2), ("ON", "ON"): (8.9, 8.9), ("OFF", "ON"): (1.4, 13.3), ("OFF", "OFF"): (10.4, 8.5)}  # (gain at (25,0), gain at (25,180)) [dBi]
STATES = list(HFSS)


def load():
    rows = list(csv.DictReader(open(PKG / "cells.csv")))
    xy = np.array([(float(r["x_mm"]) * 1e-3, float(r["y_mm"]) * 1e-3) for r in rows])
    kind = np.array([r["kind"] for r in rows])
    grp = np.array([r["group"] for r in rows])
    lib = json.loads((RES / "k2_default_lib.json").read_text())["rows"]

    def Y(mode):
        r = next(x for x in lib if x["family"] == "A" and abs(x["ls"] - 2.8) < 1e-9 and x["mode"] == mode)
        return complex(*r["Y_mS"][0]) * 1e-3

    return xy, kind, grp, {"open": Y("open"), "ON": Y("ON"), "strap": Y("strap")}, lib


def board_y(kind, grp, Yd, state, cal=1.0, db=0.0):
    """Per-cell admittance for (GA, GB) = state; contrast calibration about the ON/strap anchor and a common susceptance shift db [S]."""
    on = Yd["ON"]
    op = Yd["open"]
    y = np.where(kind == "static_open", op, np.where(kind == "static_strap", Yd["strap"], 0j)).astype(complex)
    for g, s in zip(("GA", "GB"), state):
        m = (kind == "switched") & (grp == g)
        y[m] = on if s == "ON" else op
    y = on + cal * (y - on) + 1j * db
    return y


def run(args):
    cal, db, state, grid = args
    xy, kind, grp, Yd, _ = load()
    y = board_y(kind, grp, Yd, state, cal, db)
    r = S.run_state(ENV, xy, y, lambda kr: 0j, (0j, 0j), P, 1.0, 6.5e-3, n_iter=6, grid=grid)
    if not r["ok"]:
        return {"ok": False, "cal": cal, "db": db, "state": state}
    th, ph, u = r["th"], r["ph"], r["u"]
    thd, phd = np.degrees(th), np.degrees(ph)
    g = lambda t, p_: 10 * math.log10(4 * math.pi * u[np.argmin(abs(thd - t)), np.argmin(abs(phd - p_))] + 1e-30)
    TH = np.meshgrid(th, ph, indexing="ij")[0]
    w = u * np.sin(TH)
    out = {"ok": True, "cal": cal, "db": db, "state": state, "beta_k0": r["beta_k0"], "alpha_loss": float(r["alpha_loss"]), "P_rad": r["P_rad_unit"], "P_edge": r["P_edge"],
           "G0": g(25, 0), "G180": g(25, 180), "frac_lt30": float(w[thd <= 30].sum() / w.sum()), "frac_gt60": float(w[thd >= 60].sum() / w.sum())}
    i = np.unravel_index(np.argmax(u), u.shape)
    out["peak_dir"] = [float(thd[i[0]]), float(phd[i[1]])]
    out["G_peak"] = 10 * math.log10(4 * math.pi * u.max())
    try:
        a = SLD.analyze(thd, phd, u)
        out["SLL_db"] = a["SLL_db"]
    except Exception:
        out["SLL_db"] = None
    return out


def report(res, eta_db=None):
    for r in res:
        if not r["ok"]:
            print(r["state"], "no root")
            continue
        h = HFSS[r["state"]]
        print("%s cal %.2f db %+.1f mS | beta %.3f Prad %.3f edge %.3f | G(25,0) %.1f (HFSS %.1f) G(25,180) %.1f (HFSS %.1f) | peak %s Gpk %.1f | SLL %s | <30 %.2f >60 %.2f" % (
            "/".join(r["state"]), r["cal"], r["db"] * 1e3, r["beta_k0"], r["P_rad"], r["P_edge"], r["G0"], h[0], r["G180"], h[1], r["peak_dir"], r["G_peak"],
            None if r["SLL_db"] is None else round(r["SLL_db"], 1), r["frac_lt30"], r["frac_gt60"]), flush=True)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "raw"
    if what == "raw":
        with ProcessPoolExecutor(4) as ex:
            res = list(ex.map(run, [(1.0, 0.0, s, (2.0, 5.0)) for s in STATES]))
        report(res)
    if what == "grid":
        combos = [(c, d * 1e-3) for c in (0.8, 0.9, 1.0, 1.1) for d in (-1.0, -0.5, 0.0, 0.5, 1.0)]
        jobs = [(c, d, s, (2.0, 5.0)) for c, d in combos for s in STATES]
        with ProcessPoolExecutor(16) as ex:
            res = list(ex.map(run, jobs))
        (RES / "ab_calib_grid.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
        for c, d in combos:
            rr = [r for r in res if r["cal"] == c and abs(r["db"] - d) < 1e-12]
            if not all(r["ok"] for r in rr):
                print("cal %.1f db %+.1f: no root in some state" % (c, d * 1e3))
                continue
            mod = np.array([[r["G0"], r["G180"]] for r in rr])
            hf = np.array([HFSS[r["state"]] for r in rr])
            eta = float(np.mean(mod - hf))
            resid = mod - eta - hf
            iso_a = (mod[0, 0] - mod[0, 1]) - (hf[0, 0] - hf[0, 1])
            iso_b = (mod[2, 1] - mod[2, 0]) - (hf[2, 1] - hf[2, 0])
            print("cal %.1f db %+.1f mS | eta %.2f dB | max|resid| %.2f | resid %s | iso err A %+.2f B %+.2f | beta %.3f" % (c, d * 1e3, eta, np.abs(resid).max(), np.round(resid, 1).tolist(), iso_a, iso_b, rr[0]["beta_k0"]), flush=True)
