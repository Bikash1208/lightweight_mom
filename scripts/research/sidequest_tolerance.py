"""Tolerance study of the RO3006 1.27 mm r 2.2 family (approved cell): disc radius +-0.05 mm, eps_r 6.15 +-0.15, diode slit gap 0.4 +-0.03 mm.  MODEL_ONLY.
Cases: centre, 6 one-at-a-time, 8 vertex corners.  Per case: pair + static library (MoM, mesh 64), direct A0 (dY x 3.13/2.63, a calibration from the baseline HFSS, NOT first-principles), beta/k0, X, then the board (design (i), <= 514 diodes) with the launch 0.71 (bare) and 0.96 (pad, unvalidated upper bound).
The disc-disc gap (6.4 - 2 r) moves with r: r +-0.05 -> gap +-0.10 mm, so a +-0.03 mm gap is covered by r +-0.015 (inside the r corners).
    python sidequest_tolerance.py mom | board  ->  results/sidequest_hybrid/k2_tolerance*.json"""
import itertools
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
F = 10.5e9
H, TD = 1.27, 0.002
FR = (0.0, 1.0 / 2.8, 1.8 / 2.8, 2.2 / 2.8, 2.6 / 2.8)
DYCAL = 3.13 / 2.63
R0, E0, G0 = 2.2, 6.15, 0.4
CASES = {"centre": (R0, E0, G0)}
for k, (dr, de, dg) in {"r-": (-0.05, 0, 0), "r+": (0.05, 0, 0), "e-": (0, -0.15, 0), "e+": (0, 0.15, 0), "g-": (0, 0, -0.03), "g+": (0, 0, 0.03)}.items():
    CASES[k] = (round(R0 + dr, 4), round(E0 + de, 4), round(G0 + dg, 4))
for sr, se, sg in itertools.product((-1, 1), repeat=3):
    CASES["c" + "".join("+" if s > 0 else "-" for s in (sr, se, sg))] = (round(R0 + 0.05 * sr, 4), round(E0 + 0.15 * se, 4), round(G0 + 0.03 * sg, 4))


def mom_job(a):
    warnings.simplefilter("ignore")
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.hologram.coded_cell_eval import sheet_tensor
    from src.materials import Substrate
    name, kind, fr = a
    r, e, g = CASES[name]
    gc.SUB = Substrate(name, e, TD, H * 1e-3)
    L.R = r
    L.G = g
    lm = round(r - 0.2, 3)
    ls = round(lm * fr, 3)
    y = sheet_tensor(L.cell(lm if kind == "ON" else ls, min(0.8, r - 0.4), kind), F, 64)
    return {"fam": name, "kind": kind, "fr": fr, "Y_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3]}


def stage_mom():
    jobs = []
    for n in CASES:
        jobs += [(n, "open", 1.0), (n, "ON", 1.0), (n, "none", 0.0)] + [(n, "open", f) for f in FR[1:]] + [(n, "strap", 2.2 / 2.8), (n, "strap", 1.0)]
    with ProcessPoolExecutor(40) as ex:
        res = list(ex.map(mom_job, jobs))
    Path("results/sidequest_hybrid/k2_tolerance_lib.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    import sidequest_a0_direct as A
    import sidequest_leaf_cell as L1
    from src.materials import Substrate
    out = {}
    for n, (r, e, g) in CASES.items():
        rows = [x for x in res if x["fam"] == n]

        def gg(k, f, rows=rows):
            return next(x for x in rows if x["kind"] == k and abs(x["fr"] - f) < 1e-9)["Y_mS"]

        yo, yn = complex(*gg("open", 1.0)), complex(*gg("ON", 1.0))
        mid, half = (yo.imag + yn.imag) / 2, (yo.imag - yn.imag) / 2
        lev = {"solid": (complex(*gg("none", 0.0)).imag - mid) / half, **{f"open_{f:.2f}": (complex(*gg("open", f)).imag - mid) / half for f in FR[1:]},
               "strap_0.79": (complex(*gg("strap", 2.2 / 2.8)).imag - mid) / half, "strap_1.0": (complex(*gg("strap", 1.0)).imag - mid) / half}
        dy = abs(yo - yn) * 1e-3
        a0m, bk = A.a0_direct(mid * 1e-3, dy, e, H * 1e-3, F, tand=TD)
        a0c, _ = A.a0_direct(mid * 1e-3, dy * DYCAL, e, H * 1e-3, F, tand=TD)
        sub = Substrate(n, e, TD, H * 1e-3)
        xo, xn = L1.xs_from_y(yo * 1e-3, sub, F, 3e8 / F / 4)["X_lid_ohm"], L1.xs_from_y(yn * 1e-3, sub, F, 3e8 / F / 4)["X_lid_ohm"]
        out[n] = {"stack": "RO3006_1.27", "r": r, "er": e, "G": g, "Ybar_mS": mid, "dY_mS": dy * 1e3, "A0_mom": float(a0m), "A0_cal": float(a0c), "beta_k0": float(bk), "X_off": xo, "X_on": xn, "dX": xo - xn,
                  "levels": lev, "static_range": [min(lev.values()), max(lev.values())], "stack_par": [e, TD, H], "Rs_on": 2.44 * 1.0 / g}
        print(n, "r %.3f er %.2f G %.2f  dY %.2f A0cal %.2f beta %.3f dX %.0f" % (r, e, g, dy * 1e3, a0c, bk, xo - xn), flush=True)
    Path("results/sidequest_hybrid/k2_tolerance_summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")


def stage_board():
    os.environ["FAM_SUMMARY"] = "results/sidequest_hybrid/k2_tolerance_summary.json"
    import sidequest_family_design as FD
    jobs = [(n, l, t, 0) for n in CASES for l, t in ((0.71, "bare"), (0.96, "pad"))]
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(FD.board_job, jobs))
    Path("results/sidequest_hybrid/k2_tolerance_boards.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        print(r["fam"], r["tag"], "L %.2f lam %.3f diodes %d gain %s conv %.3f edge %.3f" % (r["launch"], r["lam"], r["diodes"], [round(x, 2) for x in r["gain_model"]], r["conv"], r["edge_surv"]), flush=True)


if __name__ == "__main__":
    for s in sys.argv[1:]:
        {"mom": stage_mom, "board": stage_board}[s]()
