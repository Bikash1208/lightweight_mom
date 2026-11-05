"""Re-scaled cell families -> direct A0 -> boards (design (i): one group, one line, default-OFF two-beam, <= 514 diodes).  MODEL_ONLY.
Stages: (1) MoM of the pair (open / ON) and the static library per family; (2) direct A0 (first order, `sidequest_a0_direct`), dY of the MoM corrected by the HFSS/MoM ratio of the baseline cell (3.13/2.63);
(3) launcher (`sidequest_launcher`, validated on the baseline probe: 0.611 vs HFSS 0.627); (4) engine boards.
    python sidequest_family_design.py <stage> [...]   stages: mom | launch | board"""
import json
import math
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
F = 10.5e9
FAM = {"RT5880_1.575_r3.15": ("RT5880_1.575", 2.2, 0.0009, 1.575, 3.15), "RT5880_3.175_r3.15": ("RT5880_3.175", 2.2, 0.0009, 3.175, 3.15),
       "RO3006_1.27_r2.2": ("RO3006_1.27", 6.15, 0.002, 1.27, 2.2), "RO3006_1.27_r2.25": ("RO3006_1.27", 6.15, 0.002, 1.27, 2.25), "baseline_r3.0": ("RT5880_1.575", 2.2, 0.0009, 1.575, 3.0)}
FR = (0.0, 1.0 / 2.8, 1.8 / 2.8, 2.2 / 2.8, 2.6 / 2.8)  # slit half-length / (r - 0.2) of the baseline statics
DYCAL = 3.13 / 2.63


def mom_job(a):
    warnings.simplefilter("ignore")
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.hologram.coded_cell_eval import sheet_tensor
    from src.materials import Substrate
    fam, kind, fr = a
    _n, e, t, h, r = FAM[fam]
    gc.SUB = Substrate(fam, e, t, h * 1e-3)
    L.R = r
    lm = round(r - 0.2, 3)
    ls = round(lm * fr, 3)
    y = sheet_tensor(L.cell(lm if kind == "ON" else ls, min(0.8, r - 0.4), kind), F, 64)
    return {"fam": fam, "kind": kind, "fr": fr, "ls": ls, "Y_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3]}


def stage_mom():
    jobs = []
    for fam in FAM:
        jobs += [(fam, "open", 1.0), (fam, "ON", 1.0), (fam, "none", 0.0)] + [(fam, "open", f) for f in FR[1:]] + [(fam, "strap", 2.2 / 2.8), (fam, "strap", 1.0)]
    with ProcessPoolExecutor(40) as ex:
        res = list(ex.map(mom_job, jobs))
    Path("results/sidequest_hybrid/k2_family_lib.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    print("mom", len(res), flush=True)


def fam_summary():
    import sidequest_a0_direct as A
    import sidequest_leaf_cell as L1
    from src.materials import Substrate
    res = json.loads(Path("results/sidequest_hybrid/k2_family_lib.json").read_text())
    out = {}
    for fam, (stk, e, t, h, r) in FAM.items():
        rows = [x for x in res if x["fam"] == fam]

        def g(k, f, rows=rows):
            return next(x for x in rows if x["kind"] == k and abs(x["fr"] - f) < 1e-9)["Y_mS"]

        yo, yn = complex(*g("open", 1.0)), complex(*g("ON", 1.0))
        mid, half = (yo.imag + yn.imag) / 2, (yo.imag - yn.imag) / 2
        lev = {"solid": (complex(*g("none", 0.0)).imag - mid) / half, **{f"open_{f:.2f}": (complex(*g("open", f)).imag - mid) / half for f in FR[1:]},
               "strap_0.79": (complex(*g("strap", 2.2 / 2.8)).imag - mid) / half, "strap_1.0": (complex(*g("strap", 1.0)).imag - mid) / half}
        dy = abs(yo - yn) * 1e-3
        a0m, bk = A.a0_direct(mid * 1e-3, dy, e, h * 1e-3, F, tand=t)
        a0c, _ = A.a0_direct(mid * 1e-3, dy * DYCAL, e, h * 1e-3, F, tand=t)
        sub = Substrate(stk, e, t, h * 1e-3)
        H = 3e8 / F / 4
        xo, xn = L1.xs_from_y(yo * 1e-3, sub, F, H)["X_lid_ohm"], L1.xs_from_y(yn * 1e-3, sub, F, H)["X_lid_ohm"]
        out[fam] = {"stack": stk, "r": r, "Ybar_mS": mid, "dY_mS": dy * 1e3, "A0_mom": float(a0m), "A0_cal": float(a0c), "beta_k0": float(bk), "X_off": xo, "X_on": xn, "dX": xo - xn,
                    "levels": lev, "static_range": [min(lev.values()), max(lev.values())], "stack_par": [e, t, h]}
    return out


def launch_job(a):
    warnings.simplefilter("ignore")
    import sidequest_launcher as S
    name, e, tdl, h, t, ys = a
    k0 = 2 * math.pi * F / S.C0
    t_ = t * 1e-3
    i0 = math.sin(k0 * t_)
    r = S.run(F, e, tdl, h * 1e-3, t_, lambda z: math.sin(k0 * (t_ - z)) if z > 0 else i0, i0, ys=ys * 1e-3)
    return {"name": name, "t_mm": t, "ys_mS": ys, **{k: float(v) for k, v in r.items()}}


def stage_launch():
    fs = fam_summary()
    Path("results/sidequest_hybrid/k2_family_summary.json").write_text(json.dumps(fs, indent=1), encoding="utf-8", newline="\n")
    jobs = []
    for fam in FAM:
        e, tdl, h = fs[fam]["stack_par"]
        for ys in (0.0, fs[fam]["Ybar_mS"]):
            for t in (4.0, 6.0):
                jobs.append((fam, e, tdl, h, t, ys))
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(launch_job, jobs))
    Path("results/sidequest_hybrid/k2_launcher.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    for r in res:
        print(r["name"], "t", r["t_mm"], "Ys", round(r["ys_mS"], 2), "launch %.3f R_in %.1f beta/k0 %.3f" % (r["launch"], r["R_in_ohm"], r["beta_over_k0"]), flush=True)


def board_job(a):
    warnings.simplefilter("ignore")
    fam, launch, tag, rim_db = a
    fs = json.loads(Path("results/sidequest_hybrid/k2_family_summary.json").read_text())[fam]
    os.environ["K4_A0"] = "%.4f" % fs["A0_cal"]
    os.environ["K4_BETA_K0"] = "%.4f" % fs["beta_k0"]
    import numpy as np
    import sidequest_coded_k4 as K
    import sidequest_k2_default_design as D
    K.A0 = fs["A0_cal"]
    K.BETA_K0 = fs["beta_k0"]
    K.LAUNCH = launch
    K.AA_BG = 0.15 * (fs["stack_par"][1] / 0.0009) ** 0.5 if fs["stack_par"][1] > 0.0009 else 0.15  # background loss scaled with sqrt(tan d): flagged
    levels = sorted({round(float(v), 4) for v in fs["levels"].values()} | {1.0})
    D.setup(F, 150e-3, D.MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    opt = D.options(levels, [(1.0, -1.0, 1.0)], [(0, 1)])
    rim = np.zeros(bd.n)
    if rim_db:
        a_ = math.log(10 ** (rim_db / 10)) / (2 * 12.8e-3)
        rim = np.where(bd.rho > bd.redge - 12.8e-3, a_, 0.0)

    def ap(bd_, m, loss):
        alpha = K.A0 * np.mean(K.vec2(bd_, m), axis=0)
        aa = K.AA_BG + 2 * K.AA0 * np.mean(loss, axis=0) + rim
        integ, ie = bd_.profile(alpha, aa)
        amp = np.sqrt(2 * K.A0 / bd_.rho) * np.exp(-integ)
        w = 2 * alpha * np.exp(-2 * integ) / bd_.rho * D.st.PITCH ** 2 / (2 * np.pi)
        return amp, float(w.sum()), float(np.mean(np.exp(-2 * ie)))

    K.amp_profile = ap
    eps = [0.0, np.radians(180)]

    def run(lam):
        return D.run_design(bd, opt, eps, lam, 2.0, mu=30.0)

    lo, hi, r = 0.02, 1.5, None
    for _ in range(7):
        mid = 0.5 * (lo + hi)
        rr = run(mid)
        if int(rr["nd"].sum()) > 514:
            lo = mid
        else:
            hi = mid
            r = rr
    if r is None:
        r = run(hi)
    iso = D.isolation(r["states"])
    return {"fam": fam, "tag": tag, "launch": launch, "rim_db": rim_db, "lam": float(r["lam"]), "diodes": int(r["nd"].sum()), "gain_model": [s["gain_dbi"] for s in r["states"]],
            "sll_model": [s["sll_db"] for s in r["states"]], "conv": r["conv"], "edge_surv": r["edge"], "iso": iso}


def stage_board():
    la = json.loads(Path("results/sidequest_hybrid/k2_launcher.json").read_text())
    cal = 0.627 / next(x for x in la if x["name"] == "baseline_r3.0" and x["ys_mS"] == 0.0 and x["t_mm"] == 6.0)["launch"]
    jobs = [("baseline_r3.0", 0.627, "baseline", 0)]
    for fam in ("RT5880_1.575_r3.15", "RT5880_3.175_r3.15", "RO3006_1.27_r2.2", "RO3006_1.27_r2.25"):
        bare = next(x for x in la if x["name"] == fam and x["ys_mS"] == 0.0 and x["t_mm"] == 6.0)["launch"] * cal
        pad = next(x for x in la if x["name"] == fam and x["ys_mS"] > 0 and x["t_mm"] == 6.0)["launch"] * cal
        jobs += [(fam, 0.627, "launch0.627", 0), (fam, min(bare, 0.99), "launch_bare", 0), (fam, min(pad, 0.99), "launch_pad", 0), (fam, min(pad, 0.99), "launch_pad+rim10", 10)]
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(board_job, jobs))
    Path("results/sidequest_hybrid/k2_family_boards.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        print(r["fam"], r["tag"], "L %.3f lam %.3f diodes %d gain %s conv %.3f edge %.3f" % (r["launch"], r["lam"], r["diodes"], [round(x, 2) for x in r["gain_model"]], r["conv"], r["edge_surv"]), flush=True)


if __name__ == "__main__":
    for st in sys.argv[1:]:
        {"mom": stage_mom, "launch": stage_launch, "board": stage_board}[st]()
