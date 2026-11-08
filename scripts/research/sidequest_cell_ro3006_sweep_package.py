"""HFSS cell package, RO3006 1.27 mm radius sweep: four STEPs r 2.15 / 2.20 / 2.25 / 2.30 (rim 0.2, cross slit 0.4, lattice 6.4, diode sheets 0.4 x 1.0 at 0.8 mm), same setup as `hfss_exports/sidequest_cell_RO3006_r2p2`.
Eigenmode at 110 deg along x (the board's per-cell phase, beta/k0 about 1.3) for all four radii, plus 90 deg for r 2.20 only (ties to the HFSS cell check already run: ON 10.3617, OFF 10.0193 GHz).
Predictions MODEL_ONLY: loaded Bloch tensor, MoM mesh 64 (0.1 mm pixel: radius steps are below the pixel, see K2_POINTING.md), transverse resonance with the 30 mm PEC lid; plus an 'HFSS-offset' column = raw + the r 2.20 / 90 deg HFSS offsets (assumption).
    python sidequest_cell_ro3006_sweep_package.py [geom|pred|all]  ->  hfss_exports/sidequest_cell_RO3006_sweep/  (+ sync C:\\hmsa_sync\\data\\sidequest\\cell_RO3006_sweep)"""
import csv
import json
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

OUT = ROOT / "hfss_exports/sidequest_cell_RO3006_sweep"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\cell_RO3006_sweep")
RADII = (2.15, 2.20, 2.25, 2.30)
CASES = [(r, 110.0) for r in RADII] + [(2.20, 90.0)]
HFSS_90_R220 = {"OFF": 10.0193, "ON": 10.3617}
RAW_90_R220 = {"OFF": 9.8989, "ON": 10.3558}       # package prediction for r 2.20 at 90 deg (the run HFSS checked)
OFFSET = {k: HFSS_90_R220[k] - RAW_90_R220[k] for k in HFSS_90_R220}


def tag(r):
    return "r%s" % ("%.2f" % r).replace(".", "p")


def set_r(r):
    P.R = r
    P.LS = round(r - 0.2, 6)


def geom_all():
    OUT.mkdir(parents=True, exist_ok=True)
    allres = {}
    for r in RADII:
        set_r(r)
        parts = P.build()
        asm = cq.Assembly()
        cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "DX_SHEET": (0.9, 0.1, 0.1), "DY_SHEET": (0.1, 0.1, 0.9)}
        for nm, shp in parts.items():
            asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
        path = OUT / f"cell_RO3006_{tag(r)}.step"
        asm.export(str(path), exportType="STEP")
        raw = path.read_text(encoding="utf-8")
        it = P.intended()
        snapped, stats = snap_step_text(raw, it, description=f"sidequest RO3006 r{r} diode cell, mm, 0.001 mm grid")
        path.write_text(snapped, encoding="utf-8", newline="\n")
        bad = inexact_literals(path.read_text(encoding="utf-8"), it)
        vs = verify_snap(raw, snapped)
        res, ok = P.checks(path)
        exact = {"ok": len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001}
        res["EXACT_DIMS"] = exact
        ok &= exact["ok"]
        allres[tag(r)] = {"r_mm": r, "all_pass": bool(ok), "checks": res}
        print(tag(r), {k: v["ok"] for k, v in res.items()}, "ALL", ok, flush=True)
        if not ok:
            sys.exit("checks failed for r %s" % r)
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEPs", "cases": allres}, indent=1, default=float), encoding="utf-8", newline="\n")


def pjob(a):
    import warnings as w
    w.simplefilter("ignore")
    import sidequest_cell_search as C
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.materials import Substrate
    r, ph, mode, f = a
    gc.SUB = Substrate("ro3006", P.ER, P.TAND, P.H * 1e-3)
    L.R = r
    beta = np.radians(ph) / (P.P * 1e-3)
    return C.bloch_loaded(L.cell(round(r - 0.2, 3), P.POS, mode), f, 64, beta, 0.0)


def fgrid(ph):
    return np.arange(8.4, 13.01, 0.4) * 1e9 if ph > 100 else np.arange(8.4, 11.61, 0.4) * 1e9


def predictions():
    jobs = [(r, ph, m, f) for r, ph in CASES for m in ("open", "ON") for f in fgrid(ph)]
    with ProcessPoolExecutor(40) as ex:
        res = list(ex.map(pjob, jobs))
    it = iter(res)
    rows = []
    for r, ph in CASES:
        beta = np.radians(ph) / (P.P * 1e-3)
        fll = 299792458.0 * beta / (2 * np.pi) / 1e9
        fg = fgrid(ph)
        fr = {}
        for mode in ("open", "ON"):
            yt = np.array([next(it) for _ in fg])
            cs = CubicSpline(fg, yt[:, 0, 0].imag)
            fine = np.linspace(fg[0], fg[-1], 2000)
            v = np.array([P.stack_b(beta, f, P.ER, P.H) + float(cs(f)) for f in fine])
            ix = np.where(np.diff(np.sign(v)) != 0)[0]
            roots = [brentq(lambda f, cs=cs: P.stack_b(beta, f, P.ER, P.H) + float(cs(f)), fine[i], fine[i + 1], xtol=1e3) / 1e9 for i in ix]
            fr[mode] = [x for x in roots if x < fll]
        off = fr["open"][0] if fr["open"] else None
        on = fr["ON"][0] if fr["ON"] else None
        row = {"r_mm": r, "phase_deg": ph, "beta_over_k0_at_10.5GHz": float(beta / (2 * np.pi * 10.5e9 / 299792458.0)), "light_line_GHz": round(float(fll), 3), "f_OFF_GHz_raw": None if off is None else round(off, 4), "f_ON_GHz_raw": None if on is None else round(on, 4),
               "all_roots_OFF_GHz": [round(x, 4) for x in fr["open"]], "all_roots_ON_GHz": [round(x, 4) for x in fr["ON"]]}
        if off and on:
            row["shift_raw_GHz"] = round(on - off, 4)
            row["f_OFF_GHz_hfss_offset_est"] = round(off + OFFSET["OFF"], 4)
            row["f_ON_GHz_hfss_offset_est"] = round(on + OFFSET["ON"], 4)
            row["shift_hfss_offset_est_GHz"] = round(on + OFFSET["ON"] - off - OFFSET["OFF"], 4)
        rows.append(row)
        print(row, flush=True)
    return {"label": "MODEL_ONLY: loaded Bloch tensor (MoM mesh 64), transverse resonance with the 30 mm PEC lid; OFF = open slit, ON = 2.44 ohm only (Rs 6.1 ohm/sq on 0.4 x 1.0 sheets)",
            "hfss_offset_assumption": {"OFF_GHz": round(OFFSET["OFF"], 4), "ON_GHz": round(OFFSET["ON"], 4), "from": "HFSS r 2.20 at 90 deg: OFF 10.0193 / ON 10.3617 against the MoM 9.8989 / 10.3558; applied to every case (assumption: the OFF error is the MoM open-slit model; its radius / phase dependence is unknown)"},
            "pixel_caveat": "MoM mesh 64 has a 0.1 mm pixel: radius steps of 0.05 mm are staircased, a 0.05 mm step is a half pixel; the predicted slope d f / d r is the quantity this HFSS sweep measures", "rows": rows}


def docs(pred):
    def line(rw):
        return "| %.2f | %d | %.2f | %s / %s | %s / %s | %s / %s |" % (rw["r_mm"], rw["phase_deg"], rw["beta_over_k0_at_10.5GHz"], rw["f_OFF_GHz_raw"], rw["f_ON_GHz_raw"], rw.get("f_OFF_GHz_hfss_offset_est"), rw.get("f_ON_GHz_hfss_offset_est"), rw.get("shift_raw_GHz"), rw.get("shift_hfss_offset_est_GHz"))
    table = "\n".join(line(r) for r in pred["rows"])
    setup = f"""# HFSS: RO3006 1.27 mm diode cell, radius sweep r 2.15 / 2.20 / 2.25 / 2.30, eigenmode OFF/ON. EXPORT STANDARD. Predictions MODEL_ONLY (`predictions.json`).

**Same setup as `cell_RO3006_r2p2` (already run at r 2.20, 90 deg).** Four STEPs, identical except the disc radius: `cell_RO3006_r2p15.step`, `r2p20`, `r2p25`, `r2p30` (each `METAL_PEC`, `SUBSTRATE`, `AIR`, `DX_SHEET`, `DY_SHEET`); `geometry_checks.json` (all PASS on the re-imported STEPs, EXACT_DIMS PASS, 0.001 mm grid), `states.csv`, `predictions.json`, `results_template.csv`.
**1. Assign material `pec` to `METAL_PEC`** (patch + ground = one named part, 2 DC-isolated lumps; the STEP imports as vacuum). `SUBSTRATE`: RO3006 eps_r {P.ER}, tan d {P.TAND}, z -{P.H}..0. (Your eps_r 6.00 / 6.30 runs at r 2.22 = change the material only; r 2.22 itself = edit `r` in a copy: patch r 2.22, slit half-length 2.02, rim 0.2.)
**2. `AIR`** = ONE vacuum box z {-P.H-P.GT:.3f} ... +{P.LID:g} mm (ground bottom to the lid), {P.P} x {P.P}, one face per side wall, overlapping the solids on purpose, FIRST in the tree. **Primary/Secondary (Lattice Pair) on the four side faces of `AIR`** (they cover the substrate). **3. PEC lid:** Perfect E on the top face of `AIR` (z = +{P.LID:g}).
**Geometry (mm):** disc r (see file name; 35 um copper, z 0..0.035), cross slit {P.G} wide, half-length r - 0.2 (rim 0.2); ground 35 um (z {-P.H-P.GT:.3f}..{-P.H}); `DX_SHEET` {P.G} x {P.WD} at (0, +{P.POS}), `DY_SHEET` {P.WD} x {P.G} at (+{P.POS}, 0), z = {P.T/2}; sheets 0.14 mm apart.
**4. States:** OFF = **no boundary on the sheets**. ON = **ONE Impedance boundary on both sheets: Rs = {P.RS:g} ohm/sq = 2.44 ohm x W/L, W/L = {P.WD}/{P.G} = {P.WD/P.G:g}, reactance 0.**
**5. Eigenmode:** Lattice Pair 1 (x faces, Phase Delay `$ph_x`), Pair 2 (y faces, `$ph_y`). **Runs: `$ph_x` = 110deg, `$ph_y` = 0deg for all four radii; additionally `$ph_x` = 90deg for r 2.20 only (ties to the run already done: ON 10.3617, OFF 10.0193 GHz).** Min frequency 5 GHz, 8 modes, delta f 0.1 %; mesh: substrate 0.8 mm, patch 0.2 mm, sheets 0.1 mm. Clean Up Solutions before each run. The bound mode is the lowest TM0-like mode below the light line (f = c beta/2pi: 14.3 GHz at 110 deg, 11.7 GHz at 90 deg); identify it by the Ez fraction.
**Prediction (GHz, OFF / ON; MoM raw, then raw + the HFSS offsets found at r 2.20 / 90 deg [OFF {OFFSET['OFF']:+.3f}, ON {OFFSET['ON']:+.3f}, an assumption for the other cases]; shift ON - OFF raw / offset-corrected):**
| r [mm] | phase [deg] | beta/k0 at 10.5 GHz | raw OFF / ON | offset-corrected OFF / ON | shift raw / corrected |
|---|---|---|---|---|---|
{table}
**Use:** the quantity wanted is the slope of f_OFF and f_ON against r at 110 deg (it fixes dbeta/dr and A0(r), which the MoM cannot resolve below the 0.1 mm pixel) and the OFF-ON shift (the leakage swing). Falsifier: a shift below 0.1 GHz, or ON below OFF.
"""
    (OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
    with open(OUT / "states.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["state", "DX_SHEET", "DY_SHEET", "note"])
        w.writerow(["OFF", "no boundary", "no boundary", "diode OFF = open"])
        w.writerow(["ON", f"Impedance boundary Rs = {P.RS:g} ohm/sq, X = 0", "same boundary (ONE boundary on both sheets)", f"Rs = R_on * W / L = 2.44 * {P.WD} / {P.G}; W/L = {P.WD/P.G:g}"])
    with open(OUT / "results_template.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("r_mm,state,phase_deg,run,f_GHz,Q,Ez_fraction,notes\n" + "".join(f"{r},{s},110,x,,,,\n" for r in RADII for s in ("OFF", "ON")) + "2.20,OFF,90,x,,,,\n2.20,ON,90,x,,,,\n")
    SYNC.mkdir(parents=True, exist_ok=True)
    for f in OUT.iterdir():
        shutil.copy2(f, SYNC / f.name)
    print("synced", sorted(p.name for p in SYNC.iterdir()), flush=True)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("geom", "all"):
        geom_all()
    if what in ("pred", "all"):
        pr = predictions()
        (OUT / "predictions.json").write_text(json.dumps(pr, indent=1, default=float), encoding="utf-8", newline="\n")
        docs(pr)
