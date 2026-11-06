"""Predictions (eigenmode with the 30 mm lid), HFSS_SETUP.md, states.csv, results_template.csv, sync for the via-cell package.  MODEL_ONLY."""
import csv
import json
import math
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import brentq

import sidequest_cellB_model as M
import sidequest_via_cell as V
import sidequest_via_hw as W
import sidequest_via_cell_package as PK

OUT = ROOT / "hfss_exports/sidequest_cell_via_RT5880"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\cell_via_RT5880")
RES = ROOT / "results/sidequest_hybrid"
RB = dict(W.RB, R_b=4700.0, C_p=8e-15, C_A=10e-15, cathode="via")
C_S = 0.216e-12
PPITCH = 6.4e-3
LID = 0.030


def patches():
    rows = json.loads((RES / "cellB_rt5880_a.json").read_text())
    out = {}
    for psi in (0.0, 45.0):
        out[psi] = M.Patch([r for r in rows if r["lam"] == "RT5880" and r["mesh"] == 64 and r["psi"] == psi and abs(r["gap"] - 0.4) < 1e-6])
    return out


def y_air_lid(kr, f):
    w = 2 * math.pi * f
    k0 = w / M.C0
    kz = k0 * np.sqrt(1 - kr * kr + 0j)
    if kz.imag > 0:
        kz = -kz
    return -1j * w * M.EPS0 / (kz * np.tan(kz * LID))


def beta_lid(P, state, f, **kw):
    sc = f / V.F0
    yb = W.branch2(state, f, C_S, RB, **kw)
    k0 = 2 * math.pi * f / M.C0
    F = lambda kr: y_air_lid(kr, f) + V.y_slab(kr, f) + sc * P.y_sw(min(kr, 1.5)) + yb
    ks = np.linspace(1.0005, 1.44, 700)
    v = np.array([F(k).imag for k in ks])
    ix = np.where(np.diff(np.sign(v)) != 0)[0]
    if len(ix) == 0:
        return None
    kr = brentq(lambda k: F(k).imag, ks[ix[0]], ks[ix[0] + 1], xtol=1e-10)
    return kr * k0


def eig_f(P, state, phase_deg, **kw):
    target = math.radians(phase_deg) / PPITCH
    g = lambda f: (lambda b: -1e9 if b is None else b - target)(beta_lid(P, state, f, **kw))
    fs = np.linspace(6.0e9, 13.5e9, 76)
    v = np.array([g(f) for f in fs])
    ix = np.where(np.diff(np.sign(v)) != 0)[0]
    return None if len(ix) == 0 else brentq(g, fs[ix[0]], fs[ix[0] + 1], xtol=1e3) / 1e9


def main():
    Ps = patches()
    z = PK.zsheet()
    pred = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "lumped branch (cascade) + MoM patch grid; lid 30 mm"], "cell": {"C_s_pF": C_S * 1e12, "R_b_ohm": 4700, "island_fF": 10, "resistor_fF": 8, "diode": "MA4AGP907 measured (2.44 ohm / 28.16 fF)"}}
    rows = []
    for name, psi, phases in (("x", 0.0, (90.0, 110.0)), ("45deg", 45.0, (99.0, 120.0))):  # phase = beta*p along the propagation direction; 45 deg: Bloch phases (phi cos45, phi sin45)
        for st in ("OFF", "ON"):
            for ph in phases:
                f = eig_f(Ps[psi], st, ph)
                rows.append(dict(direction=name, state=st, beta_p_deg=ph, bloch_phases_deg=(ph, 0.0) if name == "x" else (round(ph * math.cos(math.radians(45)), 2), round(ph * math.sin(math.radians(45)), 2)), f_GHz=None if f is None else round(f, 4)))
    pred["eigenmode_with_lid_GHz"] = rows
    # sensitivity at 90 deg x
    sens = {}
    for nm, kw in (("via L -8%", dict(dL=0.92)), ("via L +8%", dict(dL=1.08)), ("diode C +20%", dict(dC=1.2)), ("C_s +20%", dict(dCs=1.2)), ("kappa 0.7", dict(kappa=0.7))):
        sens[nm] = {st: round(eig_f(Ps[0.0], st, 90.0, **kw), 4) for st in ("OFF", "ON")}
    pred["sensitivity_90deg_x_GHz"] = sens
    # Q estimate from the loss: Q = omega / (2 alpha v_g)
    qs = {}
    for st in ("OFF", "ON"):
        f1, f2 = eig_f(Ps[0.0], st, 90.0), eig_f(Ps[0.0], st, 110.0)
        vg = 2 * math.pi * (f2 - f1) * 1e9 / (math.radians(20.0) / PPITCH)
        yb = W.branch2(st, V.F0, C_S, RB)
        r = V.root(Ps[0.0], yb)
        a = r[1]
        qs[st] = {"alpha_loss_Np_per_m": float(a), "vg_over_c": vg / M.C0, "Q_lossy_estimate": float(2 * math.pi * f1 * 1e9 / (2 * a * vg))}
    pred["Q_from_loss_estimate_x90"] = qs
    pred["bands"] = "f +-0.15 GHz (+-0.3 GHz for the 45 deg runs and for the first run of the hardware: the hardware inductances and pad capacitances are lumped); Q within a factor 2 (copper loss and the anti-pad leakage are not in it)"
    pred["falsifiers"] = ["ON above OFF by less than 0.15 GHz or ON below OFF (the cell does not switch as modelled)", "the shift ON - OFF at 90 deg differs from the predicted by more than 0.3 GHz (kappa / hardware parasitics wrong: re-tune C_s with the rule below)",
                          "x and y runs differ by more than 0.2 GHz without a hardware explanation", "a second mode between OFF and ON frequencies near 10.5 GHz (a resonance of the choke / pads)"]
    pred["retune_rule"] = "if the HFSS beta_mean (average of the OFF and ON beta at the same frequency, equivalently the mean sheet) is off by d, change C_s by about d / 0.0027 per 5 % (C_s +5 % = beta_mean +0.0024); the step moves -1.5 % per +5 % C_s"
    pred["sheet_impedances_ohm_per_sq_at_10p5GHz"] = {k: [round(v.real, 3), round(v.imag, 3)] for k, v in z.items()}
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    with open(OUT / "states.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["state", "CS_SHEET (R, X ohm/sq)", "D_SHEET (R, X ohm/sq)", "RB_SHEET (R, X ohm/sq)"])
        w.writerow(["OFF (default)", f"{z['CS'].real:.3f}, {z['CS'].imag:.3f}", f"{z['D_OFF'].real:.3f}, {z['D_OFF'].imag:.3f}", f"{z['RB'].real:.3f}, {z['RB'].imag:.3f}"])
        w.writerow(["ON", f"{z['CS'].real:.3f}, {z['CS'].imag:.3f}", f"{z['D_ON'].real:.3f}, {z['D_ON'].imag:.3f}", f"{z['RB'].real:.3f}, {z['RB'].imag:.3f}"])
    with open(OUT / "results_template.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["direction", "state", "phase_x_deg", "phase_y_deg", "boundary_top", "f_GHz_real", "Q", "mode_id_note"])
        for d, ph in (("x", (90, 0)), ("x", (110, 0)), ("y", (0, 90)), ("y", (0, 110)), ("45deg", (70, 70)), ("45deg", (85, 85))):
            for st in ("OFF", "ON"):
                w.writerow([d, st, ph[0], ph[1], "PEC lid 30 mm", "", "", ""])
        for st in ("OFF", "ON"):
            w.writerow(["x", st, 90, 0, "Radiation (open top)", "", "", ""])
    g = json.loads((OUT / "geometry_checks.json").read_text())
    f90 = {(r["direction"], r["state"], r["beta_p_deg"]): r["f_GHz"] for r in rows}
    setup = f"""# HFSS_SETUP - via-switched GB cell WITH the real bias hardware, RT5880 1.575 mm, eigenmode OFF / ON (x, y, 45 deg) + radiation-top run. MODEL_ONLY. One page.
**1. Materials**: `METAL_PEC` (FOUR DC-isolated lumps: patch+via+land | ground+K pad+blind via | island+A_up | bias pad+trunk; ONE named part) -> `pec`. `SUBSTRATE`: RT5880 eps_r 2.2, tan d 0.0009, z -1.575..0. `PLY`: bond ply eps_r 3.52, tan d 0.004 (RO4450F-class, Rogers datasheet values via mirrors; compatibility with PTFE to be confirmed), z -1.711..-1.610. `AIR`: vacuum, 6.4 x 6.4, z -4.746 (3 mm under the back layer) .. +30, ONE face per side wall, FIRST in the tree (overlaps the solids on purpose).
**Files**: `cell_via_RT5880.step`, `states.csv`, `predictions.json`, `results_template.csv`, `geometry_checks.json` (all PASS on the re-imported STEP, EXACT_DIMS PASS: {g['checks']['EXACT_DIMS']['n_snapped']} literals snapped, max {g['checks']['EXACT_DIMS']['max_snap_mm']} mm).
**Geometry** (mm): patch 6.0 x 6.0 (z 0..0.035), via r 0.15 through the substrate, the anti-pad slot in L2 (x -0.55..1.15, y -0.45..1.25), the ply and the back layer L3 (z -1.746..-1.711): land pad 0.5 x 0.5 at the via, series-C gap `CS_SHEET` (x 0.25..0.55), island A (x 0.55..0.95, y -0.2..0.2) with the pad `A_up` (y 0.75..1.15) inside the slot, diode gap `D_SHEET` (x 0.95..1.25), cathode pad K (x 1.25..1.85, y -0.3..0.3) with a blind via r 0.1 at (1.55, 0) to L2, ballast gap `RB_SHEET` (y 1.15..1.65, x 0.65..0.85), bias pad (y 1.65..2.05), connection (y 1.85) to the trunk at x 2.5 that runs the full pitch (continues into the neighbour cell). Sheets are at z -1.728 (inside the L3 copper); `CS_SHEET` and `D_SHEET` are 0.3 x 0.3 (W/L 1), `RB_SHEET` is 0.2 x 0.5 (W/L 0.4).
**2. Boundaries**: Primary/Secondary (Lattice Pairs) on the four side faces of `AIR` (they also cover the substrate, ply, ground and trunk); Perfect E on the BOTTOM face (z -4.746; the back cavity couples through the anti-pad only at the 1e-4 level); top face: **Perfect E (lid 30 mm) for the eigenmode runs, Radiation (open top) for the two radiation-top runs**.
**3. Sheets as Impedance boundaries (ohm/sq, set at 10.5 GHz; they are frequency-independent in HFSS, so an eigenfrequency 0.5 GHz away carries a small error; use Lumped RLC instead if your eigenmode solver supports it: CS C = 0.216 pF, D R = 2.44 ohm (ON) / open (OFF) with C = 28.16 fF, RB R = 4.7 kohm with C = 8 fF)**: `CS_SHEET` R {z['CS'].real:.2f}, X {z['CS'].imag:.2f} (0.216 pF + 0.25 nH ESL + 0.2 ohm ESR); `RB_SHEET` R {z['RB'].real:.1f}, X {z['RB'].imag:.1f} (4.7 kohm || 8 fF, scaled by W/L 0.4); `D_SHEET` **OFF: R 0, X {z['D_OFF'].imag:.1f}** (28.16 fF) and **ON: R {z['D_ON'].real:.2f}, X 0**. One variable switches only `D_SHEET`.
**4. Eigenmode**: x runs (phase x, y) = (90, 0) and (110, 0) deg; y runs (0, 90) and (0, 110); 45 deg runs (70, 70) and (85, 85) deg (beta p = 99 and 120 deg along the diagonal); OFF and ON each: 12 runs. Min frequency 5 GHz, 8 modes, delta f 0.1 %; mesh: substrate 0.8 mm, ply 0.05 mm, patch 0.3 mm, sheets 0.05 mm, vias 0.05 mm; Clean Up Solutions before each run. The bound mode is the lowest TM0-like mode below the light line; identify it by the Ez fraction. **The cell is NOT x/y symmetric (the hardware lies along x and y): report x and y separately.**
**5. Radiation-top run**: top = Radiation, x run at 90 deg, OFF and ON: complex f -> Q; alpha = omega / (2 Q v_g) with v_g from the 90 / 110 deg frequencies. Record Q and the electric-field energy fraction in the back cavity (leakage through the anti-pad).
**Prediction (MODEL_ONLY, bands in `predictions.json`)**: x, 90 deg, lid: OFF {f90[('x','OFF',90.0)]} GHz, ON {f90[('x','ON',90.0)]} GHz; x, 110 deg: OFF {f90[('x','OFF',110.0)]}, ON {f90[('x','ON',110.0)]}; 45 deg (99 deg): OFF {f90[('45deg','OFF',99.0)]}, ON {f90[('45deg','ON',99.0)]}. Cell numbers: step 17.9 mS, beta/k0 mean 1.062 / OFF 1.19 / ON 1.015, ballast and parasitics included. **Pass**: f within +-0.05 GHz of the prediction band centre after the C_s re-tune rule in `predictions.json` and alpha within +-25 %.
**Falsifiers** (`predictions.json`): ON not above OFF by about the predicted shift; a second mode near 10.5 GHz (choke / pad resonance); x and y differing by more than 0.2 GHz without a hardware cause.
"""
    (OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
    SYNC.mkdir(parents=True, exist_ok=True)
    for f in OUT.iterdir():
        shutil.copy2(f, SYNC / f.name)
    print(json.dumps(rows, default=float)[:1500])
    print(pred["sensitivity_90deg_x_GHz"], pred["Q_from_loss_estimate_x90"])


if __name__ == "__main__":
    main()
