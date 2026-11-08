"""SUB6 slot-disc HFSS package, part 2: predictions.json, states.csv, results_template.csv, HFSS_SETUP.md (reads hfss_pred/pred_table.json).  MODEL_ONLY.
    python sub6_slotdisc_hfss_docs.py
"""
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import brentq

import sub6_slotdisc_hfss_package as pk

P, H, LID, C0, EPS0, F0, CELLS, OUT, RES = pk.P, pk.H, pk.LID, pk.C0, pk.EPS0, pk.F0, pk.CELLS, pk.OUT, pk.RES
PHASES_SCAN = range(60, 111, 5)


def bare_f(phase_deg, lid_mm):
    """Bare grounded RT5880 slab (no pattern): transverse resonance at beta = phase/p; open air (lid_mm None) or PEC lid lid_mm above the substrate. GHz."""
    from src.dispersion import sheet_model as sm

    beta = math.radians(phase_deg) / (P * 1e-3)
    eps_c = pk.EPS_R * (1 - 1j * pk.TAND)

    def b(f):
        om, k0 = 2 * math.pi * f, 2 * math.pi * f / C0
        if lid_mm is None:
            ya = sm._y_air(beta, om, k0)
        else:
            kz = np.sqrt(k0**2 - beta**2 + 0j)
            if kz.imag > 0:
                kz = -kz
            ya = -1j * om * EPS0 / (kz * np.tan(kz * lid_mm * 1e-3))
        return (ya + sm._y_slab(beta, om, eps_c, k0, H * 1e-3)).imag

    flight = beta * C0 / (2 * math.pi)  # light-line frequency; the bound mode lies below it
    fs = np.linspace(0.45 * flight, 0.9995 * flight, 3000)
    v = np.array([b(f) for f in fs])
    roots = []
    for i in range(len(fs) - 1):
        if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i] * v[i + 1] < 0:
            fr = brentq(b, fs[i], fs[i + 1])
            if abs(b(fr)) < 1e-3 * min(abs(v[i]), abs(v[i + 1])):  # a true zero, not a cot pole (poles flip sign with |b| -> inf)
                roots.append(fr / 1e9)
    return roots[-1] if roots else None  # TM0 = the highest-frequency root below the light line (no cutoff; lowest beta/k0 branch)


def main():
    tab = json.loads((RES / "hfss_pred" / "pred_table.json").read_text())
    rows = [json.loads(line) for line in (RES / "rows_RT5880.jsonl").read_text().splitlines() if line.strip()]
    mesh = {m["id"]: m for m in json.loads((RES / "mesh_check_RT5880.json").read_text())}
    eq = {"OFF": pk.parallel_equiv(pk.z_off(F0), F0), "ON": pk.parallel_equiv(pk.z_on(F0), F0)}
    eq_f = {f"{f:.1f}": {"OFF": pk.parallel_equiv(pk.z_off(f * 1e9), f * 1e9), "ON": pk.parallel_equiv(pk.z_on(f * 1e9), f * 1e9)} for f in (5.2, 5.5, 5.8, 6.1, 6.4)}
    pred = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PREDICTED_BEFORE_HFSS_RUN"], "mesh_phase_vs_f": 64, "floquet_order": 128, "p_mm": P,
            "series_models": {"OFF": "(0.15 pF || 5 kohm) + 0.6 nH", "ON": "3.0 ohm + 0.6 nH"},
            "parallel_equivalents_5p8GHz": eq, "parallel_equivalents_vs_f": eq_f, "cells": {}, "runs": [],
            "method": "directional bound-mode condition: beta_d from the TM sheet admittance Y_dd(kt = beta_d) alone (d = x along the slot, y across), Y_dd interpolated linearly in kt between adaptive nodes; phase = beta_d * p; open air, no lid; f(phase) = quadratic fit of the 5 frequencies 5.2-6.4 GHz"}
    ids = {"CELL1": "RT5880_p11.6_g0.73_q0.8_W2.9", "CELL2": "RT5880_p11.6_g1_q0.9_W2.9"}
    fits = {}
    for name, c in CELLS.items():
        grid = next(r for r in rows if r["row"]["id"] == ids[name])
        entry = dict(params=c, grid_row_id=ids[name], mean_sheet_beta_over_k0_5p8GHz_mesh64={st: grid["states"][st]["beta_over_k0"] for st in ("NOCUT", "OFF_BAR50", "ON_BAR50")}, states={})
        if ids[name] in mesh:
            m = mesh[ids[name]]
            entry["mean_sheet_beta_over_k0_5p8GHz_mesh96"] = {st: m["beta96"][st] for st in ("NOCUT", "OFF_BAR50", "ON_BAR50")}
            entry["mesh96_note"] = "mesh 96 exists only at 5.8 GHz and only as the mean-sheet beta (grid definition); the directional phase-vs-f below is mesh 64"
        else:
            entry["mesh96_note"] = "no mesh-96 run for this cell (not among the 3 best rows): mesh 64 only"
        for st in ("NOCUT", "OFF_BAR50", "ON_BAR50"):
            for d in ("x", "y"):
                tr_all = tab[f"{name}|{st}|{d}"]
                tr = [r for r in tr_all if r["nearest_node_dist"] < 0.02]  # drop unconverged points (CELL2 ON y 6.4 GHz: beta not within 0.098 of a node, 4-node cap)
                dropped = [r["f_GHz"] for r in tr_all if r["nearest_node_dist"] >= 0.02]
                f = np.array([r["f_GHz"] for r in tr])
                ph = np.array([r["phase_deg"] for r in tr])
                fit = np.polyfit(ph, f, 2)
                fits[(name, st, d)] = fit
                phi58 = float(np.interp(5.8, f, ph))
                phic = int(5 * round(phi58 / 5))
                entry["states"][f"{st}|{d}"] = dict(
                    beta_over_k0_5p8GHz=float(np.interp(5.8, f, [r["beta_over_k0"] for r in tr])), phase_deg_5p8GHz=phi58,
                    phase_vs_f=[dict(f_GHz=float(a), beta_over_k0=r["beta_over_k0"], phase_deg=float(b)) for a, b, r in zip(f, ph, tr)],
                    f_at_phase_GHz={str(q): float(np.polyval(fit, q)) for q in PHASES_SCAN}, phase_range_of_fit_deg=[float(ph.min()), float(ph.max())],
                    recommended_phases_deg=[phic - 5, phic, phic + 5], max_nearest_node_dist=max(r["nearest_node_dist"] for r in tr), dropped_unconverged_f_GHz=dropped)
        pred["cells"][name] = entry
    pred["sanity_bare_substrate_GHz"] = {str(q): dict(open_air=bare_f(q, None), lid_150mm=bare_f(q, LID), lid_100mm=bare_f(q, 100.0)) for q in (80, 90, 100)}
    names = {"NOCUT": "NOCUT", "OFF_BAR50": "OFF", "ON_BAR50": "ON"}
    n = 0
    for name in CELLS:
        for st in ("NOCUT", "OFF_BAR50", "ON_BAR50"):
            for d in ("x", "y"):
                e = pred["cells"][name]["states"][f"{st}|{d}"]
                for i, q in enumerate(e["recommended_phases_deg"]):
                    n += 1
                    pred["runs"].append(dict(run_id=f"R{n:02d}", cell=name, state=names[st], direction="along_x" if d == "x" else "across_y", ph_x_deg=q if d == "x" else 0, ph_y_deg=0 if d == "x" else q,
                                             tier=1 if i == 1 else 2, f_MoM_GHz=float(np.polyval(fits[(name, st, d)], q))))
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    with open(OUT / "states.csv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["cell", "state", "step_file", "DIODE_SHEET", "R_ohm", "C_pF", "L_nH", "series_model", "note"])
        for name in CELLS:
            dx = CELLS[name]["diode_x"]
            w.writerow([name, "NOCUT", f"{name}_nocut.step", "none (no object)", "", "", "", "", "closed slot, no notch"])
            w.writerow([name, "OFF", f"{name}_notch.step", "Lumped RLC parallel", f"{eq['OFF']['R_ohm']:.4g}", f"{eq['OFF']['C_pF']:.4g}", "unticked", "(0.15 pF || 5 kohm) + 0.6 nH", f"5.8 GHz parallel equivalent; current line ({dx},-0.5,0.0175)->({dx},0.5,0.0175)"])
            w.writerow([name, "ON", f"{name}_notch.step", "Lumped RLC parallel", f"{eq['ON']['R_ohm']:.4g}", "unticked", f"{eq['ON']['L_nH']:.4g}", "3.0 ohm + 0.6 nH", "5.8 GHz parallel equivalent; same geometry as OFF"])
    with open(OUT / "results_template.csv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["run_id", "cell", "state", "direction", "ph_x_deg", "ph_y_deg", "tier", "f_MoM_GHz", "f_HFSS_GHz", "Q", "Ez_decays_from_patch_YN", "note"])
        for r in pred["runs"]:
            w.writerow([r["run_id"], r["cell"], r["state"], r["direction"], r["ph_x_deg"], r["ph_y_deg"], r["tier"], f"{r['f_MoM_GHz']:.3f}", "", "", "", ""])
    setup(pred, eq)
    print("written; runs", len(pred["runs"]))


def setup(pred, eq):
    c1, c2 = CELLS["CELL1"], CELLS["CELL2"]
    bare = pred["sanity_bare_substrate_GHz"]["90"]
    light90 = C0 * math.radians(90) / (P * 1e-3) / (2 * math.pi) / 1e9
    txt = f"""# HFSS eigenmode: SUB6 slot-disc cells CELL1 / CELL2 (RT5880 3.175 mm, p 11.6, 5.8 GHz). EXPORT STANDARD (2026-11-05). Predictions MODEL_ONLY (`predictions.json`), NOT_HFSS_VALIDATED.

**1. Assign material `pec` to `METAL_PEC` (patch + ground). STEP imports as vacuum; without this nothing works.** (If HFSS splits it into two parts, assign both.)
**2. Air box: `AIR` = one vacuum box z = -3.210 ... +{LID:g} mm, ONE face per side wall.** It overlaps `SUBSTRATE` and `METAL_PEC` on purpose (nothing subtracted); keep `AIR` first in the tree. **Periodic pairs go on the four side faces of `AIR` only** (x = -5.8 / +5.8 and y = -5.8 / +5.8).
**3.** `SUBSTRATE`: RT5880, eps_r 2.2, tan d 0.0009 (z -3.175 ... 0). **4. PEC lid:** Perfect E on the top face of `AIR` (z = +{LID:g}). The bound mode decays as exp(-kappa z), kappa = k0 sqrt((beta/k0)^2 - 1) = 39-56 /m for beta/k0 1.05-1.21, so a 150 mm lid leaves exp(-kappa z) of about 3e-3 (1.05) ... 2e-4 (1.21) (the 10 GHz packages used 30 mm). Lid effect on the BARE substrate (MoM transverse resonance, 90 deg): open air {bare['open_air']:.3f} GHz, lid 150 mm {bare['lid_150mm']:.3f} GHz, lid 100 mm {bare['lid_100mm']:.3f} GHz; the MoM predictions below are open air.
**5.** Project Variables `$ph_x`, `$ph_y` (deg). Lattice Pair 1: x faces of `AIR`, Phase Delay `$ph_x`; Pair 2: y faces, `$ph_y`. **Two run types: ALONG the slot (x): `$ph_x = phi`, `$ph_y = 0`; ACROSS the slot (y): `$ph_x = 0`, `$ph_y = phi`.**
**6. Files / states** (`states.csv`; the geometry differs only NOCUT vs notch):
| state | file | `DIODE_SHEET` (Lumped RLC, PARALLEL) |
|---|---|---|
| NOCUT | `CELL1_nocut.step` / `CELL2_nocut.step` | none (closed slot, no notch, no sheet) |
| OFF | `CELL1_notch.step` / `CELL2_notch.step` | **R = {eq['OFF']['R_ohm'] / 1e3:.3f} kohm, C = {eq['OFF']['C_pF']:.4f} pF, L unticked** |
| ON | same notch file | **R = {eq['ON']['R_ohm']:.1f} ohm, L = {eq['ON']['L_nH']:.4f} nH, C unticked** |

Parallel equivalents are exact at 5.8 GHz only (OFF series model (0.15 pF || 5 kohm) + 0.6 nH = {eq['OFF']['Z_series'][0]:.2f} {eq['OFF']['Z_series'][1]:+.2f}j ohm; ON 3.0 ohm + 0.6 nH = {eq['ON']['Z_series'][0]:.2f} {eq['ON']['Z_series'][1]:+.2f}j ohm). If an eigenfrequency lands more than 0.15 GHz from 5.8, the equivalents at that frequency are in `predictions.json` (`parallel_equivalents_vs_f`).
**Sheet `DIODE_SHEET`: 1.00 (across, y) x 0.35 (along, x) mm at z = +0.0175, centred in the rim between the slot end and the disc edge, long edges on the notch walls y = -0.5 / +0.5.** Current line across the notch: CELL1 ({c1['diode_x']:.3f}, -0.5, 0.0175) -> ({c1['diode_x']:.3f}, +0.5, 0.0175); CELL2 ({c2['diode_x']:.3f}, -0.5, 0.0175) -> ({c2['diode_x']:.3f}, +0.5, 0.0175).
**7.** Eigenmode: min frequency 4.5 GHz, 6 modes, delta f 0.1 %, max 15 passes; initial mesh 1 mm on `SUBSTRATE`, 0.5 mm on the patch faces, 0.1 mm on `DIODE_SHEET` and the notch. **Clean Up Solutions before each run.** Patch 35 um (z 0 ... 0.035), ground 35 um (z -3.210 ... -3.175).
**Mode: mode 1 (lowest, TM-like bound mode below the light line; the light line at 90 deg is {light90:.3f} GHz).** Record f and Q; check Ez is large over the patch and decays upward (a higher mode with Ez ~ 0 circulating around the patch is TE-like: do not use). Runs are in `results_template.csv` (also `predictions.json` `runs`): tier 1 = the phase where the model puts the mode at 5.8 GHz (rounded to 5 deg); tier 2 = +-5 deg, to bracket 5.8 GHz if tier 1 is more than about 0.15 GHz off.
**Sanity.** Bare substrate (ground only, no metal on top): model {bare['open_air']:.3f} GHz open air / {bare['lid_150mm']:.3f} GHz with the 150 mm lid at 90 deg. With `METAL_PEC` assigned the NOCUT run must be clearly lower than the bare substrate. **Expected null result:** ALONG x the model puts NOCUT, OFF and ON within 0.03 GHz of each other at the same phase (the notch/diode does not load the x-directed current); ACROSS y they differ by 0.1-0.7 GHz at the same phase (`predictions.json`). If HFSS shows OFF = ON across y, the sheet or the RLC is not active.
Geometry: disc r {c1['r']:.3f} / {c2['r']:.3f} mm, slot L {c1['L']:.3f} / {c2['L']:.3f} x W 2.900 mm along x, notch 1.00 mm from the slot end to the disc edge at the +x end (CELL1 / CELL2). Not run in AEDT by the author.
"""
    (OUT / "HFSS_SETUP.md").write_text(txt, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
