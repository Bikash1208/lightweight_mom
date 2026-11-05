"""Docs of the K2-LP board package: HFSS_SETUP.md, predictions.json, assign_rlc_groups.py (AEDT helper), then syncs the package. Reads _design_summary.json / geometry_checks.json written by
sidequest_board_k2_lp_package.py.     python scripts/research/sidequest_board_k2_lp_docs.py   (MODEL_ONLY; HFSS cost = scaling ESTIMATE, nothing here is an HFSS result)"""
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "hfss_exports/sidequest_board_k2_lp"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\board_K2_LP")
chk = json.loads((OUT / "geometry_checks.json").read_text())
dsum = OUT / "_design_summary.json"
d = json.loads(dsum.read_text()) if dsum.exists() else json.loads((OUT / "design_summary.json").read_text())
S = d["summary"]
F = 10.54
launch, conv, edge = 0.627, d["conv"], d["edge"]
states = []
for k, nm in enumerate("AB"):
    sd = S["states_detail"][k]
    states.append({"state": nm, "R_A_ohm": 2.44 if nm == "A" else 1e9, "R_B_ohm": 1e9 if nm == "A" else 2.44, "groups_ON": "GA" if nm == "A" else "GB",
                   "target_theta_deg": 25.0, "target_phi_deg": 0.0 if nm == "A" else 180.0, "peak_gain_dbi": sd["gain_dbi"], "peak_theta_deg": sd["peak_theta_deg"], "peak_phi_deg": sd["peak_phi_deg"],
                   "sll_db": sd["sll_db"], "xpd_db": sd["xpd_db"], "polarisation": "LP phi-hat (LP-V)",
                   "s11_db": -15.0,
                   "s11_note": "bare calibrated probe (Phase 26 board, HFSS S11 -15.0 dB, same stack and probe); the hologram's loading of the probe is NOT modelled",
                   "launch_fraction": launch, "sw_to_beam_fraction": round(conv, 3), "efficiency_of_accepted_power_into_beam": round(launch * conv, 3),
                   "edge_power_pct_of_surface_wave": round(100 * edge, 1), "aperture_efficiency_pct": S["ap_eff_pct"][k]})
pred = {"labels": ["MODEL_ONLY", "cell HFSS-validated (C/C 10.670, R/R 11.006 GHz at 90 deg; dY -3.13 mS, alpha 0.844 Np/m); board NOT HFSS-validated",
                   "gain = ideal-aperture directivity x launch 0.627 (bare probe) x surface-wave-to-beam conversion x scalar-cell element factor"],
        "design": {"frequency_ghz": F, "K": 2, "cells": 1430, "switched_cells": 360, "diode_sheets": 720, "static_open": 526, "static_strap": 544, "static_pct": S["static_pct"],
                   "lines": "1 (+ inverter): GA on the line, GB on its complement", "source": "results/sidequest_hybrid/scalar_board.json, row LP-V eps_B 180 bit lam 0.28 (720 diodes)"},
        "states": states,
        "falsifiers": "peak gain < 12 dBi, peak more than 3 deg from (25, 0) / (25, 180), SLL above -8 dB, or the A and B patterns not mirror images (phi 0 vs 180)",
        "NA": ["board S11 with the hologram", "HFSS gain", "direct radiation of the probe (not patterned)"]}
(OUT / "predictions.json").write_text(json.dumps(pred, indent=1), encoding="utf-8", newline="\n")

setup = f"""# HFSS_SETUP - K = 2 LP board (25, 0) / (25, 180), scalar rim-disc cell, 150 mm hexagon, {F} GHz. MODEL_ONLY (cell HFSS-validated, board not). One page.
**1. Assign material `pec` to `METAL_PEC` (patches, ground and probe: ONE named part, 1432 DC-isolated lumps; the STEP imports as vacuum).** `SUBSTRATE`: RT5880 eps_r 2.2, tan d 0.0009. `AIR`: vacuum, first in the tree (it overlaps the solids on purpose; later objects take priority).
**Files** (mm): `board_K2_LP.step` ({chk['step_mb']} MB: `METAL_PEC`, `SUBSTRATE`, `AIR`, `PORT`, 720 sheets `GA_<ix>_<iy>_DX|DY` / `GB_...`), `states.csv`, `sheets.csv` (current-line coordinates), `cells.csv`, `assign_rlc_groups.py` (AEDT helper, untested), `predictions.json`, `geometry_checks.json` (all PASS on the re-imported STEP).
**Geometry**: substrate hexagon (0, +-154), (+-133.368, +-77), 1.575 thick, probe hole r 0.635; ground 35 um (z -1.610..-1.575) with an opening r 2.05; 1430 cells on the 6.4 lattice (integer index x 6.4): disc r 3.0 with a cross slit 0.4 x 5.6 (rim 0.2), 35 um (z 0..0.035), cells chosen by their centres inside the 150 mm hexagon, no cell for rho < 6.5; probe r 0.635, z -1.575..+6.0 (6.0 above the substrate); `PORT` annulus r 0.635..2.05 at z -1.575.
**Cells**: 526 static open (+1), 544 static strapped (-1, metal in both footprints = a permanently ON diode with R -> 0), 360 switched. A switched cell has sheets DX (0, +0.8) and DY (+0.8, 0), 0.4 x 1.0 at z 0.0175, both in ONE group: `GA_*` is ON in state A, `GB_*` ON in state B.
**2. Boundaries**: Radiation on all six faces of `AIR` (z -1.610..+20, i.e. 14 mm above the probe top, > lambda/4 = 7.1 mm; x/y pad 10 mm; one face per side). PML only if licence and RAM allow; say which was used.
**3. Port**: `PORT` -> Lumped Port 50 ohm, integration line (0.635, 0, -1.575) -> (2.05, 0, -1.575), renormalise 50, no de-embedding.
**4. Variables** `R_A`, `R_B` (ohm). Lumped RLC on every `GA_*` sheet: Parallel, R = `R_A`, C = 28.16 fF, no L; on every `GB_*` sheet: R = `R_B`; current line per sheet from `sheets.csv` (DX along x, DY along y). **State A: R_A = 2.44, R_B = 1e9. State B: R_A = 1e9, R_B = 2.44.** One project, two solves, no re-import (`states.csv`).
**5. Solution**: ONE frequency **{F} GHz** (the design frequency of the model and of the Phase 26 board), Driven Modal, first-order basis, max passes 6, max delta S 0.02, no sweep. Initial mesh: substrate max length 1.0 mm, patch faces 0.4 mm, diode sheets 0.15 mm, probe and port 0.3 mm, lambda refinement on. Direct solver with HPC (see cost). Mesh reuse A -> B (no Clean Up) is allowed: say so, and S11 must agree within 0.5 dB between a reused and a fresh B.
**6. Export (per state, exactly this)**: far-field infinite sphere theta 0..180, phi 0..360 step **2 deg**, **complex rETheta and rEPhi in V** (not mV), one CSV per state (`rE_A.csv`, `rE_B.csv`); complex S11 at {F} GHz; antenna parameters (peak realized gain, directivity, radiation efficiency, accepted power); passes, delta S, tets, peak RAM, solve time per state.
**Model predictions (MODEL_ONLY)**: both states peak {S['gain'][0]:.2f} dBi at (25, 0) / (25, 180), SLL {S['SLL'][0]} dB, XPD {S['XPD'][0]} dB, aperture efficiency {S['ap_eff_pct'][0]} %, S11 -15 dB (bare probe only), {100 * launch * conv:.1f} % of the accepted power into the beam (launch {launch} x conversion {conv:.2f}), {100 * edge:.0f} % of the surface-wave power reaches the edge. Falsifiers: `predictions.json`.
**Cost (ESTIMATE, scaling only; the Phase 26 board's tets/RAM/time were not recorded = NA)**: `results/sidequest_hybrid/BOARD_K2_LP_HFSS_COST.md`; 1430 cells + 720 sheets, about 9-25 M tets, 100-300 GB, 20-50 h for state A on 32-64 cores, 7-16 h for state B.
"""
(OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
helper = '''"""AEDT (IronPython) helper - UNTESTED draft: assigns one parallel-RLC boundary per diode sheet (R = R_A for GA_*, R_B for GB_*, C 28.16 fF) from sheets.csv.
Run in HFSS: Tools > Run Script. Create the project variables R_A and R_B first. Check the first boundary by hand (current-line direction) before trusting the rest."""
import csv
import os

oDesign = oProject.GetActiveDesign()
oMod = oDesign.GetModule("BoundarySetup")
here = os.path.dirname(os.path.abspath(__file__))
for r in csv.DictReader(open(os.path.join(here, "sheets.csv"))):
    var = "R_A" if r["group"] == "GA" else "R_B"
    z = r["z_mm"] + "mm"
    oMod.AssignLumpedRLC(["NAME:" + r["sheet"], "Objects:=", [r["sheet"]],
                          ["NAME:CurrentLine", "Start:=", [r["current_line_start_x"] + "mm", r["current_line_start_y"] + "mm", z], "End:=", [r["current_line_end_x"] + "mm", r["current_line_end_y"] + "mm", z]],
                          "RLC Type:=", "Parallel", "UseResist:=", True, "Resistance:=", var, "UseInduct:=", False, "UseCap:=", True, "Capacitance:=", "28.16fF"])
'''
(OUT / "assign_rlc_groups.py").write_text(helper, encoding="utf-8", newline="\n")
print("HFSS_SETUP words", len(setup.split()))
if dsum.exists():
    dsum.rename(OUT / "design_summary.json")
SYNC.mkdir(parents=True, exist_ok=True)
for f in OUT.iterdir():
    shutil.copy2(f, SYNC / f.name)
print("synced", sorted(p.name for p in SYNC.iterdir()))
