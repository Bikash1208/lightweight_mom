"""Docs of the K2-default board packages: HFSS_SETUP.md, predictions.json, assign_impedance_D.py (AEDT helper, untested), HFSS cost note; then syncs the package.
    python sidequest_board_k2_default_docs.py <outname> <syncname> <R_mm>      (MODEL_ONLY; the cost is a scaling ESTIMATE, nothing here is an HFSS result)"""
import csv
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTNAME, SYNCNAME, R_MM = sys.argv[1], sys.argv[2], int(sys.argv[3])
OUT = ROOT / "hfss_exports" / OUTNAME
SYNC = Path(r"C:\hmsa_sync\data\sidequest") / SYNCNAME
chk = json.loads((OUT / "geometry_checks.json").read_text(encoding="utf-8"))
ds = json.loads((OUT / "design_summary.json").read_text(encoding="utf-8"))
S = ds["summary"]
n_cells, n_sw = ds["n_cells"], S["switched_cells"]
n_sheets = 2 * n_sw
RS = 2.44 * 1.0 / 0.4
F = 10.5
LAUNCH = 0.627
ISO = S["isolation"]
states = []
for k, nm in enumerate(("1", "2")):
    states.append({"state": nm, "name": "default (unpowered, diodes OFF)" if nm == "1" else "circuit ON (diodes ON)", "D_sheets": "no boundary" if nm == "1" else f"one Impedance boundary, Rs = {RS} ohm/sq",
                   "target_theta_deg": 25.0, "target_phi_deg": 0.0 if nm == "1" else 180.0, "peak_gain_dbi": S["gain_dbi"][k], "peak_dir_theta_phi_deg": S["peak_dir"][k], "sll_db": S["sll_db"][k], "xpd_db": S["xpd_db"][k],
                   "aperture_efficiency_pct": S["ap_eff_pct"][k], "polarisation": "LP phi-hat (LP-V)" if S["pol"] == "LP-V" else "LP theta-hat (LP-H)",
                   "other_beam_level_dbi": ISO["beam2_in_state1_dbi"] if nm == "1" else ISO["beam1_in_state2_dbi"], "other_beam_rel_to_own_peak_db": ISO["beam2_in_state1_rel_db"] if nm == "1" else ISO["beam1_in_state2_rel_db"],
                   "s11_db": -15.0, "s11_note": "bare calibrated probe (Phase 26 board, HFSS S11 -15.0 dB, same stack and probe); the hologram's loading of the probe is NOT modelled"})
pred = {"labels": ["MODEL_ONLY", "cell HFSS-validated for the switching (C/C 10.670, R/R 11.006 GHz at 90 deg; dY -3.13 mS, alpha 0.844 Np/m, beta/k0 1.0731 back-out at 10.54 GHz, not re-derived at 10.5); analog static cells and the board are NOT HFSS-validated",
                   "gain = ideal-aperture directivity x launch 0.627 (bare probe) x surface-wave-to-beam conversion x scalar-cell element factor"],
        "design": {"id": "(i) family A only, one group D_*", "frequency_ghz": F, "hexagon_circumradius_mm": R_MM, "cells": n_cells, "switched_cells": n_sw, "diode_sheets": n_sheets, "diodes": S["diodes"], "static_pct_cells": S["static_pct_cells"],
                   "cells_by_type": ds["cells_by_type"], "suppression_weight_mu": S["mu"], "pol": S["pol"], "eps_B_deg": S["eps_B_deg"], "lam": S["lam"], "g": S["g"]},
        "states": states,
        "falsifiers": "peak gain more than 4 dB below the prediction, peak more than 3 deg from (25, 0) / (25, 180), the other beam within 8 dB of the own beam (predicted -15 dB), or the two patterns not mirror images",
        "NA": ["board S11 with the hologram", "HFSS gain", "direct radiation of the probe (not patterned)", "tets / RAM / time of the Phase 26 150 mm board (not recorded)"]}
(OUT / "predictions.json").write_text(json.dumps(pred, indent=1), encoding="utf-8", newline="\n")

# ------------------------------------------------------------------ cost (ESTIMATE): scaled from the earlier K2_LP estimate (1430 cells + 720 sheets: 9-25 M tets, 100-300 GB, 20-50 h state A) and the 90 mm test board (572 cells, 5+ h)
eq = lambda c, s: c + 0.25 * s  # noqa: E731  (sheets weigh a quarter of a cell in the tet count)
r_t = eq(n_cells, n_sheets) / eq(1430, 720)
tets = (9 * r_t, 25 * r_t)
ram = (10 * tets[0], 15 * tets[1])
t_a = tuple(5.0 * (n_cells / 572) ** e * f for e, f in ((1.3, 1.2), (1.7, 2.0)))
cost = {"tets_M": [round(x, 1) for x in tets], "ram_gb": [round(x) for x in ram], "time_state1_h": [round(x) for x in t_a], "time_state2_h_mesh_reused": [round(x / 3) for x in t_a], "fits_256gb": bool(ram[1] <= 256), "cells": n_cells, "sheets": n_sheets}
(OUT / "hfss_cost.json").write_text(json.dumps(cost, indent=1), encoding="utf-8", newline="\n")

setup = f"""# HFSS_SETUP - K2-default board (i): default = beam 1 (25, 0), circuit ON = beam 2 (25, 180), ONE diode group `D_*`, {S['pol']} linear, exactly {F} GHz, {R_MM} mm hexagon circumradius. MODEL_ONLY. One page.
**1. Assign material `pec` to `METAL_PEC` (patches, ground and probe: ONE named part, {chk['checks']['metal_single_part_lumps']['lumps']} DC-isolated lumps; the STEP imports as vacuum).** `SUBSTRATE`: RT5880 eps_r 2.2, tan d 0.0009. `AIR`: vacuum, first in the tree (overlaps the solids on purpose; later objects take priority).
**Files** (mm): `board_K2_default.step` ({chk['step_mb']} MB: `METAL_PEC`, `SUBSTRATE`, `AIR`, `PORT`, {n_sheets} sheets `D_<ix>_<iy>_DX|DY`), `states.csv`, `sheets.csv` (current-line coordinates), `cells.csv`, `assign_impedance_D.py` (AEDT helper, untested), `predictions.json`, `hfss_cost.json`, `geometry_checks.json` (all PASS on the re-imported STEP).
**Geometry**: substrate hexagon circumradius {R_MM} (vertices (0, +-{R_MM}), (+-X, +-{R_MM / 2:g}), X rounded to 1 um), 1.575 thick, probe hole r 0.635; ground 35 um (z -1.610..-1.575) with an opening r 2.05; {n_cells} cells on the 6.4 lattice (integer index x 6.4), discs r 3.0 (35 um) fully inside the outline (>= 0.2 mm to the edge), no cell for rho < 6.5; probe r 0.635, z -1.575..+6.0; `PORT` annulus r 0.635..2.05 at z -1.575.
**Cells**: {n_sw} switched cells (rim disc, cross slit 0.4 x 5.6, rim 0.2) with two sheets each, DX (0, +0.8) and DY (+0.8, 0), 0.4 x 1.0 at z 0.0175, ALL in the one group `D_*`; {n_cells - n_sw} static analog cells (disc with a cross slit of half-length 0 / 1.0 / 1.8 / 2.2 / 2.6 / 2.8, or with the footprints strapped): `cells.csv`. Static cells carry no sheet.
**2. Boundaries**: Radiation on all six faces of `AIR` (z -1.610..+20, i.e. 14 mm above the probe top, > lambda/4 = 7.14 mm; x/y pad 10 mm; one face per side). PML only if licence and RAM allow; say which was used.
**3. Port**: `PORT` -> Lumped Port 50 ohm, integration line (0.635, 0, -1.575) -> (2.05, 0, -1.575), renormalise 50, no de-embedding.
**4. The two states (same mesh, one project)**: **State 1 (default, diodes OFF): NO boundary on any `D_*` sheet (OFF = open).** **State 2 (circuit ON): ONE Impedance boundary assigned to ALL `D_*` sheets, Rs = {RS} ohm/sq, reactance 0** (Rs = R_on x width / gap = 2.44 ohm x 1.0 mm / 0.4 mm from the real sheet size, verified on the re-imported STEP). No RLC boundaries, no variables.
**5. Solution**: ONE frequency **{F} GHz**, Driven Modal, first-order basis, max passes 6, max delta S 0.02, no sweep. Initial mesh: substrate max length 1.0 mm, patch faces 0.4 mm, `D_*` sheets 0.15 mm, probe and port 0.3 mm, lambda refinement on. Direct solver with HPC (see cost). Mesh reuse state 1 -> 2 (no Clean Up) is allowed: say so; S11 of a fresh state 2 must agree within 0.5 dB.
**6. Export (per state, exactly this)**: far-field infinite sphere theta 0..180, phi 0..360 step **2 deg**, **complex rETheta and rEPhi in V** (full sphere), one CSV per state (`rE_state1.csv`, `rE_state2.csv`); complex **S11** at {F} GHz; antenna parameters (peak realised gain, directivity, radiation efficiency, accepted power); passes, delta S, tets, peak RAM, solve time per state.
**Model predictions (MODEL_ONLY)**: state 1 peak {S['gain_dbi'][0]} dBi at (25, 0), SLL {S['sll_db'][0]} dB, XPD {S['xpd_db'][0]} dB; state 2 {S['gain_dbi'][1]} dBi at (25, 180), SLL {S['sll_db'][1]} dB, XPD {S['xpd_db'][1]} dB; **isolation**: beam 2 in state 1 {ISO['beam2_in_state1_dbi']} dBi ({ISO['beam2_in_state1_rel_db']} dB), beam 1 in state 2 {ISO['beam1_in_state2_dbi']} dBi ({ISO['beam1_in_state2_rel_db']} dB); S11 -15 dB (bare probe only). Falsifiers: `predictions.json`.
**Cost (ESTIMATE)**: {cost['tets_M'][0]}-{cost['tets_M'][1]} M tets, {cost['ram_gb'][0]}-{cost['ram_gb'][1]} GB, {cost['time_state1_h'][0]}-{cost['time_state1_h'][1]} h state 1 on 32-64 cores, {cost['time_state2_h_mesh_reused'][0]}-{cost['time_state2_h_mesh_reused'][1]} h state 2 (mesh reused): `hfss_cost.json`; fits 256 GB: {cost['fits_256gb']}.
"""
(OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
names = [r["sheet"] for r in csv.DictReader(open(OUT / "sheets.csv", encoding="utf-8"))]
helper = f'''"""AEDT (IronPython) helper - UNTESTED draft: ONE Impedance boundary (Rs = {RS} ohm/sq, reactance 0) on all D_* sheets = state 2. State 1 = delete this boundary (no boundary on the sheets).
Run in HFSS: Tools > Run Script. Check the boundary after the run (Objects list, Rs)."""
oDesign = oProject.GetActiveDesign()
oMod = oDesign.GetModule("BoundarySetup")
names = {json.dumps(names)}
oMod.AssignImpedance(["NAME:Imp_D", "Objects:=", names, "Resistance:=", "{RS}", "Reactance:=", "0", "InfiniteGround:=", False])
'''
(OUT / "assign_impedance_D.py").write_text(helper, encoding="utf-8", newline="\n")
SYNC.mkdir(parents=True, exist_ok=True)
for f in OUT.iterdir():
    shutil.copy2(f, SYNC / f.name)
print(OUTNAME, "words", len(setup.split()), cost, "synced", sorted(p.name for p in SYNC.iterdir()))
