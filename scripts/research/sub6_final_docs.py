"""SUB6 FINAL BOARD: regenerate HFSS_SETUP.md, predictions.json, results_template.csv for the board package from the synthesis output (final_design) and the layout summary.  PRELIMINARY_EXTRAPOLATED_CELL / MODEL_ONLY / NOT_HFSS_VALIDATED.
    python sub6_final_docs.py          (SUB6_MESH=96 for the mesh-96 design)
"""
import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sub6_final_design as FD
import sub6_final_hfss_package as PK
import sub6_final_layout as L

OUT = PK.OUT
MESH = os.environ.get("SUB6_MESH", "64")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fd = FD.final_design("cal")
    m = {k: (float(v) if not isinstance(v, str) else v) for k, v in fd["metrics"].items()}
    meta = json.load(open(L.OUT / "fab" / "layout_summary.json"))
    rob = json.load(open(L.OUT / "robustness.json")) if (L.OUT / "robustness.json").exists() else {"required_pass": False}
    clr_ok = bool(meta["verification"].get("all_ge_0p15", False))
    final = bool(rob["required_pass"] and clr_ok and MESH == "96")
    STATUS = "FINAL_CANDIDATE" if final else "PRELIMINARY"
    for fn in (L.OUT / "fab" / "layout_summary.json", OUT / "geometry_checks.json"):
        if fn.exists():
            t = fn.read_text(encoding="utf-8").replace("PRELIMINARY_EXTRAPOLATED_CELL / MODEL_ONLY", STATUS + " / MODEL_ONLY")
            fn.write_text(t, encoding="utf-8", newline="\n")
    pred = dict(label=STATUS + " / MODEL_ONLY / NOT_HFSS_VALIDATED", robustness_required_pass=rob["required_pass"], min_clearances_mm=meta["verification"].get("min_clearances_mm"), mesh_of_switched_cell=MESH, N_switched=240, gates=dict(balance_dB_le=3, AR_dB_le=3, B_suppression_dB_ge=12, SLL_dB_le=-13), model_metrics=m, phase_definition="RHCP/LHCP at theta 25, phi 0 (beam A) and 180 (beam B)", strings=meta["strings"])
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1), encoding="utf-8", newline="\n")
    with open(OUT / "results_template.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["state", "freq_GHz", "S11_dB", "gain_RHCP_beamA_dBic", "gain_RHCP_beamB_dBic", "axial_ratio_A_dB", "axial_ratio_B_dB", "SLL_dB", "radiation_efficiency", "note"])
        for st in ("OFF", "ON"):
            for fq in (5.6, 5.8, 6.0):
                w.writerow([st, fq, "", "", "", "", "", "", "", "model: see predictions.json" if fq == 5.8 else ""])
    md = f"""# HFSS_SETUP.md: SUB6 final board, N 240 ({STATUS} / MODEL_ONLY / NOT_HFSS_VALIDATED; switched-cell table = mesh {MESH})

**STATUS: {STATUS}.** Switched-cell table = mesh {MESH} (AEM-LAB1, compute/sub6-island-meshcheck; NODEFIX 566f8fe: the kt-node extrapolation changed the values by only ~0.1 %, no re-run). FINAL_CANDIDATE requires: all trace/copper clearances >= 0.15 mm ({clr_ok}) and the +-15 % alpha robustness (joint x0.85 / x1.15, no re-optimisation) passing all gates ({rob['required_pass']}). Not validated by HFSS; mesh 128 (if it changes the cell by > 10 %) would need a re-run.

## Regeneration (all outputs are script-generated)
`SUB6_MESH=<64|96> python sub6_final_cal.py run` (steps 1-2), `python sub6_final_export.py` (step 3 routing + DXF, drill, cell table, string table), `python sub6_final_hfss_package.py` (board STEP), `python sub6_final_docs.py` (this file, predictions, template).

## Files
BOARD_N240.step (bodies METAL_PEC, SUBSTRATE, SUBSTRATE_FR4, AIR, DIODE_SHEET, PORT_ANNULUS), diode_sheets.csv (240 sheets: centre, current direction, string), predictions.json, results_template.csv.
Fabrication (results/sub6_final_board/): fab/sub6_final_board.dxf, fab/drill.drl/.csv, fab/cell_table.csv, fab/string_table.csv.

## Geometry (mm, board frame, feed at the origin, +z up)
Top copper 35 um at z 0 ... 0.035 (758 slotted discs, 240 with islands); RT5880 (eps_r 2.2, tan d 0.0009) z -3.175 ... 0; ground z -3.210 ... -3.175 with a 0.5 mm radius anti-pad per via and a 2.05 mm opening at the probe;
FR4 0.8 (eps_r 4.4, tan d 0.02 ASSUMED typical) z -4.010 ... -3.210; bottom copper 35 um z -4.045 ... -4.010 (bias traces 0.15 mm, via pads 0.3 mm radius, 24 connector pads). Vias: drill 0.3 (barrel r 0.15) from the top copper to the bottom pad.
Probe: PEC pin r 0.635 from z -4.045 to z +{PK.Z_PIN_TOP} (12.1 mm above the ground top, no hat), PORT_ANNULUS sheet r 0.635 ... 2.05 at z {PK.Z_FR4B} (coax opening). Air box: hexagon bounding box + 15 mm laterally, 15 mm below, 40 mm above the top copper.

## HFSS steps
1. Import BOARD_N240.step (units mm). Assign materials: SUBSTRATE RT5880 as above, SUBSTRATE_FR4 FR4 (eps_r 4.4, tan d 0.02), AIR vacuum, METAL_PEC -> Perfect E (or copper 35 um with surface roughness off; PEC is the model assumption).
2. PORT_ANNULUS: wave port (coax, Z0 from the pin and opening radii; renormalise to 50 ohm) with the integration line from the pin to the opening edge. If the wave port at the FR4 plane is not accepted use a lumped port on the same sheet and say so in the results.
3. AIR faces: radiation boundary on all outer faces (>= lambda/4 at 5.8 GHz = 12.9 mm from the board: 15 mm used).
4. DIODE_SHEET (240 sheets): lumped RLC boundary with the current direction from diode_sheets.csv (current_dir_x, y per sheet, along the local +y of each cell, i.e. across the 1.00 mm diode gap). The sheet is 1.00 x 0.35 mm.
   OFF set: R 3876 ohm parallel C 0.1703 pF (all 240 sheets).   ON set: R 162.4 ohm parallel L 0.6113 nH (all 240 sheets).   (parallel equivalents at 5.8 GHz of the BAR50 model)
   Bias traces, anti-pads and vias are geometry only: no bias network or DC block is modelled.
5. Run order: OFF first (all sheets RLC OFF), then ON (re-assign the same sheets to the ON boundary); same mesh for both (Setup, adaptive solution at 5.8 GHz, max delta S 0.02, 2-3 passes; lambda-based mesh refinement of the substrate and the top copper layer).
   Discrete sweep 5.6 / 5.8 / 6.0 GHz optional; far field at 5.8 GHz.
6. Far field: infinite sphere, theta 0 ... 90 deg step 1 deg, phi 0 ... 360 deg step 1 deg (cuts phi 0 and 180 are beams A and B at theta 25 deg). Use the standard HFSS GainRHCP / GainLHCP quantities. Report GainRHCP at (25, 0) and (25, 180) in OFF and ON, axial ratio at those angles, SLL, and the radiated efficiency; fill results_template.csv.

## Model prediction (MODEL_ONLY, N 240, calibrated, mesh {MESH} switched cell)
OFF: beam A {m['off_GA']:.2f} dBic, beam B {m['off_GB']:.2f} dBic, balance {m['bal']:.2f} dB, AR(-3 dB beam) A/B {m['off_ARA3']:.2f}/{m['off_ARB3']:.2f} dB, SLL {m['sll_off']:.1f} dB, eta {m['eta_off']:.2f}.
ON: beam A {m['on_GA']:.2f} dBic, B suppression {m['supp']:.1f} dB, AR {m['on_ARA3']:.2f} dB, SLL {m['sll_on']:.1f} dB, eta {m['eta_on']:.2f}.   Gate margin {m['margin']:.2f} dB: {m['gate']}.
Expect HFSS to differ (MoM leakage only, directional beta approximation, bias network not modelled). The board STEP has no tolerance exploration; run OFF first.
"""
    (OUT / "HFSS_SETUP.md").write_text(md, encoding="utf-8", newline="\n")
    print("docs written", OUT)


if __name__ == "__main__":
    main()
