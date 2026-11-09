"""SUB6 FINAL BOARD, 13-pad variant: README_FAB.md, HFSS_SETUP.md, predictions.json, results_template.csv from the gate files.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python sub6_final_13pad_docs.py
"""
import csv
import json
import os
import shutil
import sys
from pathlib import Path

os.environ["SUB6_VARIANT"] = "13pad"
os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sub6_final_design as FD
import sub6_final_layout as L

ROOT = L.ROOT
D = L.OUT                                  # results/sub6_final_board/mesh96_13pad
FAB = D / "fab"
H = ROOT / "hfss_exports" / "sub6_final_board_m96_13pad"
SRC = ROOT / "results" / "sub6_final_board" / "mesh96"


def jl(p, default=None):
    try:
        return json.load(open(p))
    except Exception:
        return default


def main():
    H.mkdir(parents=True, exist_ok=True)
    conn, drc, rf = jl(FAB / "connectivity.json", {}), jl(FAB / "drc.json", {}), jl(FAB / "rf_unchanged.json", {})
    thick, thin = jl(H / "step_checks_thick.json", {}), jl(H / "step_checks_thin.json", {})
    gb_t, gb_n = jl(H / "geometry_build_thick.json", {}), jl(H / "geometry_build_thin.json", {})
    rob = jl(SRC / "robustness.json", {"required_pass": False})
    fd = FD.final_design("cal")
    m = {k: (float(v) if not isinstance(v, str) else v) for k, v in fd["metrics"].items()}
    kicad = jl(D / "kicad" / "kicad_gates.json", {})
    gates = dict(CONNECTIVITY=conn.get("state"), DRC=drc.get("state"), RF=rf.get("state"), KICAD_DRC=kicad.get("drc"), NO_OVERLAP_thick=thick.get("NO_OVERLAP"), NO_OVERLAP_thin=thin.get("NO_OVERLAP"), SHEETS_thick=thick.get("SHEETS"), SHEETS_thin=thin.get("SHEETS"))
    ok = all(v in ("CONNECTIVITY_OK", "DRC_OK", "RF_UNCHANGED", "NO_OVERLAP", "SHEETS_OK", "KICAD_DRC_0_ERRORS") for v in gates.values()) and rob.get("required_pass")
    STATUS = "FINAL_CANDIDATE_13PAD" if ok else "PRELIMINARY_13PAD"
    pins = list(csv.DictReader(open(FAB / "connector_pinout.csv")))
    strs = list(csv.DictReader(open(FAB / "string_table.csv")))
    mc = drc.get("min_clearances_mm", {})
    pred = dict(label=STATUS + " / MODEL_ONLY / NOT_HFSS_VALIDATED", gates=gates, mesh_of_switched_cell="96", N_switched=240, gate_requirements=dict(balance_dB_le=3, AR_dB_le=3, B_suppression_dB_ge=12, SLL_dB_le=-13), model_metrics=m,
                note="RF (cells, vias, anti-pads, ground, probe, diode sites) unchanged by the 13-pad bias re-route; the model metrics are those of the mesh-96 design (the bias network is not modelled).", strings=strs)
    (H / "predictions.json").write_text(json.dumps(pred, indent=1), encoding="utf-8", newline="\n")
    with open(H / "results_template.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["state", "freq_GHz", "S11_dB", "gain_RHCP_beamA_dBic", "gain_RHCP_beamB_dBic", "axial_ratio_A_dB", "axial_ratio_B_dB", "SLL_dB", "radiation_efficiency", "note"])
        for st in ("OFF", "ON"):
            for fq in (5.6, 5.8, 6.0):
                w.writerow([st, fq, "", "", "", "", "", "", "", "model: see predictions.json" if fq == 5.8 else ""])
    z = (gb_t.get("meta") or {}).get("z", {})
    ab = (gb_t.get("meta") or {}).get("air_box", {})
    zn = (gb_n.get("meta") or {}).get("z", {})
    pin_tbl = "\n".join(f"| {r['pad_index_from_top']} | {r['net']} | {r['string']} | {r['x_mm']} | {r['y_mm']} | {r['voltage_at_10mA_V']} | {r['current_mA']} |" for r in pins)
    # ------------------------------------------------------------ README_FAB.md
    readme = f"""# Fabrication README: SUB6 final board, 13-pad bias layout ({STATUS} / MODEL_ONLY / NOT_HFSS_VALIDATED)
RF is unchanged from the 24-pad board (cells, discs, islands, diode sites, vias, anti-pads, ground, probe): {rf.get('state')}. Only the B.Cu bias routing and the edge connector changed.
The B.Cu bias layer was routed with KiCad 10.0.8 + Freerouting 2.5.0 (B.Cu only, 0/45 degree segments), checked by KiCad DRC and by an independent script.

## Stack-up (top to bottom)
| layer | material | thickness mm |
|---|---|---|
| L1 F.Cu | 35 um Cu: discs, slots, islands, diode lands (SC-79 0.35 x 0.35, 1.35 pitch) | 0.035 |
| core | Rogers RT5880 (eps_r 2.2, tan d 0.0009) | 3.175 |
| L2 In1.Cu | 35 um Cu ground, anti-pad 1.0 mm at every via, 4.1 mm probe opening | 0.035 |
| dielectric | FR4 (eps_r 4.4, tan d 0.02 ASSUMED) | 0.8 |
| L3 B.Cu | 35 um Cu bias traces, via pads, 13 connector pads | 0.035 |
(KiCad file: 4 copper layers because KiCad needs an even count; In2.Cu is empty and not fabricated.) No top solder mask, no silkscreen.
Via drill 0.30, pad 0.60, anti-pad 1.00 (480 vias), probe hole 1.27. Min trace/gap 0.15.

## Connector (13 pads, 1.0 x 1.0 mm, pitch 2.54, B.Cu, near the left flat edge, x = {drc.get('connector_x', -179.705)} mm; RETURN in the middle)
| pad (from top) | net | string | x mm | y mm | V at 10 mA | I mA |
|---|---|---|---|---|---|---|
{pin_tbl}
Each drive net S<string> feeds the ANODE end of its string (20 diodes in series, about 19.0 V at 10 mA); the 12 cathode ends (last diode of each string) are joined on B.Cu into the single net RETURN (12 x 10 mA = 120 mA).
Net naming: S<string> follows the string id of diode_sheets.csv (D<string>_<pos>); the pad order along the connector is a planar fan order (not 1..12): use the table.

## Files
fab/gerber_kicad/ (kicad-cli export of the routed board): sub6_board_routed-F_Cu.gtl, -In1_Cu.g1, -B_Cu.gbl, -F_Paste.gtp, -Edge_Cuts.gm1, sub6_board_routed.drl (0.30 plated vias + 1.27 probe hole, one Excellon file, MixedPlating).
fab/gerber/ (script-generated RF layers, byte-identical to the 24-pad files of 526f98de): sub6_final-F_Cu.gtl, sub6_final-In1_Cu.g2, sub6_final-PTH.drl, sub6_final-Edge_Cuts.gm1. Use these for F_Cu / In1_Cu / drill if byte-identical RF is wanted: the KiCad F_Cu / In1_Cu are the same geometry with circles written as 2 degree polygons (max deviation 0.0008 mm).
fab/assembly/: F_Paste (480 lands), F_Assy, diode_placement.csv, DXF/PDF assembly drawing, README_ASSEMBLY.md (polarity: anode on the island, cathode on the disc, current along the string).
fab/sub6_final_board.dxf, cell_table.csv, string_table.csv, connector_pinout.csv, drill.csv, B_Cu_full.png, B_Cu_connector_zoom.png; kicad/sub6_board_routed.kicad_pcb (+ .kicad_pro, .kicad_dru) editable in KiCad >= 10.

## Rules and results
Min clearances (copper edge to copper edge, mm): trace-trace {mc.get('trace_to_trace_gap_mm')}, trace to foreign anti-pad edge {mc.get('trace_edge_to_foreign_antipad_edge_mm')} (required 0.172), trace to foreign via pad {mc.get('trace_edge_to_foreign_via_pad_edge_mm')}, trace to connector pad {mc.get('trace_edge_to_foreign_connector_pad_edge_mm')}, trace to board outline {mc.get('trace_edge_to_board_outline_mm')}, trace to probe opening {mc.get('trace_edge_to_probe_opening_edge_mm')} (all >= 0.15).
Widths: signal 0.25 mm (short neck-downs to {drc.get('min_trace_width_mm')} mm near pads where the clearance does not allow 0.25), RETURN {drc.get('RETURN_widths_mm')} segments of 0.5 / 0.374 mm (min {drc.get('min_RETURN_width_mm')} mm, {drc.get('RETURN_length_by_width_mm')} mm of length). 35 um Cu carries about 0.8 A at 0.15 mm width (IPC-2221, 20 K rise), so 120 mA needs no minimum width; the 0.5 mm RETURN is a margin.

## SMA connector (flat on the bottom, no legs)
Body 6.35 x 6.35 mm soldered to a B.Cu ground pad 7.35 x 7.35 mm (body + 0.5 mm per side, opening r 2.05 around the pin); bias keep-out 10.35 x 10.35 mm; 12 blind GND vias (B.Cu to In1.Cu, FR4 only, R 2.7 mm, pitch 1.41 mm) in fab/gerber_kicad/sub6_board_routed-in1-back.drl (T1 0.30, layer pair In1-B.Cu) and the main .drl; recommendation and the through-via alternative: see HFSS_SETUP.md. Nearest bias trace edge {drc.get('min_trace_edge_chebyshev_distance_to_origin_mm')} mm from the centre (keep-out half size 5.175).

## Points to check with the manufacturer (not resolved here)
1. Via aspect ratio: 0.30 mm drill through about 4.08 mm = 13.6:1, above the usual 8 to 10:1 plating limit; mixed RT5880 / FR4 stack with a bonded ground layer. Ask for feasibility or a larger drill.
2. The RETURN bus reaches the middle pad from the board-edge side (it passes the pads S12..S08 on the edge side): a single-layer fan cannot reach a middle RETURN pad otherwise. With RETURN at the first or last pad the bus would enter directly.
3. BAR50-type diode assumed for the model; check the SC-79 land pattern of the real part (0.35 x 0.35 at 1.35 pitch used). No solder mask / silkscreen / paste layers beyond F_Paste are generated; FR4 properties assumed.
"""
    (FAB / "README_FAB.md").write_text(readme, encoding="utf-8", newline="\n")
    # ------------------------------------------------------------ HFSS_SETUP.md
    gb_line = lambda g: f"{g.get('n_bodies')} bodies, {g.get('step_size_MB')} MB, 0 inexact literals ({g.get('audit_inexact_literals')}), duplicate names {g.get('duplicate_names')}"
    H_SMA = 6.0
    sma = f"""
## SMA connector (flat on the FR4 side, NO legs; user add-on v2)
Assumed dimensions (no part number given): body 6.35 x 6.35 mm, height {H_SMA} mm below the board (ASSUMED), centre pin 1.27 mm (the existing probe hole, unchanged), flat body soldered to the B.Cu ground pad.
B.Cu ground pad: solid square 7.35 x 7.35 mm (body + 0.5 mm per side), net GND, with a circular clearance r 2.05 mm around the pin (same as the opening in the In1.Cu ground). Bias keep-out on B.Cu: square 10.35 x 10.35 mm (pad + 1.5 mm each side); nearest bias trace edge {drc.get('min_trace_edge_chebyshev_distance_to_origin_mm')} mm from the centre (keep-out half size 5.175 mm: clearance {drc.get('trace_edge_to_SMA_keepout_mm')} mm). The previous route had one drive lead (S06) 5.07 mm from the centre, inside the new keep-out: it was re-routed (Freerouting, whole B.Cu re-run), nothing else crossed it.
Ground tie: ring of 12 vias (0.30 drill / 0.60 pad) at radius 2.7 mm, pitch 1.41 mm (<= 2 mm), pad -> In1.Cu ground.
Evaluation of the via type: (a) FR4-only BLIND vias B.Cu -> In1.Cu (depth 0.87 mm, aspect 2.9:1; FR4 core drilled and plated before bonding to the RT5880 core, i.e. a sequential-lamination step, or controlled-depth drilling to the 35 um ground layer; extra process step, cost roughly +30 to 60 % on a prototype run, the ground layer must be the target pad of the blind via: ask the fabricator); the RT5880 stays untouched, RF unchanged. (b) THROUGH vias to the top: standard drilling, cheapest, but 12 holes through the RT5880 and a PEC via fence at r 2.7 mm around the feed pin (a coax-like shield of the probe region, changes the feed environment; needs an HFSS check). Top-layer items hit by (b): NONE: the ring (r 2.7 + pad 0.3 = 3.0 mm from the centre) is 3.0 mm inside the nearest top copper (the cells nearest to the centre have their discs 6.0 mm from it; nearest top vias / diodes are > 30 mm away). **Recommendation: (a) blind vias** (RF identical to the 24-pad design). No cell was moved. The KiCad board and both STEP files contain (a); the variant (b) is the same script with SMA_VIA=through.
In the STEP the SMA body is a PEC block 6.35 x 6.35 x {H_SMA} mm with a coax bore r 2.05 mm (fused with the B.Cu pad, the 12 via rods and the ground plate into METAL_PEC_GROUND), the centre pin runs through the bore to the bottom of the block, and the coax port sheet PORT_ANNULUS (pin r 0.635 to r 2.05) stays at the FR4 bottom plane as before. The bore is open at the bottom (an open coax stub of {H_SMA} mm below the port plane, below cut-off); the SMA dielectric and the real connector interior are not modelled.
"""
    md = f"""# HFSS_SETUP.md: SUB6 final board, N 240, 13-pad bias layout ({STATUS} / MODEL_ONLY / NOT_HFSS_VALIDATED; switched-cell table = mesh 96)

**STATUS: {STATUS}.** Gates: {json.dumps(gates)}.
The model contains the RF structure and the B.Cu bias traces (as PEC copper); no bias network / DC block / diode package is modelled. RF is unchanged from the previous (24-pad) package; the intersection errors of the previous STEP are fixed (see Bodies).

## Files (this folder)
BOARD_N240_clean.step (thick: 35 um copper solids), BOARD_N240_thin_clean.step (zero-thickness copper and ground as sheets), hfss_board_setup.py, set_state.py, diode_sheets.csv (240 sheets), geometry_build_thick/thin.json, step_checks_thick/thin.json, predictions.json, results_template.csv.
Both STEP files: thick {gb_line(gb_t)}; thin {gb_line(gb_n)}.

{sma}
## Which STEP to use
**Recommended: BOARD_N240_thin_clean.step.** Same plan coordinates, 4 x smaller file and a much lighter mesh (no 35 um solids: copper, ground and traces are Perfect-E sheets). Use the thick STEP only to check the effect of the copper thickness (expected small at 35 um).
Thin stack (contiguous): RT5880 z -3.175 ... 0, ground sheet z -3.175, FR4 z -3.975 ... -3.175 (0.8), bottom-copper sheets z -3.975, top-copper sheets z 0, diode sheets z 0, port sheet z -3.975. Barrels stay 0.3 mm rods (holes in both dielectrics), pin r 0.635 to z {zn.get('pin_top')}.
Thick stack: top Cu z 0 ... 0.035, RT5880 -3.175 ... 0, ground -3.210 ... -3.175, FR4 -4.010 ... -3.175 (the 0.035 band of the ground plane is FR4 outside the ground plate, so the 480 anti-pad voids are filled with FR4 and no enclosed air pockets exist), bottom Cu -4.045 ... -4.010, diode sheets z 0.0175, port sheet z -4.010, pin to z {z.get('pin_top')}.
Air box (both): x +-{ab.get('ax')}, y +-{ab.get('ay')} mm; thick z {z.get('air_lo')} ... {z.get('air_hi')}; thin z {zn.get('air_lo')} ... {zn.get('air_hi')}.

## Bodies (every name unique, no two 3D bodies share volume)
METAL_PEC_<net>: one fused solid per electrical net (cell copper + via barrels + via pads + B.Cu traces + connector pad), METAL_PEC_C<site> for each of the 518 static cells, METAL_PEC_GROUND (ground plate with 480 anti-pads and the probe opening), METAL_PEC_PIN. The metal is NOT one single solid because it is not one connected conductor (ground, pin and the nets are electrically separate); hfss_board_setup.py unites all METAL_PEC_* into one object METAL_PEC (multi-lump) when HFSS accepts it, otherwise it keeps them and assigns PEC to all.
SUBSTRATE (RT5880) and SUBSTRATE_FR4 have holes for every barrel and the probe pin (and the ground plate in the FR4); AIR = box minus all solids; diode sheets D<string>_<pos> (240, matching diode_sheets.csv, 1.00 x 0.35 mm between the island edge and the disc edge, touching both); PORT_ANNULUS (r 0.635 ... 2.05).
Gate NO_OVERLAP (pairwise OCC common volume > 1e-9 mm^3 = fail, all bounding-box-overlapping pairs plus AIR against every solid and the AIR volume balance): thick {thick.get('NO_OVERLAP')} ({thick.get('pairs_checked')} pairs, worst common volume {thick.get('worst_common_volume_mm3')}, AIR worst {thick.get('AIR_worst_common_volume_mm3')}, AIR balance {thick.get('AIR_balance_mm3')} mm^3), thin {thin.get('NO_OVERLAP')} ({thin.get('pairs_checked')} pairs, worst {thin.get('worst_common_volume_mm3')}, AIR worst {thin.get('AIR_worst_common_volume_mm3')}, balance {thin.get('AIR_balance_mm3')} mm^3).
Gate SHEETS: thick {thick.get('SHEETS')}, thin {thin.get('SHEETS')} (each sheet touches the island body and the disc body of its cell and has zero overlap area with every solid except AIR; a sheet inside the AIR region of the open slot is normal in HFSS).

## HFSS steps (scripted)
1. Copy the two scripts, the STEP and diode_sheets.csv into one folder and set DIR / THIN at the top of hfss_board_setup.py. In AEDT 2026 R1: Tools > Run Script > hfss_board_setup.py (IronPython, no PyAEDT). It: new Driven Modal design BOARD_N240; imports the STEP (mm); materials RT5880 (2.2 / 0.0009), FR4_assumed (4.4 / 0.02), vacuum, pec; runs the validation check and prints the number of intersect errors (0 expected; if not 0 it subtracts clones of the metal from the dielectrics and of everything from AIR and re-checks); Perfect E / pec for the metal; 240 lumped RLC boundaries D<string>_<pos> (the current line runs across the 1.00 mm gap = across the narrow dimension of the slot, from the anode edge on the island to the cathode edge on the disc; the board-frame direction is (current_dir_x, current_dir_y) in diode_sheets.csv, NOT the unit-cell variable $ph_x; the script checks sheet 1 (perpendicular to the slot axis, 1.000 mm long) and stops if not); lumped port 50 ohm on PORT_ANNULUS (integration line pin -> opening edge; flag USE_WAVE_PORT for a wave port); radiation boundary on the 6 outer AIR faces; Setup1 5.8 GHz, max delta S 0.02, max passes 10 (optional sweep 5.6 / 5.8 / 6.0); infinite sphere theta 0 ... 90, phi 0 ... 360, step 1 deg.
2. State variable: all 240 RLC boundaries use the project variables $Rd, $Ld, $Cd. OFF: R 3876 ohm || C 0.1703 pF (L = 100 uH, about open). ON: R 162.4 ohm || L 0.6113 nH (C = 10 aF, about open). set_state.py (edit STATE) switches them in the open project without rebuilding the geometry or the mesh; re-run the analysis.
3. Run OFF first, then ON. Far field: GainRHCP / GainLHCP at (theta 25, phi 0) and (25, 180); fill results_template.csv.
The scripts follow recorded-script syntax and were NOT tested inside AEDT (no AEDT licence in the build session): the API names (Import options, ValidateDesign, GetMessages, AssignLumpedRLC, AssignLumpedPort, AssignRadiation, InsertSetup, InsertFarFieldSphereSetup) should be compared with a recorded GUI session if a call is rejected.

## Expected RAM / time (ORDER OF MAGNITUDE ONLY, not measured: no HFSS run available)
Thin STEP: about 3 to 8 million tetrahedra after 5 to 8 passes, direct solver 200 to 600 GB RAM (use domain decomposition, DDM, on a cluster or a large-memory node), roughly 1 to 2 h per pass, 10 to 20 h for 10 passes per state.
Thick STEP: 2 to 3 x more tetrahedra (35 um solids and their surface meshes), 0.5 to 1.5 TB, 2 to 3 days: not recommended.
Recommended mesh settings: lambda-based refinement ON (lambda/4 in the dielectrics), max delta S 0.02, percent refinement 20 to 30, at most 10 passes, initial mesh on the thin STEP; length-based operations: top-copper sheets (METAL_PEC_N*, S*, C*) max element length 3 mm, diode sheets 0.2 mm, ground sheet 5 mm, PORT_ANNULUS 0.3 mm, METAL_PEC_PIN 0.5 mm; discrete sweep only after the 5.8 GHz solution has converged; keep the mesh of OFF for ON (set_state.py does not re-mesh) and check the ON delta S once.

## Model prediction (MODEL_ONLY, N 240, calibrated, mesh 96 switched cell)
OFF: beam A {m.get('off_GA', 0):.2f} dBic, beam B {m.get('off_GB', 0):.2f} dBic, balance {m.get('bal', 0):.2f} dB, AR(-3 dB beam) A/B {m.get('off_ARA3', 0):.2f}/{m.get('off_ARB3', 0):.2f} dB, SLL {m.get('sll_off', 0):.1f} dB, eta {m.get('eta_off', 0):.2f}.
ON: beam A {m.get('on_GA', 0):.2f} dBic, B suppression {m.get('supp', 0):.1f} dB, AR {m.get('on_ARA3', 0):.2f} dB, SLL {m.get('sll_on', 0):.1f} dB, eta {m.get('eta_on', 0):.2f}. Gate margin {m.get('margin', 0):.2f} dB: {m.get('gate')}.
Expect HFSS to differ (MoM leakage only, directional beta approximation, bias network not modelled).
"""
    (H / "HFSS_SETUP.md").write_text(md, encoding="utf-8", newline="\n")
    print("docs written; STATUS", STATUS, gates)


if __name__ == "__main__":
    main()
