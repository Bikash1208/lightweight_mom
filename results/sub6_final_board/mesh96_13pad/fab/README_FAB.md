# Fabrication README: SUB6 final board, 13-pad bias layout (FINAL_CANDIDATE_13PAD / MODEL_ONLY / NOT_HFSS_VALIDATED)
RF is unchanged from the 24-pad board (cells, discs, islands, diode sites, vias, anti-pads, ground, probe): RF_UNCHANGED. Only the B.Cu bias routing and the edge connector changed.
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

## Connector (13 pads, 1.0 x 1.0 mm, pitch 2.54, B.Cu, near the left flat edge, x = -179.705 mm; RETURN in the middle)
| pad (from top) | net | string | x mm | y mm | V at 10 mA | I mA |
|---|---|---|---|---|---|---|
| 1 | S11 | 11 | -179.705 | 15.24 | 19.0 | 10 |
| 2 | S12 | 12 | -179.705 | 12.7 | 19.0 | 10 |
| 3 | S10 | 10 | -179.705 | 10.16 | 19.0 | 10 |
| 4 | S09 | 9 | -179.705 | 7.62 | 19.0 | 10 |
| 5 | S07 | 7 | -179.705 | 5.08 | 19.0 | 10 |
| 6 | S08 | 8 | -179.705 | 2.54 | 19.0 | 10 |
| 7 | RETURN | all 12 (cathode ends) | -179.705 | 0.0 | 0 (return) | 120 |
| 8 | S06 | 6 | -179.705 | -2.54 | 19.0 | 10 |
| 9 | S01 | 1 | -179.705 | -5.08 | 19.0 | 10 |
| 10 | S05 | 5 | -179.705 | -7.62 | 19.0 | 10 |
| 11 | S04 | 4 | -179.705 | -10.16 | 19.0 | 10 |
| 12 | S03 | 3 | -179.705 | -12.7 | 19.0 | 10 |
| 13 | S02 | 2 | -179.705 | -15.24 | 19.0 | 10 |
Each drive net S<string> feeds the ANODE end of its string (20 diodes in series, about 19.0 V at 10 mA); the 12 cathode ends (last diode of each string) are joined on B.Cu into the single net RETURN (12 x 10 mA = 120 mA).
Net naming: S<string> follows the string id of diode_sheets.csv (D<string>_<pos>); the pad order along the connector is a planar fan order (not 1..12): use the table.

## Files
fab/gerber_kicad/ (kicad-cli export of the routed board): sub6_board_routed-F_Cu.gtl, -In1_Cu.g1, -B_Cu.gbl, -F_Paste.gtp, -Edge_Cuts.gm1, sub6_board_routed.drl (0.30 plated vias + 1.27 probe hole, one Excellon file, MixedPlating).
fab/gerber/ (script-generated RF layers, byte-identical to the 24-pad files of 526f98de): sub6_final-F_Cu.gtl, sub6_final-In1_Cu.g2, sub6_final-PTH.drl, sub6_final-Edge_Cuts.gm1. Use these for F_Cu / In1_Cu / drill if byte-identical RF is wanted: the KiCad F_Cu / In1_Cu are the same geometry with circles written as 2 degree polygons (max deviation 0.0008 mm).
fab/assembly/: F_Paste (480 lands), F_Assy, diode_placement.csv, DXF/PDF assembly drawing, README_ASSEMBLY.md (polarity: anode on the island, cathode on the disc, current along the string).
fab/sub6_final_board.dxf, cell_table.csv, string_table.csv, connector_pinout.csv, drill.csv, B_Cu_full.png, B_Cu_connector_zoom.png; kicad/sub6_board_routed.kicad_pcb (+ .kicad_pro, .kicad_dru) editable in KiCad >= 10.

## Rules and results
Min clearances (copper edge to copper edge, mm): trace-trace 0.167, trace to foreign anti-pad edge 0.1885 (required 0.172), trace to foreign via pad 0.3885, trace to connector pad 0.1669, trace to board outline 2.0211, trace to probe opening 3.309 (all >= 0.15).
Widths: signal 0.25 mm (short neck-downs to 0.186 mm near pads where the clearance does not allow 0.25), RETURN {'0.374': 4, '0.5': 81} segments of 0.5 / 0.374 mm (min 0.374 mm, {'0.374': 4.4, '0.5': 1055.8} mm of length). 35 um Cu carries about 0.8 A at 0.15 mm width (IPC-2221, 20 K rise), so 120 mA needs no minimum width; the 0.5 mm RETURN is a margin.

## SMA connector (flat on the bottom, no legs)
Body 6.35 x 6.35 mm soldered to a B.Cu ground pad 7.35 x 7.35 mm (body + 0.5 mm per side, opening r 2.05 around the pin); bias keep-out 10.35 x 10.35 mm; 12 blind GND vias (B.Cu to In1.Cu, FR4 only, R 2.7 mm, pitch 1.41 mm) in fab/gerber_kicad/sub6_board_routed-in1-back.drl (T1 0.30, layer pair In1-B.Cu) and the main .drl; recommendation and the through-via alternative: see HFSS_SETUP.md. Nearest bias trace edge 5.359 mm from the centre (keep-out half size 5.175).

## Points to check with the manufacturer (not resolved here)
1. Via aspect ratio: 0.30 mm drill through about 4.08 mm = 13.6:1, above the usual 8 to 10:1 plating limit; mixed RT5880 / FR4 stack with a bonded ground layer. Ask for feasibility or a larger drill.
2. The RETURN bus reaches the middle pad from the board-edge side (it passes the pads S12..S08 on the edge side): a single-layer fan cannot reach a middle RETURN pad otherwise. With RETURN at the first or last pad the bus would enter directly.
3. BAR50-type diode assumed for the model; check the SC-79 land pattern of the real part (0.35 x 0.35 at 1.35 pitch used). No solder mask / silkscreen / paste layers beyond F_Paste are generated; FR4 properties assumed.
