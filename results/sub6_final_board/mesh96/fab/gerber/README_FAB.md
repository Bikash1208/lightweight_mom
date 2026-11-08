# Fabrication README: SUB6 final board (hexagon, N 240, mesh-96 switched-cell design)
Status: MODEL_ONLY / NOT_HFSS_VALIDATED; design status in ../layout_summary.json (FINAL_CANDIDATE).  Files are script-generated (scripts/research/sub6_final_gerber.py); Gerber RS-274X, format 4.6 mm, Excellon METRIC.

## Stack-up (top to bottom)
| layer | material | thickness |
|---|---|---|
| L1 top copper | 35 um Cu (cells: discs, slots, islands) | 0.035 |
| core | Rogers RT5880 (eps_r 2.2, tan d 0.0009) | 3.175 |
| L2 ground | 35 um Cu, anti-pads 1.0 mm at every via, 4.1 mm opening at the probe | 0.035 |
| dielectric | FR4 (eps_r 4.4, tan d 0.02 ASSUMED) | 0.8 |
| L3 bottom copper | 35 um Cu (bias traces, via pads, connector pads) | 0.035 |
Total about 4.08 mm.  Board outline: regular hexagon, circumradius 211.547 mm (flat-to-flat 366.410 mm), vertices up (+y), origin = probe centre.

## Files
sub6_final-F_Cu.gtl (L1), sub6_final-In1_Cu.g2 (L2 ground), sub6_final-B_Cu.gbl (L3), sub6_final-Edge_Cuts.gm1 (outline), sub6_final-PTH.drl (Excellon); ../cell_table.csv, ../string_table.csv, ../drill.csv, ../sub6_final_board.dxf.
Polarity: dark = copper; the top layer uses clear polarity for the slots, so apply the file in order.

## Drill and rules
Vias: 480 plated holes, drill 0.30 mm, pad 0.6 mm (top land = the disc or island itself), anti-pad 1.0 mm in L2.  Probe hole 1.27 mm (pin r 0.635) at the origin, 4.1 mm opening in L2, 2.0 mm pad on L3.
Min trace 0.15 mm, min copper-to-copper gap on L3 0.274 mm (all trace/pad/anti-pad/outline clearances >= 0.15 mm, see ../layout_summary.json).  Connector: 24 pads 1.0 x 1.0 mm, pitch 2.54 mm (one column near the left flat edge, 3.5 mm from the pad centres to the outline).

## Points to check with the manufacturer (not resolved here)
1. Via aspect ratio: 0.3 mm drill through about 4.08 mm = 13.6 : 1, above the usual 8-10 : 1 plating limit; mixed RT5880 / FR4 stack with a bonded ground layer.  Ask for feasibility or a larger drill.
2. The BAR50 diode footprint is NOT defined: the Gerbers leave the 1.00 mm gap between island and rim open (the HFSS model uses a 1.00 x 0.35 mm RLC sheet).  Pads for the real diode package must be added; diode orientation: anode on the island via, cathode on the disc via of the same cell (series string island(i) - diode - disc(i) - trace - island(i+1)).
3. No solder mask, silkscreen or paste layers are generated; no copper thickness or finish is specified beyond 35 um.
4. FR4 properties are assumed typical values.  The bias strings (12 x 20 diodes, 19.0 V at 10 mA each) end on the 24 connector pads (see ../string_table.csv for the pad of each string start/end).
