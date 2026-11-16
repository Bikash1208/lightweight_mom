# BOARD_FINAL_SPREAD_ELS.aedt read-back (UNSOLVED, clean)
Project reopened after the final save (AEDT 2025.1, PyAEDT 1.9.0, validate_simple only). File saved 16 Nov 10:57; reopen and read-back did not change it.

- Designs: OFF_ELS10, ON_ELS, OFF_ELS10_ABS, ON_ELS_ABS (no extra designs).
- Layout: R = 198 mm cell area, N = 200 switched cells, uniform (area-uniform) placement, cap 0.3 Np/m, 856 sites, 656 static, rot180 frame. Board radius 228 mm = 198 + 30 mm absorber ring (D = 456 mm).
- Absorber (ABS designs only): inner 198, outer 228 mm, z 0 to 30.035 mm, bbox [-228.0, -228.0, 0.0, 228.0, 228.0, 30.035], er = ur = 1.5 - j1.0 (idealised); material read-back eps 1.5, tan d 0.666666666666667, mu 1.5, tan dm 0.666666666666667.
- Sheets (read back): OFF_ELS10 Rs 21.6 Xs -297.5; ON_ELS Rs 2.8 Xs 14.6; 200 impedance sheets per design.
- Setup1: 5.8GHz, MaximumPasses 20, MaxDeltaS 0.02, MinimumPasses 6, basis order 1, solver Direct Solver; interpolating sweep 5.0-6.5 GHz, 0.05 GHz step; spheres Sphere1 (full, 2 deg) and Lower (theta 90-180).
- Mesh operations: none in any design (automatic mesh only). No mesh or solution data saved (results folder holds geometry cache only). validate_simple = 1 for all four designs.
- Port: Port1 lumped, 50 ohm, renormalised; radiation boundary on AIR.
- Solver choice: direct. No mesh estimate was made; the coherence-layout mesh (1.41 M tets, about 60 GB direct at pass 1, growing to about 110 GB) is the reference, and this board has 856 sites instead of 1072, so the direct solver is expected to stay under 120 GB. If the first pass exceeds that, switch to the iterative solver.
- Exporter checks: 11 of 11 PASS (`checks.md`).
- Model prediction for this layout (MODEL_ONLY, no HFSS offsets applied): balance A-B +1.59 dB in the two-beam state; one-beam SLL about -13.9 dB (about 1.1 dB short of -15); two-beam SLL about -16.6 dB; converted power 0.204 / 0.081. Cell values used: alpha 2.16 Np/m, r = 0.262 (from the Ls 5.5 cell gate: Delta_OFF +0.0273, Delta_ON +0.1040).
