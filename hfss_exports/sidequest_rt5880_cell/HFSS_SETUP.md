# HFSS: RT5880 coded cell, 4 diode states, eigenmode 90 deg along x and y (one page). EXPORT STANDARD (2026-11-05). Predictions MODEL_ONLY (`predictions.json`).

**1. Assign material `pec` to `METAL_PEC` (and the ground). STEP imports as vacuum; without this nothing works.** (`METAL_PEC` = patch + ground; if HFSS splits it, assign both parts.)
**2. Air box from the ground bottom to the lid: `AIR` = one vacuum box z = -1.610 ... +30 mm, ONE face per side wall.** It overlaps `SUBSTRATE` and `METAL_PEC` on purpose (nothing is subtracted); keep `AIR` first in the tree so the later objects take priority (check the mesh view: substrate dielectric, patch conductor). **Periodic pairs go on the four side faces of `AIR` only.** Pairs on single faces of the air above the substrate leave the substrate walls PEC (the earlier "no switching" runs).
**Files:** `cell_geometry.step` (mm; `METAL_PEC`, `SUBSTRATE`, `AIR`, sheets `DX_SHEET`, `DY_SHEET`), `states.csv` (states differ only in RLC values), `geometry_checks.json` (all PASS on the re-imported STEP).
Geometry = the MoM model: disc r 3.000, slits 0.4 x 5.6 (**rim 0.2 mm at all four ends**), x diode across the vertical slit at y = +0.8, y diode across the horizontal slit at x = +0.8, sheets 0.4 x 1.0 at z = +0.0175, edges on the slit walls, 0.14 mm apart. Patch 35 um (z 0..0.035), ground 35 um (-1.610..-1.575), substrate -1.575..0.
**3.** `SUBSTRATE`: RT5880 (eps_r 2.2, tan d 0.0009). **4. PEC lid:** Perfect E on the top face of `AIR` (z = +30). **5.** Project Variables `$ph_x`, `$ph_y` (deg); Lattice Pair 1: x = -3.2 / +3.2 side faces of `AIR` (Phase Delay `$ph_x`); Pair 2: y faces (`$ph_y`).
**6. RLC** (Lumped RLC) on `DX_SHEET`, `DY_SHEET`, current line across the gap: DX (-0.2, 0.8, 0.0175) -> (0.2, 0.8, 0.0175); DY (0.8, -0.2, 0.0175) -> (0.8, 0.2, 0.0175). **ON: R = 2.44 ohm only. OFF: C = 28.16 fF only** (`states.csv`).
**7.** Eigenmode: min frequency 5 GHz, 8 modes, delta f 0.1 %, mesh 1 mm substrate, 0.1 mm slits/sheets. x run: `$ph_x = 90deg`, `$ph_y = 0deg`; y run: swapped. **Clean Up Solutions before each run.** Bound mode = TM0 below the light line 11.711 GHz.

**Sanity.** Bare substrate (no metal): model 11.42 GHz (the earlier 11.07 was measured with PEC substrate walls: re-measure). With `METAL_PEC` assigned, OFF/OFF must be clearly lower (HFSS 10.670).

**HFSS result (2026-11-05, x run 90 deg, GHz) against the model:**
| state (DX/DY) | HFSS | model | model, OFF = open slit |
|---|---|---|---|
| C/C | **10.670** | 10.525 | 10.662 |
| R/C | **10.665** | 10.990 | |
| C/R | **10.690** | 10.528 | |
| R/R | **11.006** | 11.014 | 11.014 |
Both diodes together switch (+0.336 GHz); each alone does not. C/C equals the model with no diode capacitance; R/R equals the model (within 0.1 %). The MoM cannot reproduce the single-diode behaviour (`results/sidequest_hybrid/CELL_INDEPENDENCE.md`). **HFSS diagnostics:** (1) in R/C list every mode 10.4-11.2 GHz with its Ez fraction: is the TM0 band tracked in each state? (2) R/C with R = 0.01 ohm; (3) R/C with DY removed; (4) C/C with C = 282 fF (must fall clearly below 10.5).
**Back out alpha_beam** (the both-diode scalar switch): `python scripts/research/hfss_alpha_backout.py results.csv` -> dY -3.13 mS, alpha_beam 0.844 Np/m, beta/k0 1.073.
