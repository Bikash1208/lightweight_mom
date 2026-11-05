# HFSS: RT5880 cell, 4 diodes (both halves of each slit bridged), independent-axes check (one page). EXPORT STANDARD. Predictions MODEL_ONLY (`predictions.json`).

**1. Assign material `pec` to `METAL_PEC` (and the ground). STEP imports as vacuum; without this nothing works.**
**2. `AIR` = one vacuum box from the ground bottom to the lid (z -1.610 ... +30), one face per side wall; it overlaps the other solids on purpose. Lattice pairs on the four side faces of `AIR` only** (`$ph_x`, `$ph_y`), keep `AIR` first in the tree. **3.** `SUBSTRATE`: RT5880 (2.2, tan d 0.0009). **4.** PEC lid: Perfect E on the top face of `AIR`.
**Files:** `cell_geometry_4diode.step` (`METAL_PEC`, `SUBSTRATE`, `AIR`, sheets `DX_TOP`, `DX_BOTTOM`, `DY_RIGHT`, `DY_LEFT`), `states.csv`, `geometry_checks.json` (all PASS on the re-imported STEP). Geometry = the 2-diode cell of `../sidequest_rt5880_cell/` (disc r 3.000, slits 0.4 x 5.6, rim 0.2 at all four ends) plus the mirrored sheets: 0.4 x 1.0 at (0, +-0.8) and (+-0.8, 0), z +0.0175.
**5. RLC** (one boundary per sheet), current line across each gap: DX_TOP (-0.2, 0.8, 0.0175)->(0.2, 0.8, 0.0175); DX_BOTTOM (-0.2, -0.8, 0.0175)->(0.2, -0.8, 0.0175); DY_RIGHT (0.8, -0.2, 0.0175)->(0.8, 0.2, 0.0175); DY_LEFT (-0.8, -0.2, 0.0175)->(-0.8, 0.2, 0.0175).
**x diodes (DX_TOP + DX_BOTTOM) always carry the same state, y diodes (DY_RIGHT + DY_LEFT) the same state: ON = R 2.44 ohm only, OFF = C 28.16 fF only** (`states.csv`).
**6.** Eigenmode: min frequency 5 GHz, 8 modes, delta f 0.1 %, mesh 1 mm substrate, 0.1 mm slits/sheets. x run: `$ph_x = 90deg`, `$ph_y = 0deg`; y run: swapped. **Clean Up Solutions before each run.** Track the TM0 band below the light line 11.711 GHz (list the modes 10.2-11.6 GHz with their Ez fraction in every state).
**Why this check:** the 2-diode cell switched only with both diodes ON (HFSS C/C 10.670, R/C 10.665, C/R 10.690, R/R 11.006 GHz) while the MoM predicts independent switching. Here each slit is bridged on both halves. **Decisive runs: x ON only (R/C) and y ON only (C/R).**

**Prediction, x-run / y-run, GHz (MoM mesh 64, +-0.15; OFF modelled as open slit (HFSS C/C ~ open) | as 28.16 fF):**
| state (x/y diodes) | OFF = open | OFF = 28.16 fF |
|---|---|---|
| xOFF_yOFF | 10.662 / 10.662 | 10.370 / 10.370 |
| xON_yOFF | 11.110 / 10.648 | 11.105 / 10.368 |
| xOFF_yON | 10.648 / 11.110 | 10.368 / 11.105 |
| xON_yON | 11.105 / 11.105 | 11.105 / 11.105 |
Independent in the MoM: x diodes alone move the x run by +0.45 GHz (dY_xx -4.4 mS at 10.8 GHz), the y diodes shift it by < 0.015 GHz (cross coupling 0.1 mS); alpha_beam 1.6 Np/m (OFF = open) / 3.7 (OFF = 28 fF). If HFSS again shows no single-diode switching, the cause is in the HFSS lumped boundaries/modes, not the half-bridging.
