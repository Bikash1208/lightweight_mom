# HFSS: RT5880 coded cell, 4 diode states, eigenmode 90 deg along x and y (one page). Predictions MODEL_ONLY (`predictions.json`).

**1. Assign material `pec` to `METAL_PEC` (and the ground). STEP imports as vacuum; without this nothing works.** (`METAL_PEC` = patch + ground in one object; if HFSS splits it, assign both parts.)
**Files:** `cell_geometry.step` (one geometry; mm; objects `METAL_PEC`, `SUBSTRATE`, `AIR` = box minus metal, sheets `DX_SHEET`, `DY_SHEET`), `states.csv` (the states differ only in RLC values), `geometry_checks.json`.
Geometry = the MoM model: disc r 3.000, slits 0.4 wide x 5.6 long (**rim 0.2 mm at all four ends**), x diode across the vertical slit at y = +0.8, y diode across the horizontal slit at x = +0.8, sheets 0.4 x 1.0 at z = +0.0175 (mid-plane of the 35 um patch), edges on the slit walls, 0.14 mm apart, nothing else in the slits. Patch 35 um (z 0..0.035), ground 35 um (-1.610..-1.575), substrate z -1.575..0, air 0..30.
**Checks (before export and on the re-imported STEP): all PASS** - bodies/names, no intersecting solids, air = box minus metal, patch = MoM mask, rims at 4 ends, symmetry, sheet edges on patch faces, sheet-sheet 0.14 mm, sheet-other slit wall 0.10 mm, EXACT_DIMS.
**2.** `SUBSTRATE`: RT5880 (eps_r 2.2, tan d 0.0009); `AIR`: vacuum. **3.** **PEC lid:** assign Perfect E to the top face of `AIR` (z = +30). **4.** Project Variables `$ph_x`, `$ph_y` (deg); Lattice Pair 1: x = -3.2 / +3.2 side faces of `SUBSTRATE` and `AIR` (Phase Delay `$ph_x`), Lattice Pair 2: y faces (`$ph_y`).
**5. RLC** on `DX_SHEET` and `DY_SHEET` (Lumped RLC), current line across the gap: DX from (-0.2, 0.8, 0.0175) to (0.2, 0.8, 0.0175); DY from (0.8, -0.2, 0.0175) to (0.8, 0.2, 0.0175). **ON: R = 2.44 ohm only. OFF: C = 28.16 fF only.** Set per `states.csv` (4 states).
**6.** Eigenmode: min frequency 5 GHz, 8 modes, delta f 0.1 %, mesh 1 mm substrate, 0.1 mm slits/sheets. Runs: x: `$ph_x = 90deg`, `$ph_y = 0deg`; y: `$ph_x = 0deg`, `$ph_y = 90deg`. **Clean Up Solutions before each run** (change of RLC values/variables otherwise reuses old results). The bound mode is the TM0 band below the light line 11.711 GHz.

**Sanity check.** Bare substrate (no metal assigned) at 90 deg: you measured **11.07 GHz** (model 11.42, 3 % higher: unexplained, treat as an offset). With `METAL_PEC` assigned, **OFF/OFF must be clearly lower** (model 10.53 GHz, 0.9 GHz below the model bare value): if it is within 0.3 GHz of 11.07, the metal is not assigned.

**Prediction A (as designed, rim 0.2 mm), GHz and mS, lid included, +-0.15 GHz:**
| state | f x-run | f y-run | beta/k0 (x / y) | Y_axis (x / y) |
|---|---|---|---|---|
| xOFF_yOFF | 10.525 | 10.525 | 1.113 / 1.113 | 11.33 / 11.33 |
| xON_yOFF | 10.989 | 10.528 | 1.066 / 1.112 | 7.11 / 11.31 |
| xOFF_yON | 10.528 | 10.989 | 1.112 / 1.066 | 11.31 / 7.11 |
| xON_yON | 11.014 | 11.014 | 1.063 / 1.063 | 6.86 / 6.86 |
Diode on the propagation axis ON: +0.46 GHz; crossed diode: < 0.003 GHz; ON/ON within 0.03 GHz of one ON. Metal assigned but sheets inactive (slits open): 10.66 GHz.
**Prediction B (no rim, four pieces, diode the only bridge; weak):** OFF/OFF 11.33 GHz, on-axis ON 11.20 (-0.13), ON/ON 11.13 (x-run; Y_axis 2.3 / 4.6 / 5.6 mS). A shift of only 0.1-0.2 GHz.
**Back out alpha_beam:** fill `results_template.csv` (state names `OFF_OFF` etc.) and run `python scripts/research/hfss_alpha_backout.py results.csv` (lid-aware transverse resonance; formulas in the script).
