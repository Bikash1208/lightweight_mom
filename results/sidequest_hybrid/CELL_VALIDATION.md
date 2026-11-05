# Cell validation after the first HFSS runs (2026-11-05)
**Result:** HFSS (90 deg, x run) gave OFF 11.059 / ON 11.055 GHz (y sheet removed: 11.086 / 11.081); the model predicted 10.53 / 11.00. **Root cause (user): the STEP carries no material, the patch imported as vacuum**, so HFSS solved a bare grounded substrate (user-measured bare value 11.07 GHz; the identical OFF/ON numbers are what an absent metal gives). The first package also had geometry errors; all are fixed in the new package `hfss_exports/sidequest_rt5880_cell/`.

## Geometry errors found and fixed (automatic checks `scripts/research/check_step_cell.py`, `tests/test_hfss_cell_package.py`)
| item | old package | new package |
|---|---|---|
| sheet-to-sheet contact | footprint 0.4 x 1.2 at 0.8 mm from the centre: inner edge on the slit-crossing edge, **the two sheets touched at the corner (0.2, 0.2)** (distance 0.000 mm) | footprint 0.4 x **1.0** (centre kept at 0.8 mm): sheets 0.141 mm apart, 0.100 mm from the other slit's wall |
| rims | STEP re-opened with OCC: rim present at all four ends, 4-fold/mirror symmetric (0 asymmetric points of 289); the missing top/right rim reported from HFSS could not be reproduced in the STEP, and the new check would flag it | rims present at 4 ends, symmetry 0/289, MoM-mask mismatch 0 of 9159 points |
| bodies | 4 STEP files, `Patch`/`Ground`/`Air` separate, no material hint | **one** `cell_geometry.step`: `METAL_PEC` (patch + ground), `SUBSTRATE`, `AIR` (box minus metal, common volume with the metal 0), `DX_SHEET`, `DY_SHEET`; states in `states.csv` |
| material | not stated | **step 1 of the sheet: assign `pec` to `METAL_PEC`** (STEP imports as vacuum) |
Checks run on the shapes before export and on the re-imported STEP (all PASS): bodies_and_names, no_intersections (pairwise common volume 0), air_is_box_minus_metal, patch_vs_MoM, rims_at_4_ends, symmetry, sheet_edge_coincidence (long edges and end points on the patch faces, centre in the gap, z = 0.0175, 0.4 x 1.0), sheet_to_sheet_clearance (0.141 >= 0.1), sheet_to_other_slit_clearance (0.100 >= 0.1), EXACT_DIMS (0 inexact literals; snap <= 2e-16 mm). The checks catch the old errors (unit tests: footprint 1.2 -> sheets touch; slits through the disc -> rim missing; shifted slit -> rim/symmetry fail).

## Model changes
* The HFSS model has a PEC lid 30 mm above the patch plane. The transverse resonance now uses the lid-aware air admittance (-j omega eps0 / (kz tan(kz L))); effect on the bare slab at 90 deg: 11.46 -> 11.42 GHz, on the cell predictions < 0.02 GHz.
* **Bare substrate: model 11.42 GHz, HFSS (user) 11.07 GHz: 3 % lower, unexplained** (not the lid; candidate: laminate/ground/boundary details of the HFSS model). Treat as an offset when comparing absolute frequencies; differences between states are the discriminating quantity.
* The wider footprint change (1.2 -> 1.0 mm) leaves the prediction essentially unchanged: OFF/OFF 10.53, on-axis ON 10.99, ON/ON 11.01 GHz (x run), crossed diode < 0.003 GHz.

## Predictions for HFSS (A as designed, B without rim; `predictions.json`; GHz, +-0.15)
| | OFF/OFF | on-axis ON | ON/ON | sheets inactive, metal assigned | bare substrate |
|---|---|---|---|---|---|
| A (rim 0.2 mm) | 10.53 | 10.99 | 11.01 | 10.66 (slits open) | 11.42 (HFSS 11.07) |
| B (no rim, 4 pieces, diode the only bridge) | 11.33 | 11.20 | 11.13 | | |
B is a weak switch (ΔY 2.3 -> 4.6 mS on-axis, shift only -0.13 GHz, alpha ~0.07 Np/m, earlier through-slit result). The OFF/OFF-below-bare test discriminates "metal assigned": A predicts 0.9 GHz below the model bare value.

## Rim convergence (priority 2) and the no-rim cell (priority 3)
No collapse with mesh: Delta Y (on-axis, rim 0.2 mm, pos 0.8) = -5.45 / -4.84 / -4.53 mS at mesh 32 / 64 / 96 (1 / 2 / 3 cells across the rim), still converging toward about -4.2 mS (`CODED_CELL_SEARCH.md`). Mesh 128 (4 cells across 0.2 mm) does not fit in memory (32 k unknowns, ~17 GB dense matrix); not run. The vacuum-patch diagnosis explains the HFSS result without invoking the rim; the rim sensitivity stays open until a cell is HFSS-confirmed. B predicted above.
