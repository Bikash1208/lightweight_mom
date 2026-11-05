# Coded hologram, K = 2 with a different polarisation per state (round 8) — MODEL_ONLY

> **PENDING CELL VALIDATION (2026-11-05):** HFSS (air box, pec metal) shows the crossed two-diode rim-disc cell switches ONLY with both diodes ON (C/C 10.670, R/C 10.665, C/R 10.690, R/R 11.006 GHz at 90 deg); the model's per-axis independence is falsified. Every number in this document that depends on independent per-axis switching of this cell (alpha per axis, polarisation switching, board gain, diode counts) is on hold. See `results/sidequest_hybrid/CELL_INDEPENDENCE.md`.


State A **RHCP at (25 deg, 0 deg)**, state B **LHCP at (25 deg, 180 deg)**; 150 mm hexagon, 10.54 GHz, Phase 26 RT5880 stack, centre probe (launch 0.627), 1430 cells, one MA4AGP907 per switched cell.
Code: `scripts/research/sidequest_coded_pol_cell.py`, `sidequest_coded_pol.py`; data `coded_pol_cell.json`, `coded_pol_results.json`; tests `tests/test_sidequest_coded_pol.py`.

## 1. The rotated-load bug (own commit, `src/loads.py`, `src/ris/diodes*.py`, `tests/test_loads_rotated_diode.py`)
* **Cause.** `diodes_measured.load` / `diodes.diode_load` always set `transverse_conducts=False`. For a rotated axis this adds `OPEN_OHM (v.J)^2` over the footprint, i.e. the pointwise condition `v.J = 0` on x/y rooftops, which for the 3-4 mesh cells of a diode footprint admits no current at all. The footprint is then an open circuit in **every** diode state (ON == OFF to 1e-9). Not a quadrature or sign error: the uniform-current test of `test_loads_rotated.py` passes because uniform J is the one admissible pattern. (Known to R3: `results/phase20_r3/CELLS.md`, "a rotated diode region must conduct transversely".)
* **Fix (additive).** `loads.default_transverse_conducts(direction)`: False for x/y/0/90 (unchanged, bit-identical `load_matrix`), True for any other angle; `diode_load`, `diodes_measured.load` and `_chain_load` default to it (`transverse_conducts=None`), an explicit value is still honoured; `check_load_geometry` warns for a rotated blocked load. No solver kernel changed.
* **Checks (10 tests, 50 s).** ON != OFF at 30/45/60 deg (contrast 0.7 mS on the test cell; the old setting gives ON == OFF); rotation covariance of the sheet tensor (u-component and ON-OFF contrast within 10 % over 0/45/90/135 deg, principal axis within 1 deg, 0 vs 90 deg an exact axis swap); unloaded rotation invariance (90 deg exact, 45 deg within 20 % of the staircased patch); reciprocity (Y_xy = Y_yx to 1e-9) and mirror symmetry (30 vs 150 deg); axis-aligned `load_matrix` equal; near-axis rotated path (0.01 deg) within 5 % of the exact 0 deg value.

## 2. Cells
A rotating cell must fit the 6.4 mm lattice at every angle, so the 6.0 x 4.8 mm x-only pair of round 7 (|dY| 7.6 mS, alpha 4.0 Np/m) cannot be turned. Rotation-invariant envelope: **split disc, diameter 6.0 mm**, gap 0.4 mm, footprint 1.2 mm (mesh 32, mean over axes 0/22.5/45 deg): OFF 2.22j mS, ON 0.61 + 7.15j mS along the axis, 5.07j mS across (static), |dY| = 4.97 mS, mean sheet near isotropic (beta/k0 1.039). First-order **alpha_beam = 0.67 Np/m** (6x lower), alpha_abs 0.22 Np/m. Footprint 2.4 / 3.6 mm and gap 0.2 / 0.6 mm are worse (alpha 0.25 / 0.18 / 0.00 / 0.20). The weak contrast is the main cost of rotated cells.

## 3. Model
Cell current `J = (pi/2) A(rho) (u.rho_hat) u e^{-j beta rho}` (|J|^2 ~ cos^2 delta, delta = axis minus azimuth), far field with the scalar element factor of `sidequest_arch_study`, co-pol = conj(Jones).E. Per cell and state the diode sign s = +-1; per cell the axis chi (8 angles, 22.5 deg steps) is chosen jointly with (s_A, s_B) to maximise s_A Re(e^{j eps_A} c_A) + s_B Re(e^{j eps_B} c_B) minus lambda for s_A != s_B (a diode costs). eps_B and lambda are scanned; gain checked on the full 1 deg x 2 deg pattern. Mean cos^2 delta of the chosen axes 0.84 (the optimiser puts each axis near the local ray). R2 (two diodes: x and y elements, each independent) assumes the **same dY per axis** as the disc (not computed for a crossed cell: optimistic).

## 4. Results (full pattern; baseline = lambda 0)
| route / lambda | static | diodes | lines (inv / no inv) | gain A / B [dBi] | SLL A / B [dB] | AR A / B [dB] | XPD [dB] |
|---|---|---|---|---|---|---|---|
| R1 rotated, 0 (ideal bits) | 52.7 % | 676 | 1 / 2 | 16.6 / 16.6 | -17.4 / -17.4 | 1.9 / 1.9 | 19.3 |
| **R1, 0.5 (recommended)** | **63.8 %** | **518** | **1 / 2** | **16.4 / 16.4** | -13.8 / -13.3 | **1.6 / 1.7** | 20.8 / 20.2 |
| R1, 0.8 | 71.3 % | 410 | 1 / 2 | 16.1 / 16.1 | -9.9 / -9.8 | 1.6 / 1.4 | 21 |
| R1, 1.2 | 81.4 % | 266 | 1 / 2 | 15.5 / 15.6 | -5.8 / -6.0 | 1.6 / 1.5 | 21 |
| R2 two diodes, 0 | 51.3 % | 1393 | 1 / 2 | 17.3 / 17.3 | -14.5 / -14.1 | 0.9 / 0.8 | 25 / 27 |
| **R2, 0.5 (recommended)** | **72.0 %** | **801** | **1 / 2** | **17.0 / 16.9** | -14.1 / -13.3 | **1.2 / 0.9** | 23.5 / 26.0 |
| R2, 0.8 | 78.1 % | 626 | 1 / 2 | 16.6 / 16.7 | -10.1 / -10.8 | 0.9 / 0.9 | 26 |
Pointing is exact (25.0, 0) and (25.0, 180) in every row. The ideal-pol reference (continuous amplitude, perfect CP, same aperture) is 18.0 (R1) / 18.6 (R2) dBi. Both states use the same word set: static (+,+), (-,-) and the complementary pair (+,-), (-,+), so **1 control line plus an inverter** (2 lines without) for either route. The 1 dB budget is not the limit; the sidelobes are (R1 above ~64 % static, R2 above ~72 %).

## 5. Reading
* **Rotated 1-diode cells work** (needed the fix): per-state RHCP/LHCP, AR 1.6-1.7 dB, XPD ~20 dB, 518 diodes, 1 line + inverter. Opposite handedness needs the PB-type phase of one axis to be +2chi for A and -2chi for B, so a fixed axis serves both states only where both codes agree; the optimiser finds it by choosing chi per cell (baseline design: the 8 axes are used about equally, 144-205 cells each) and the cost is the 1-bit loss, not a lost beam.
* **Absolute gain 16.4 dBi vs 18.5 dBi for the LP K = 2 of round 7**: the disc cell has 6x lower alpha, which flattens the aperture taper (higher directivity) but converts less power (eta = launch x converted = 0.082 for R1, 0.096 for R2); this is a cell limit, not a hologram limit. R2 (+0.6 dBi, better AR/XPD) uses every cell at full depth but doubles the diodes (801 vs 518) and relies on an unverified 2-diode cell.
* **Rotated-cell staircase**: 22.5 and 45 deg patches are staircased on the 0.2 mm mesh; u-component moves +-8 % between 0/22.5/45 deg (mesh 32).

## Risks
Scalar element factor and first-order leakage for a rotated tensor (no HFSS); alpha of the disc varies 0.54-0.67 Np/m with footprint staircase and mesh; the R2 cell is assumed, not solved.

## 6. Solved crossed two-diode cell (round 9) — replaces the assumed R2 numbers above
Cell (`sidequest_coded_cross_cell.py`; data `coded_cross_scan.json`, `coded_cross_rim_scan.json`): the 6.0 mm disc cut by a vertical and a horizontal 0.4 mm slit, x diode (0.4 x 1.2 mm, along x) bridging the vertical slit at y = +1.6 mm, y diode bridging the horizontal slit at x = +1.6 mm; measured MA4AGP907, MoM mesh 64. Two slit lengths:
* **Through slits** (four islands, one floats): weak, rejected. Y_xx OFF/OFF 2.15j; ON/OFF 4.05j; dY 1.9j (y OFF) / 2.9j mS (y ON); alpha_beam 0.070 Np/m, beta/k0 1.030; ON/ON has Y_xy = -1.92j mS (not additive). Four diodes (two per axis) only reach 0.11 Np/m.
* **Rim slits (recommended cell)**: slits 5.6 mm long, leaving a 0.2 mm metal rim at each end, so the islands stay connected and the diode in ON state *removes* slit capacitance. Four states, Y in mS [xx, xy, yy]:
| x/y diode | Y_xx | Y_xy | Y_yy |
|---|---|---|---|
| OFF/OFF | 0.01 + 9.48j | 0.03j | 0.01 + 9.48j |
| ON/OFF | 0.11 + 6.25j | -0.11j | 0.01 + 9.48j |
| OFF/ON | 0.01 + 9.48j | -0.11j | 0.11 + 6.25j |
| ON/ON | 0.12 + 6.16j | 0.43j | 0.12 + 6.16j |
  dY per axis (ON - OFF) **-3.28j + 0.10 mS**, the other axis moves by < 0.1 mS (ON/ON deviates from additive by 0.09j mS in Y_xx and Y_xy = 0.43j mS). **alpha_beam 0.92 Np/m** (mesh 32: 1.15), alpha_abs 0.16 Np/m, **beta/k0 1.065** (mean sheet 7.84j mS = the Phase 26 design loading). The sign of dY only relabels which diode state is logic 1.
  Risks: the 0.2 mm rim is one mesh cell at mesh 32 (alpha +25 % vs mesh 64) and a fabrication limit; 2 diodes per cell.
K = 2 RHCP/LHCP with the solved cells (same optimiser; `sidequest_coded_pol_cross.py rim|through`, `coded_pol_cross_{rim,through}_results.json`):
| cell / lambda | static | diodes | lines (inv / no inv) | gain A / B [dBi] | SLL [dB] | AR [dB] | XPD [dB] |
|---|---|---|---|---|---|---|---|
| rim, 0 | 50.9 % | 1404 | 1 / 2 | 18.5 / 18.5 | -14.7 | 0.8 | 26.6 |
| **rim, 0.5 (recommended)** | **73.2 %** | **766** | **1 / 2** | **18.2 / 18.2** | -13.3 | 0.8 | 27.0 |
| rim, 0.8 | 79.8 % | 578 | 1 / 2 | 17.8 / 17.8 | -9.7 | 0.9 | 26 |
| through, 0.5 | 72.3 % | 792 | 1 / 2 | 7.6 / 7.6 | -13.8 | 0.9 | 25.5 |
Ideal-pol reference 19.9 dBi (rim), eta = 0.129. Pointing exact (25, 0)/(25, 180).
