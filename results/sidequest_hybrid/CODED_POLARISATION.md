# Coded hologram, K = 2 with a different polarisation per state (round 8) — MODEL_ONLY

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
Cell (`sidequest_coded_cross_cell.py`, data `coded_cross_scan.json`): the 6.0 mm disc cut by a vertical and a horizontal 0.4 mm slit (through the disc, four islands), x diode (0.4 x 1.2 mm, along x) bridging the vertical slit at y = +1.6 mm, y diode bridging the horizontal slit at x = +1.6 mm; one island floats. Measured MA4AGP907, MoM mesh 64 (mesh 32 within 5 %). Four states (Y in mS, [xx, xy; yx, yy]):
| x/y diode | Y_xx | Y_xy | Y_yy |
|---|---|---|---|
| OFF/OFF | 2.15j | -0.02j | 2.15j |
| ON/OFF | 0.13+4.05j | -0.18j | 2.16j |
| OFF/ON | 2.16j | -0.18j | 0.13+4.05j |
| ON/ON | 0.44+5.07j | -1.92j | 0.44+5.07j |
* dY per axis: x diode alone 1.90j mS (y OFF), 2.91j mS (y ON); mean **2.40j + 0.28 mS**; same for y. The two diodes are **not additive** (ON/ON has Y_xy = -1.92j mS and a larger dY): not in the model below.
* **alpha_beam 0.070 Np/m** (mesh 32: 0.073), alpha_abs 0.155 Np/m, **beta/k0 1.030** (mean sheet 3.36j mS). The disc with one rotated diode has dY 4.97 mS and 0.67 Np/m: the crossed cell is **10x weaker**, because one diode bridges only 1.2 mm of the 6 mm island row and the other island floats.
* Variants (mesh 64): four diodes (two per axis, driven in pairs): dY 2.72 mS, alpha 0.110 Np/m; half-length slits (rim connects the islands): dY of the wrong sign or ~0 (alpha <= 0.09 / 1.1 at a rim-limited sheet shift): rejected.
K = 2 RHCP/LHCP with this cell (same optimiser, `sidequest_coded_pol_cross.py`, `coded_pol_cross_results.json`):
| lambda | static | diodes | lines (inv / no inv) | gain A / B [dBi] | SLL [dB] | AR [dB] | XPD [dB] |
|---|---|---|---|---|---|---|---|
| 0 | 51.7 % | 1380 | 1 / 2 | 7.9 / 8.0 | -14.9 / -14.4 | 1.0 / 0.8 | 25 / 26 |
| **0.5** | **72.3 %** | **792** | **1 / 2** | **7.6 / 7.6** | -13.8 / -13.8 | 0.9 / 0.9 | 25.5 |
| 0.8 | 78.1 % | 626 | 1 / 2 | 7.3 / 7.3 | -10.6 / -10.4 | 1.0 / 0.9 | 25 |
Ideal-pol reference 9.3 dBi; eta = 0.011. Polarisation purity is excellent, the gain is not: **7.6 dBi vs 16.4 dBi for the rotated disc**.
