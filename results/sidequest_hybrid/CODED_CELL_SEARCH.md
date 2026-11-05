# Cell search for a stronger leakage cell (round 12): RO3010 primary vs RT5880 secondary — MODEL_ONLY

> **PENDING CELL VALIDATION (2026-11-05):** HFSS (air box, pec metal) shows the crossed two-diode rim-disc cell switches ONLY with both diodes ON (C/C 10.670, R/C 10.665, C/R 10.690, R/R 11.006 GHz at 90 deg); the model's per-axis independence is falsified. Every number in this document that depends on independent per-axis switching of this cell (alpha per axis, polarisation switching, board gain, diode counts) is on hold. See `results/sidequest_hybrid/CELL_INDEPENDENCE.md`.

Code: `scripts/research/sidequest_cell_search.py` (cell builders + loaded Bloch tensor), `..._full.py` (4 states), `..._sens.py` (eps_r, mesh), `sidequest_search_board.py`, `..._board_sens.py` (K = 2 boards, same engine as round 11 with a 3rd tensor element and a background loss). Data `cell_search_*.json`, `cell_search_board_*.json`. Tests `tests/test_sidequest_cell_search.py`.

## 0. Correction of the earlier RT5880 numbers (read first)
Rounds 9-11 extracted the cell tensor at **normal incidence** (kt = 0). The surface wave sees the tensor at **kt = beta ~ 1.07-1.3 k0** (Phase 30 lesson, `tensor_sheet.bloch_sheet_tensor`). With the loaded Bloch tensor (diode loads added to the MoM matrix, kt = beta) the same RT5880 rim-disc cell gives dY = -3.9j mS and alpha_beam **1.3 Np/m** (diode position 1.6 mm) / **1.8-1.9 Np/m** (position 0.8 mm, mesh 64; mesh 96: 1.74), against 0.7-0.9 Np/m from the normal-incidence tensor. The RT5880 board therefore reaches **21.5 dBi** (not the 18.3 dBi of round 11): the "alpha = 0.7 Np/m, 80 % of the power at the edge" statement of `CODED_IMPROVE.md` was a property of the kt = 0 extraction. Edge power is 52 %, not 79 %. Leakage strength is still the limiting parameter (section 5), but at a higher level.

## 1. RO3010 switched superpixel (primary): 2 x 2 sub-patches, 6.4 mm
Stack RO3010 0.635 mm, eps_r 11.2, tan d 0.0022, 10.05 GHz (the E1 library frequency; E1 does not work at 10.54 GHz, `HFSS_EIGEN_E1.md`), measured MA4AGP907 at 10.05 GHz (ON 2.43 - 3.72j, OFF -562j ohm). A 6.0 mm disc cell is resonant on this stack (lambda_g/2) and gave meaningless sheet values; the superpixel keeps four non-resonant 3.0 mm sub-patches (uniform 0.2 mm gaps, the E1-like sheet) and widens the two central gaps to 0.4 mm over 0.4 mm only where the diodes sit (x diode in the vertical gap at y = +1.6, y diode in the horizontal gap at x = +1.6; footprint 0.4 x 0.4 mm). Every edge is on the 0.1 mm mesh (mesh 64).
Four states, Bloch Y tensor [xx, xy; yx, yy] in mS at the cell's own beta (1.08 k0):
| x/y diode | Y_xx | Y_xy | Y_yy |
|---|---|---|---|
| OFF/OFF | 0.04 + 17.73j | -0.00j | 0.04 + 17.73j |
| ON/OFF | 2.54 + 3.37j | 0.15j | 0.04 + 17.73j |
| OFF/ON | 0.04 + 17.73j | 0.15j | 2.54 + 3.37j |
| ON/ON | 5.12 + 0.49j | -8.02j | 5.12 + 0.49j |
dY per axis (ON - OFF) = **-14.4j + 2.5 mS**; the crossed diode moves the other axis by < 0.01 mS, but **ON/ON is not additive** (Y_xy = -8.0j mS). The ON diode makes the gap neck resonate against the gap capacitance: Y collapses from 17.7 to 3.4 mS (|dY|/Y = 0.8): this is a strong scatterer, **outside the first-order regime** (E1 library: M_t <= 0.66, delta ~ 2 mS). alpha_beam / alpha_abs from the first-order formula depend on the reference sheet and cannot be pinned down:
| definition | alpha_beam [Np/m] | alpha_abs [Np/m] |
|---|---|---|
| about the cell's own mean, x and y sheets averaged (target of the item) | **7.1 (mesh 64) - 8.8 (refined beta)** | 2.5-2.8 |
| about the x-only mean (Y_xx 10.6 mS, beta ~ 1.03 k0), mesh 64 / 96 | 1.15 / 0.84 | 1.8 / 1.5 |
| about the E1 reference sheet (17 mS, beta 1.18 k0), c1 = dY/pi | about 30 (invalid) | NA |
The 5-10 Np/m target is met only by the first definition. Mesh 48 puts the 0.1 mm notch off the grid and the diode drops out (ON = OFF): mesh 64 only (96 is also off-grid for the notch: alpha 0.84 with the x-only mean).

## 2. RT5880 secondary: multi-gap / loaded-slit candidates (Bloch, mesh 64)
| cell | dY [mS] | alpha_beam [Np/m] | alpha_abs | beta/k0 |
|---|---|---|---|---|
| rim disc 6.0, diode at 1.6 mm (round 9-11 cell) | -3.93j | 1.33 | 0.11 | 1.074 |
| **rim disc, diode at 0.8 mm (chosen)** | -4.74j | **1.84** | 0.094 | 1.079 |
| 2 x 2 patches, central gap 0.4, outer 0.4 | +4.18j | 0.33 | 0.12 | 1.035 |
| 2 x 2 patches, central 0.4, outer 0.2 | +5.72j | 1.09 | 0.22 | 1.042 |
| disc + two static extra slits (offset 1.0 mm) | -22.2j | 0.46 (out of regime) | 0.35 | 1.295 |
| disc + two static extra slits (offset 2.0 mm) | -2.3j + 2.1 | 0.63 | 1.68 | 1.093 |
No multi-gap or loaded-slit variant beats the rim disc moved to position 0.8 mm (the 5-10 Np/m target is not reachable on RT5880 with this diode: best 1.8-1.9). Four states of the chosen cell: OFF/OFF 11.22j, ON/OFF 0.07 + 6.47j (xx), ON/ON 0.10 + 6.25j with Y_xy = 0.5j: additive.

## 3. K = 2 boards (diodes <= 766, 1 control line + inverter per axis; MODEL_ONLY)
RO3010: board mean beta 1.14 k0, switched cell A0 = 7.1 Np/m (own-mean definition), AA0 2.5, background loss 0.6 Np/m, static cells = E1-type tensor cells (any axis angle, |m| <= sqrt(13.1/A0), 13.1 Np/m = first-order maximum of the E1 library; 9.2 Np/m Floquet-ratio variant), launch 0.476 (PinLaunch, 6 mm pin on RO3010; the HFSS-calibrated 0.627 belongs to RT5880). RT5880: A0 = 1.84, background 0.15, launch 0.627, option sets of round 11 (I1, I2, I5).
| run | gain A/B [dBi] | SLL [dB] | AR [dB] | ap. eff. [%] | edge [%] | static [%] | SW->beam [%] | diodes | lines |
|---|---|---|---|---|---|---|---|---|---|
| RT5880_bits | [20.56, 20.56] | [-14.2, -14.2] | [0.81, 0.81] | 12.52 | 59.1 | 73.3 | 36.6 | 764 | 1 |
| RT5880_I1I2I5 | [21.48, 21.48] | [-15.5, -15.4] | [0.72, 0.74] | 15.50 | 52.5 | 73.3 | 43.8 | 765 | 1 |
| RO3010_bits | [20.77, 20.77] | [-15.5, -15.5] | [1.06, 1.06] | 14.48 | 10.4 | 73.2 | 70.0 | 766 | 1 |
| RO3010_E1 | [21.68, 21.69] | [-13.7, -13.8] | [0.5, 0.51] | 17.84 | 14.5 | 73.6 | 66.6 | 754 | 1 |
| RO3010_E1_loss | [21.66, 21.66] | [-13.5, -13.5] | [0.51, 0.52] | 17.75 | 14.7 | 73.9 | 66.5 | 746 | 1 |
| RO3010_E1_launch0627 | [22.88, 22.89] | [-13.7, -13.8] | [0.5, 0.51] | 23.50 | 14.5 | 73.6 | 66.6 | 754 | 1 |
| RO3010_E1_floquet | [21.55, 21.54] | [-14.9, -14.9] | [0.34, 0.34] | 17.29 | 16.7 | 73.4 | 63.7 | 760 | 1 |
| RO3010_E1_taper | [21.38, 21.39] | [-15.5, -15.6] | [0.81, 0.81] | 16.66 | 12.7 | 73.4 | 67.6 | 762 | 1 |
`RO3010_E1_launch0627` = same board with the RT5880 launch fraction (launcher-neutral comparison); `RO3010_E1_taper` = depth taper of the E1 cells (uniform-illumination target saturated in every sector).
Equal-diode-count ratio of aperture efficiency, RO3010 / RT5880: **17.8 / 15.5 = 1.15** (own launch fractions), **23.5 / 15.5 = 1.5** at equal launch. **The criterion (>= 2x) is not met -> RT5880 is retained.**
Not carried over to RO3010: I2 (pre-phased pairs), because only one switched geometry works on the mesh; I5 loss-aware assignment gives no gain (17.75 vs 17.84 dBi).

## 4. Sensitivity (eps_r +-0.3, mesh) of the switched-cell alpha and of the board
| cell | eps_r 10.9 / 2.5-0.3 | nominal | eps_r 11.5 / 2.5 | mesh 32 / 64 / 96 |
|---|---|---|---|---|
| RO3010 notch, alpha_beam (x-only mean) [Np/m] | 0.85 | 1.15 | 1.60 | invalid / 1.15 / 0.84 |
| RO3010 notch, dY [mS] | -14.8j | -14.2j | -13.7j | invalid / -14.25j / -13.97j |
| RT5880 disc pos 0.8, alpha_beam [Np/m] (eps_r 1.9 / 2.2 / 2.5) | 0.39 | 1.93 | 3.43 | 2.28 / 1.93 / 1.74 |
| RT5880 disc pos 0.8, dY [mS] (eps_r 1.9 / 2.2 / 2.5) | -3.1j | -4.8j | -7.8j | -5.45j / -4.84j / -4.53j |
Converged (mesh 64 -> 96, 1/n): RT5880 **alpha_beam ~ 1.5-1.7 Np/m, |dY| ~ 4.2-4.5 mS**; RO3010 notch: dY stable (-14 mS), alpha not converged (x-only 0.8-1.2; own-mean 7-9).
Board gain versus alpha (same engine): RT5880 (round-11 option sets) alpha 0.39 / 0.7 / 1.0 / 1.84 / 3.43 Np/m -> **15.6 / 17.9 / 19.2 / 21.2 / 22.8 dBi** (ap. eff. 4.0 / 6.8 / 9.2 / 14.7 / 21 %); RO3010 (E1 statics) alpha 0.85 / 1.15 / 1.6 / 3.0 / 7.1 -> **16.5 / 17.6 / 18.4 / 19.6 / 21.7 dBi** (5.4 / 7.0 / 8.3 / 11 / 17.8 %). eps_r +-0.3 moves alpha by x0.7/x1.4 on RO3010 and x0.2/x1.8 on RT5880 (the +-0.3 requested is far wider than a laminate datasheet tolerance, which was not checked here). The E1 static capacity is held at 13.1 Np/m in these runs (its own eps_r dependence is not modelled).

Board versus switched-cell alpha (rows `RT5880_A0_*`, `RO3010_A0_*`; same columns as section 3):
| run | gain A/B [dBi] | SLL [dB] | AR [dB] | ap. eff. [%] | edge [%] | static [%] | SW->beam [%] | diodes | lines |
|---|---|---|---|---|---|---|---|---|---|
| RT5880_A0_0.39 | [15.59, 15.59] | [-15.1, -15.1] | [0.59, 0.59] | 3.99 | 85.9 | 73.4 | 10.3 | 762 | 1 |
| RT5880_A0_0.7 | [17.9, 17.91] | [-15.2, -15.2] | [0.59, 0.58] | 6.80 | 78.5 | 73.3 | 17.7 | 763 | 1 |
| RT5880_A0_1.0 | [19.22, 19.23] | [-15.5, -15.5] | [0.57, 0.58] | 9.21 | 72.0 | 73.4 | 24.2 | 760 | 1 |
| RT5880_A0_1.84 | [21.24, 21.24] | [-16.3, -16.3] | [0.61, 0.62] | 14.66 | 56.6 | 73.5 | 39.7 | 758 | 1 |
| RT5880_A0_3.43 | [22.8, 22.81] | [-15.8, -15.9] | [0.58, 0.58] | 20.99 | 36.0 | 73.4 | 60.2 | 761 | 1 |
| RO3010_A0_0.85 | [16.5, 16.5] | [-11.0, -11.0] | [0.56, 0.56] | 5.41 | 54.1 | 73.4 | 16.6 | 761 | 1 |
| RO3010_A0_1.15 | [17.59, 17.59] | [-11.1, -11.1] | [0.27, 0.27] | 6.95 | 50.3 | 73.3 | 21.2 | 765 | 1 |
| RO3010_A0_1.6 | [18.37, 18.37] | [-12.1, -12.1] | [0.5, 0.51] | 8.32 | 46.7 | 73.6 | 25.7 | 755 | 1 |
| RO3010_A0_3.0 | [19.6, 19.6] | [-15.5, -15.6] | [0.44, 0.44] | 11.05 | 38.6 | 73.5 | 36.0 | 759 | 1 |
| RO3010_A0_7.125 | [21.68, 21.69] | [-13.7, -13.8] | [0.5, 0.51] | 17.84 | 14.5 | 73.6 | 66.6 | 754 | 1 |

## 5. What limits the board
Gain follows the switched-cell alpha almost dB for dB on both stacks (table above): **leakage strength is the limiting parameter**. RO3010 helps only through the E1 static cells (44 % of the cells, 13 Np/m, user-stated HFSS-validated edge 25 % on the 90 mm board; this repository has no entry for that run, `hfss_results/README.md`) and through a switched cell whose alpha is uncertain by a factor 8 (0.84-8.8 Np/m by definition). RT5880 reaches 21.5 dBi at alpha 1.8 Np/m without any new cell, once the Bloch tensor is used.

## 6. Launcher
Baseline unchanged: bare 6 mm probe (RT5880: launch 0.627, S11 -15 dB, HFSS calibrated; RO3010 6 mm pin: model 0.476, not calibrated). The 10 mm / 3 mm top-hat probe (launch 0.80, S11 -1.3 dB) is an HFSS item only, not in any number above.

## HFSS checklist
RO3010 notch superpixel (4 states, ON collapse Y 17.7 -> 3.4 mS, ON/ON non-additivity, mesh-alignment); E1 static capacity (13.1 vs 9.2 Np/m) on the actual 6.4 mm tiling; RT5880 rim disc at position 0.8 mm; both cells' alpha with a Bloch/supercell run; RO3010 6 mm probe launch fraction and S11; eps_r of the laminate (measured).
