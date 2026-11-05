# Grouped coded hologram (K fixed states, 1-bit, shared control lines) — model level, MODEL_ONLY

> **PENDING CELL VALIDATION (2026-11-05):** HFSS (air box, pec metal) shows the crossed two-diode rim-disc cell switches ONLY with both diodes ON (C/C 10.670, R/C 10.665, C/R 10.690, R/R 11.006 GHz at 90 deg); the model's per-axis independence is falsified. Every number in this document that depends on independent per-axis switching of this cell (alpha per axis, polarisation switching, board gain, diode counts) is on hold. See `results/sidequest_hybrid/CELL_INDEPENDENCE.md`.


Scripts: `sidequest_coded_cell.py` (cell pair), `sidequest_coded_model.py` (board model + greedy), `sidequest_coded_run.py`; tests `tests/test_sidequest_coded.py`; data `coded_cell_scan.json`, `coded_cell_mesh64.json`, `coded_hologram_results.json`.
Reference idea: 1-bit coded hologram with one MA4AGP907 per cell (Zhao/Pan/Shen, IEEE TAP 2026, 10.1109/TAP.2025.3621516; bibliographic-only, not opened).

## Cell pair (one diode, two impedance states, measured MA4AGP907, Phase 26 stack, 10.54 GHz)
Split patch 6.0 x 4.8 mm (neighbour gap 0.4 mm), diode gap 0.4 mm, footprint 1.2 mm, current along x (0 deg only; the rotated-load bug is avoided); screened 27 geometries (mesh 32), chosen by alpha near the available 3.4 Np/m and a low ohmic share.
| state | Y_xx [mS] (mesh 64) | Y_yy [mS] |
|---|---|---|
| OFF (28 fF) | 0.002 + 2.985j | 3.885j |
| ON (2.44 ohm) | 1.199 + 10.593j | 3.885j |
Contrast: dB = 7.6 mS, |dY| / |Y_mean| = 1.13 (mesh 32: 7.5 mS). Mean sheet beta/k0 = 1.054; first-order leakage (repo `leakage_alpha_first_order`, c1 = |dY|/pi, theta 25 deg, zeta 0) alpha_beam = **4.0 Np/m** (3.7 on mesh 32),
ohmic absorption of the ON half alpha_abs = 1.26 Np/m (mesh 32: 1.18). Along a ray at angle psi to the diode axis the depth scales as cos^2 psi (alpha_beam ~ cos^4, alpha_abs ~ cos^2); the TM share cos^2 zeta of the radiated harmonic is co-polar (theta_hat), the rest is cross-polar and not counted.
Only Y_xx switches; Y_yy is static, so a cell modulates the rays along its axis: **x-axis cells for |cos psi| >= |sin psi|, y-axis cells elsewhere (still one diode per cell)**, the same board for K = 2 and K = 4.

## Method
Centre probe, 150 mm hexagon, 1430 cells of 6.4 mm, every cell initially a diode cell. State k code s_ki = sign cos(beta rho_i - k0 sin(theta) u_k.r_i + eps_k); a cell is static if its code is the same in all states.
Gain = G_ideal x |peak amplitude ratio|^2, G_ideal = D_ideal x 0.627 (launch) x converted beam power (alpha-limited, ohmic loss included); currents (pi/2) a s e^{-j beta rho} so that the fundamental equals the ideal hologram.
Joint optimisation: phase offsets eps_k (grid search) + greedy relaxation of the cheapest cells to static words within a coherent-amplitude budget (1 dB), then merging of groups (complement pairs share a line + inverter) with the remaining budget; `frac_static` splits the budget between the two phases.
Full-pattern evaluation of every reported solution (1 deg x 2 deg grid). Polarisation: scalar TM cell = LP (theta_hat), ideal by construction.

## K = 2: (25, 0) and (25, 180), LP
| solution | static cells | diodes | control lines | gain per state | SLL | pointing |
|---|---|---|---|---|---|---|
| ideal codes (no relaxation) | 48.7 % | 734 | 1 (+ inverter; 2 without) | 18.8 / 18.8 dBi | -14.6 / -14.6 dB | (25.0, 0) / (25.0, 180) |
| **recommended: eps2 = 202 deg, 25 % of budget to static** | **71.7 %** | **405** | **1 (+ inverter)** | 18.6 / 18.5 dBi (-0.2 / -0.3 dB) | -15.5 / -13.7 dB | (25.0, 0) / (25.0, 180) |
| 50 % | 77.0 % | 329 | 1 | 18.3 / 18.3 dBi | -12.3 / -11.3 dB | same |
| static-max (1 dB) | 85.1 % | 213 | 1 | 17.8 / 17.8 dBi (-1.0 dB) | **-8.5 / -8.0 dB** | same |
Words: 01 and 10 (complementary) -> one line. The pointing is exact to the 1 deg grid. The limit is the SLL, not the 1 dB gain: beyond ~72 % static the sidelobes rise by more than 3 dB. Gain is 18.8 dBi because only 13 % of the surface-wave power ends in the beam
(ohmic loss in the ON half, cos^4 depth off-axis, 1-bit); the Phase 26 tensor-CP board (23 dBi, HFSS) is a different cell.

## K = 4: (25, 0/90/180/270), LP
| solution | static cells | diodes | control lines | gain per state | SLL | pointing |
|---|---|---|---|---|---|---|
| ideal codes | 14.8 % | 1218 | 7 (+ inverters; 14 without) | 18.8 / 18.8 / 18.8 / 18.8 dBi | -14.6 / -13.7 | exact (25.0, k 90) |
| **recommended: eps = (0, 0, 180, 180) deg, min-lines** | **53.3 %** | **668** | **3 (+ inverters; 6 without)** | 18.0 / 17.9 / 18.0 / 17.9 dBi (-0.8 dB) | -16.2 / -15.1 dB | exact |
| 25 % / 50 % | 39.3 % / 47.6 % | 868 / 750 | 3 / 3 | 18.4 / 18.2 dBi | -15.7 ... -14.7 | exact |
| static-max (1 dB) | 58.9 % | 588 | 7 (14) | 17.8 dBi (-1.0 dB) | -13.1 / -13.5 dB | exact |
Compared with the gate approach: 405 (K = 2) / 668 (K = 4) diodes and 1 / 3 lines instead of 68 gate diodes, but the whole aperture is switched and the ideal-code static fraction is only 49 % / 15 %.

## Polarisation per state (<= 5 lines)
1. Scalar one-diode cells radiate only the TM (theta_hat) polarisation; a different polarisation per state needs a cell that controls two orthogonal admittance components.
2. 2 diodes per cell (x and y axes, both 0/90 deg): the x and y modulations radiate theta_hat / phi_hat parts independently (LP-H/V per state, CP with a quadrature offset between the two codes): doubles the diodes and words (K-bit words per axis).
3. Rotated cells (axis at 45 deg or a continuously rotated slot) give CP or slanted LP with one diode per cell but need rotated lumped loads: **the Phase 28 bug (rotated loads give ON = OFF at non-90 deg) must be fixed first**; 0/90 deg designs do not need the fix.
4. Either route gives 2^K words per axis or per orientation, so the line count grows (up to 2 x (2^K - 2)).
Limits of this study: first-order leakage, cos^2 anisotropy scaling and the TM-share rule are model assumptions; mesh sensitivity of the cell (alpha +8 % from mesh 32 to 64); ideal polarisation/XPD; no edge diffraction; HFSS not run.
