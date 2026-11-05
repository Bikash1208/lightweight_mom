# Crossed rim cell: rim width robustness and mesh convergence (round 10) — MODEL_ONLY

> **PENDING CELL VALIDATION (2026-11-05):** HFSS shows the rim-disc cell does not switch (90 deg along x: OFF 11.059 / ON 11.055 GHz; model 10.53 / 11.00). Every number in this document that depends on the crossed rim-disc cell (alpha, dY, board gain, diode counts) is on hold until a cell is HFSS-confirmed. See `results/sidequest_hybrid/CELL_VALIDATION.md`.

Cell: 6.0 mm disc, two 0.4 mm slits of length 2 x (3.0 - rim), x diode at (0, +1.6) and y diode at (+1.6, 0), measured MA4AGP907. Code `sidequest_coded_rim.py`, `sidequest_coded_rim_board.py`; data `coded_rim/rim*_m*.json`, `coded_rim_board.json`, `coded_rim_extrap.json`.
dY = Y_xx(ON/OFF) - Y_xx(OFF/OFF) (other diode OFF); alpha/beta from the first-order leakage model; board = K = 2 RHCP (25, 0) / LHCP (25, 180), two diodes per cell, lambda 0.5 (eps_B scanned, <= 1 dB). Mesh edge: 0.2 / 0.1 / 0.0667 mm at mesh 32 / 64 / 96.

| rim [mm] | mesh | dY [mS] | alpha_beam [Np/m] | beta/k0 | static | diodes | gain A/B [dBi] | SLL [dB] | AR [dB] | XPD [dB] |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 32 | -3.59j | 1.044 | 1.080 | 71.9 % | 805 | 18.75 / 18.73 | -13.8 / -13.6 | 0.89 / 0.82 | 25.8 / 26.5 |
| 0.2 | 64 | -3.24j | 0.891 | 1.075 | 72.0 % | 801 | 18.13 / 18.17 | -13.3 / -14.1 | 0.89 / 0.86 | 25.8 / 26.1 |
| 0.2 | 96 | -3.06j | 0.811 | 1.073 | 72.4 % | 788 | 17.76 / 17.79 | -13.1 / -13.7 | 0.93 / 0.89 | 25.4 / 25.8 |
| 0.25 | 64 | -2.34j | 0.470 | 1.065 | 73.2 % | 766 | 15.60 / 15.60 | -13.3 | 0.77 | 27.0 |
| 0.25 | 96 | -2.47j | 0.525 | 1.066 | 73.4 % | 762 | 16.05 / 16.05 | -13.2 | 0.79 | 26.8 |
| 0.3 | 64 | -2.34j | 0.470 | 1.065 | 73.2 % | 766 | 15.60 / 15.60 | -13.3 | 0.77 | 27.0 |
| 0.3 | 96 | -2.47j | 0.525 | 1.066 | 73.4 % | 762 | 16.05 / 16.05 | -13.2 | 0.79 | 26.8 |
Mesh 32 cannot resolve a rim change (0.2 / 0.25 / 0.3 give the same metal mask: identical rows). At mesh 64 and 96, 0.25 and 0.3 mm give **identical** masks too (the slit end falls into the same cell), so they are the same data point, "rim >= 0.25 mm", not two.

## 1. Rim choice
Going from 0.2 to >= 0.25 mm costs **1.7 dB (mesh 96) to 2.5 dB (mesh 64)**, far above the 0.5 dB limit; alpha falls by 35-47 %. Interpolating the mesh-96 loss linearly over the 0.05 mm, 0.5 dB would be reached near 0.215 mm, but the available meshes cannot resolve an intermediate rim. **Chosen: 0.2 mm** (the widest rim tested that meets the limit). The board's static %, diodes, SLL, AR, XPD do not move with the rim (+-1.5 points, +-0.5 dB, +-0.1 dB); only the gain does, through alpha.

## 2. Mesh convergence of the chosen cell (rim 0.2 mm)
| mesh | dY [mS] | alpha_beam [Np/m] | alpha_abs [Np/m] | beta/k0 |
|---|---|---|---|---|
| 32 | -3.59j + 0.09 | 1.044 | 0.111 | 1.080 |
| 64 | -3.24j + 0.09 | 0.891 | 0.103 | 1.075 |
| 96 | -3.06j + 0.09 | 0.811 | 0.098 | 1.073 |
Not converged, and not a clean power law (the rim is 2-3 cells wide). Extrapolating 64 -> 96 as a + b/n or a + b/n^2: **|dY| = 2.7-2.9 mS, alpha_beam = 0.65-0.75 Np/m, alpha_abs ~ 0.095, beta/k0 ~ 1.07** (the earlier 0.92 Np/m at mesh 64 with the 4-state mean is the same cell, 4 % apart). Board with the extrapolated alpha: 72.7 % static, 782 diodes, **16.9-17.4 dBi** (vs 18.1 at mesh 64, 17.8 at mesh 96), SLL -13.3 dB, AR 0.84 dB, XPD 26 dB. Quote **17 +- 0.5 dBi** for the converged cell; the 18.2 dBi of round 9 is the mesh-64 value and is optimistic by about 1 dB.
Mesh sensitivity is the main uncertainty of this cell; HFSS decides.

## 3. HFSS package (not sent): `hfss_exports/sidequest_coded_cross_cell/`
Four STEP files (OFF/OFF, ON/OFF, OFF/ON, ON/ON; sheets `Dx_<state>`, `Dy_<state>`), `HFSS_SETUP.md`, `predictions.json` (TM0 band per state, x and y, 90 and 120 deg, from the mesh-64 uniform-sheet model). Rim 0.200 mm. EXACT_DIMS PASS (4 solids + 2 sheets, bboxes/volumes to < 1e-5, 0 inexact literals, no collapsed or merged vertices). The STEP disc is a true circle; the MoM uses a 48-gon (area -0.7 %). Expected splits are small (+0.12 GHz on the on-axis diode at 90 deg).
