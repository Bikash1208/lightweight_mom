# SUB6 board model with ISLAND switched cells (stage: sub6-board-island)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Code `scripts/research/sub6_board_model4.py` (derived from `sub6_board_model2.py`, same Phase 26 engine, real per-state alpha/beta, joint 90-degree rotations, taper, ring exclusion, placements a/b); data `points_island.csv` (216 points), `layout_*_t0.8.csv`, `provenance_run_island.json`, inputs `input_island_rows.json/csv` (copied from compute/sub6-island-cell 6d53878). Launch fraction 0.656 (standard SMA pin, sub6-launcher). 16 s on 64 workers.

## Architecture and cells
STATIC = NOCUT rows (64 rows, filters as in round 1). SWITCHED = island rows of the island scan **without the resonance flag** (16 of 48), further filtered by me to alpha_OFF >= 0.9 Np/m and sqrt(alpha_ON/alpha_OFF) <= 0.5 (near-silent in ON): **9 rows**, alpha_OFF 0.96-2.70 Np/m, beta_OFF 1.08-1.11, beta_ON 1.00-1.05 (each cell's own beta in each state). Static cells radiate beam A in both states; island cells radiate beam B in OFF. General weights, rho_d = mean signed ratio J_ON/J_OFF of the table (about -0.33), c = -rho_d/(1-rho_d): static J = (u_A + c u_B)/(1-f), switched J_OFF = u_B/((1-rho_d) f), J_ON = rho_i J_OFF with rho_i = sign(Im Delta_ON/Im Delta_OFF) sqrt(alpha_ON/alpha_OFF) of the chosen row (my choice: magnitude from alpha, sign from the reactance polarity; the complex phase of Delta_ON/Delta_OFF is not used). The island cells have no A_perp (mirror symmetric), so the rotations act on the Delta term only.
Sweep: N_sw 110/120/135 (diodes on hand) and 180/240/300 (ideal) x sides 180/200/220 mm x placement a/b x taper 0.5/0.8/1.0 x ring exclusion off/on = 216 points.

## Gates (balance <= 3 dB, AR <= 3 dB over the -3 dB beam, B suppression >= 12 dB, SLL <= -13 dB both states)
**MEETS_ALL on the stricter -3 dB-beam AR gate: 23 of 216** (peak-AR gate: 188 of 216). Best point: 220 mm, N_sw 135, placement b, taper 0.5, ring exclusion yes, margin 0.39 dB.

| side | N_sw (f) | best pl / taper / ring-excl | OFF G_A / G_B / bal | OFF AR3dB A / B | ON G_A / B supp | ON AR3dB | SLL off / on | eta off / on | rotated | MEETS_ALL (-3 dB AR / peak AR) of 12 | best gate (margin) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 220 | 110 (0.12) | b / 0.5 / yes | 12.7 / 11.2 / 1.50 | 2.81 / 2.79 | 10.8 / 25.8 | 2.12 | -14.3 / -16.0 | 0.12 / 0.04 | 9 | 5 / 12 | MEETS_ALL (0.19) |
| 220 | 120 (0.13) | a / 0.5 / yes | 13.0 / 11.8 / 1.21 | 2.77 / 2.58 | 10.5 / 15.7 | 2.23 | -14.9 / -17.7 | 0.14 / 0.04 | 8 | 3 / 12 | MEETS_ALL (0.23) |
| 220 | 135 (0.14) | b / 0.5 / yes | 13.7 / 12.7 / 1.01 | 2.61 / 2.19 | 12.4 / 28.2 | 2.12 | -13.7 / -16.3 | 0.16 / 0.06 | 9 | 3 / 11 | MEETS_ALL (0.39) |
| 220 | 180 (0.19) | a / 0.5 / no | 15.6 / 14.8 / 0.80 | 3.11 / 2.60 | 14.7 / 26.0 | 2.06 | -15.7 / -17.2 | 0.24 / 0.10 | 9 | 0 / 12 | MISSES_AR (-0.11) |
| 220 | 240 (0.26) | b / 0.8 / no | 17.2 / 16.5 / 0.64 | 3.17 / 3.06 | 17.0 / 25.9 | 1.77 | -16.6 / -17.4 | 0.34 / 0.18 | 11 | 0 / 12 | MISSES_AR (-0.17) |
| 220 | 300 (0.32) | b / 0.5 / no | 18.3 / 18.0 / 0.31 | 3.14 / 3.20 | 18.4 / 16.5 | 2.28 | -18.1 / -17.9 | 0.43 / 0.25 | 6 | 0 / 12 | MISSES_AR (-0.20) |
| 200 | 110 (0.15) | b / 0.5 / yes | 12.3 / 11.4 / 0.85 | 2.54 / 2.85 | 11.1 / 28.1 | 2.21 | -14.2 / -16.4 | 0.15 / 0.05 | 3 | 6 / 12 | MEETS_ALL (0.15) |
| 200 | 120 (0.16) | a / 0.5 / no | 12.7 / 11.9 / 0.76 | 2.54 / 2.75 | 12.0 / 30.8 | 2.15 | -13.1 / -15.7 | 0.17 / 0.06 | 7 | 4 / 10 | MEETS_ALL (0.06) |
| 200 | 135 (0.18) | a / 1 / no | 13.5 / 12.5 / 0.99 | 2.68 / 2.56 | 12.8 / 28.4 | 2.02 | -13.2 / -17.7 | 0.20 / 0.08 | 10 | 2 / 8 | MEETS_ALL (0.17) |
| 200 | 180 (0.24) | a / 0.5 / no | 15.3 / 14.5 / 0.81 | 2.76 / 3.01 | 15.4 / 28.0 | 2.40 | -15.5 / -19.2 | 0.29 / 0.14 | 6 | 0 / 11 | MISSES_AR (-0.01) |
| 200 | 240 (0.32) | b / 0.5 / yes | 16.7 / 15.7 / 0.97 | 3.08 / 3.12 | 17.3 / 27.4 | 2.19 | -14.1 / -15.9 | 0.39 / 0.24 | 5 | 0 / 12 | MISSES_AR (-0.12) |
| 200 | 300 (0.40) | b / 0.5 / yes | 17.3 / 16.9 / 0.47 | 3.11 / 3.12 | 18.0 / 21.4 | 2.42 | -13.0 / -16.5 | 0.48 / 0.30 | 8 | 0 / 12 | MISSES_AR (-0.12) |
| 180 | 110 (0.18) | b / 0.5 / no | 12.4 / 11.3 / 1.12 | 3.27 / 2.81 | 11.6 / 26.2 | 2.34 | -12.9 / -16.2 | 0.17 / 0.07 | 6 | 0 / 3 | MISSES_AR+SLL (-0.27) |
| 180 | 120 (0.19) | b / 1 / no | 12.9 / 11.9 / 1.03 | 3.38 / 2.95 | 12.1 / 26.2 | 2.21 | -13.7 / -17.2 | 0.19 / 0.08 | 7 | 0 / 8 | MISSES_AR (-0.38) |
| 180 | 135 (0.21) | a / 0.5 / no | 13.7 / 12.9 / 0.82 | 3.45 / 3.22 | 13.4 / 25.7 | 2.09 | -13.9 / -18.6 | 0.23 / 0.10 | 6 | 0 / 12 | MISSES_AR (-0.45) |
| 180 | 180 (0.29) | a / 0.8 / no | 15.3 / 14.9 / 0.41 | 3.98 / 3.43 | 15.6 / 29.0 | 2.17 | -14.3 / -18.1 | 0.34 / 0.19 | 4 | 0 / 12 | MISSES_AR (-0.98) |
| 180 | 240 (0.38) | b / 1 / yes | 15.6 / 16.0 / 0.41 | 3.35 / 2.65 | 16.1 / 20.3 | 2.96 | -13.9 / -16.6 | 0.43 / 0.26 | 4 | 0 / 10 | MISSES_AR (-0.35) |
| 180 | 300 (0.48) | b / 1 / yes | 14.6 / 14.7 / 0.11 | 3.18 / 3.36 | 15.1 / 19.3 | 2.92 | -13.2 / -13.2 | 0.34 / 0.26 | 2 | 0 / 7 | MISSES_AR (-0.36) |

| subset | points | ON supp: unrotated / rotated | ON AR_A 3dB beam | OFF AR_B 3dB beam | SLL_on | MEETS_ALL (-3 dB AR) |
|---|---|---|---|---|---|---|
| all | 216 | 23.1 / 22.8 | 2.24 | 3.15 | -17.0 | 23 |
| ring exclusion on | 108 | 21.5 / 21.5 | 2.24 | 3.14 | -16.7 | 13 |
| ring exclusion off | 108 | 24.7 / 24.2 | 2.24 | 3.16 | -17.2 | 10 |
| taper 0.5 | 72 | 22.5 / 22.3 | 2.26 | 3.16 | -16.9 | 9 |
| taper 0.8 | 72 | 22.5 / 22.3 | 2.25 | 3.16 | -16.9 | 9 |
| taper 1.0 | 72 | 24.3 / 23.9 | 2.21 | 3.13 | -17.1 | 5 |
| placement a | 108 | 22.3 / 22.0 | 2.23 | 3.17 | -17.0 | 11 |
| placement b | 108 | 23.9 / 23.7 | 2.25 | 3.13 | -16.9 | 12 |
| N_sw <= 135 (diodes on hand) | 108 | 24.1 / 23.6 | 2.24 | 3.04 | -17.2 | 23 |
| N_sw >= 180 (ideal) | 108 | 22.2 / 22.1 | 2.23 | 3.26 | -16.8 | 0 |

Rotation lowers B suppression (cap 17 dB set in the optimiser; pre-rotation values are in the table) and buys AR. Gains are referred to accepted power x model eta x launch 0.656.

## Limits
Switched table is only 9 rows (alpha_req is nearly uniform, so one or two rows are used; amplitude mismatch rms 0-8 % at the checked points). rho_d is a table mean, the per-site spread of rho_i leaves a residual B in ON. Island cells: mesh 64 with 0.3-0.6 mm clearances of 2-3 mesh cells (qualitative), no HFSS check, no bias network, floating-island coupling to neighbours not modelled, the OFF-state beta of island cells (1.08-1.11) differs from static cells (1.05-1.2). The AR over the -3 dB beam is from the power formula of `metrics`. Stopped here.
