# SUB6 board model with EDGE-DIODE switched strings (stage: sub6-board-edge)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Code `scripts/research/sub6_board_edge.py` (imports `sub6_board_island2.py` and `sub6_board_model4.py` unchanged; V3 method as the island follow-up); data `points_edge.csv` (192 rows = 48 designs x V0-V3), `provenance_run_edge.json`, input `input_edge_diode.json/csv` (compute/sub6-edge-diode-cell 34c4a37). Launch fraction 0.656. 13 s on 48 workers.

## Setup
STATIC = NOCUT rows (64); SWITCHED = edge-diode rows with the slot at 90 deg to the diode axis, real OFF/ON alpha, beta and Delta, filtered by me to alpha_OFF >= 0.9 Np/m and sqrt(alpha_ON/alpha_OFF) <= 0.5: **5 rows** ['EDGE_th90_q0.7_W2.4', 'EDGE_th90_q0.7_W2.9', 'EDGE_th90_q0.7_W3.6', 'EDGE_th90_q0.8_W2.4', 'EDGE_th90_q0.8_W2.9'] (alpha_OFF 1.28-2.96, alpha_ON 0.21, **beta_OFF 1.20-1.30, beta_ON 1.017**). The edge-diode scan has no resonance proxy. Architecture, weights (rho_d about -0.3), rotations and V0-V3 as in the island follow-up. Switched cells lie in straight lines of consecutive cells along x or y (length 10-24, no overlap, no end-to-end contact of same-axis lines), each line one series string of n-1 diodes (0.95 V each at 10 mA); placement a = random compatible lines, b = greedy + remove/replace local search on the ideal-cell projection objective (exactly N_sw cells). Sweep: sides 200 / 220 mm, N_sw 135 / 240, taper 0.5/0.8/1.0, ring exclusion off/on.

## Gates
**MEETS_ALL: 1 of 192 rows** (220 mm N240 a V3); the misses are SLL (OFF-state sidelobes of the string structure, -9 to -14 dB) and, at N 240, balance / suppression.
| variant | designs | MEETS_ALL | mean margin | best margin | gate counts |
|---|---|---|---|---|---|
| V0 | 48 | 0 | -4.52 | -0.05 | {'MISSES_SLL': 28, 'MISSES_AR+SLL': 17, 'MISSES_BALANCE+SUPPRESSION+SLL': 3} |
| V1 | 48 | 0 | -3.64 | -0.19 | {'MISSES_SLL': 35, 'MISSES_AR+SLL': 10, 'MISSES_SUPPRESSION+SLL': 3} |
| V2 | 48 | 0 | -4.03 | -0.45 | {'MISSES_SLL': 28, 'MISSES_AR+SLL': 17, 'MISSES_SUPPRESSION+SLL': 3} |
| V3 | 48 | 1 | -3.22 | 0.33 | {'MISSES_SLL': 43, 'MISSES_BALANCE+SUPPRESSION+SLL': 3, 'MEETS_ALL': 1, 'MISSES_AR+SLL': 1} |

Best V3 design per (side, N, placement):
| side | N_sw | place | best V3 config (taper / excl / w) | OFF G_A / G_B / bal | OFF AR_A / AR_B (-3 dB) | ON G_A / B supp / AR_A | SLL off / on | eta off / on | rotated | lines (x) len min-max | string V min-max | diodes | gate (margin) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 220 | 135 | a | 0.5 / no / 0 | 12.9 / 13.4 / 0.53 | 1.82 / 1.89 | 11.4 / 17.6 / 1.75 | -11.5 / -12.9 | 0.20 / 0.06 | 19 | 8 (7) 11-21 | 9.5-19.0 V | 127 | MISSES_SLL (-1.48) |
| 220 | 135 | b | 0.5 / yes / 0 | 13.6 / 12.6 / 0.95 | 1.70 / 1.75 | 12.1 / 19.9 / 2.12 | -12.8 / -17.5 | 0.17 / 0.06 | 14 | 7 (4) 10-24 | 8.6-21.9 V | 128 | MISSES_SLL (-0.17) |
| 220 | 240 | a | 0.5 / no / 0.75 | 16.7 / 16.2 / 0.52 | 1.79 / 1.82 | 16.8 / 14.7 / 1.85 | -14.7 / -13.3 | 0.36 / 0.17 | 15 | 19 (4) 10-20 | 8.6-18.1 V | 221 | MEETS_ALL (0.33) |
| 220 | 240 | b | 0.5 / yes / 0.5 | 15.9 / 16.0 / 0.05 | 1.72 / 1.81 | 16.0 / 17.1 / 2.37 | -9.1 / -12.6 | 0.34 / 0.18 | 17 | 13 (9) 10-24 | 8.6-21.9 V | 227 | MISSES_SLL (-3.87) |
| 200 | 135 | a | 1 / yes / 0 | 12.2 / 11.8 / 0.42 | 1.72 / 2.41 | 11.7 / 15.8 / 2.18 | -11.3 / -15.9 | 0.20 / 0.07 | 8 | 11 (6) 10-15 | 8.6-13.3 V | 124 | MISSES_SLL (-1.67) |
| 200 | 135 | b | 0.5 / yes / 0 | 13.3 / 11.8 / 1.46 | 2.09 / 2.32 | 12.0 / 18.6 / 2.61 | -9.9 / -10.3 | 0.20 / 0.08 | 8 | 7 (4) 10-24 | 8.6-21.9 V | 128 | MISSES_SLL (-3.09) |
| 200 | 240 | a | 1 / no / 0.25 | 14.5 / 16.6 / 2.14 | 1.79 / 1.84 | 15.1 / 15.3 / 2.16 | -10.7 / -8.9 | 0.45 / 0.20 | 14 | 19 (2) 10-21 | 8.6-19.0 V | 221 | MISSES_SLL (-4.11) |
| 200 | 240 | b | 1 / no / 0 | 13.9 / 17.7 / 3.79 | 2.01 / 2.03 | 16.1 / 10.3 / 2.06 | -14.3 / -11.0 | 0.41 / 0.25 | 9 | 13 (2) 10-24 | 8.6-21.9 V | 227 | MISSES_BALANCE+SUPPRESSION+SLL (-2.00) |

## Beta-shift effect on beam A and the slot-angle caveat
| side | N_sw | place | G_A ON, beta shift OFF (OFF phase) / ON | cost (dB) | sw d_beta (ON-OFF) | static dphi rms / max (deg) | theta(slot, line) mean (deg), fraction >= 67.5 | theta-aware amplitude factor (mean, fraction < 0.5) |
|---|---|---|---|---|---|---|---|---|
| 220 | 135 | a | 11.9 / 11.4 | 0.46 | -0.188 | 31.6 / 123 | 40, 0.19 | 0.66, 0.19 |
| 220 | 135 | b | 12.3 / 12.1 | 0.26 | -0.188 | 19.4 / 69 | 42, 0.24 | 0.68, 0.21 |
| 220 | 240 | a | 16.1 / 16.8 | -0.71 | -0.188 | 37.8 / 260 | 45, 0.26 | 0.67, 0.23 |
| 220 | 240 | b | 16.1 / 16.0 | 0.02 | -0.188 | 33.4 / 124 | 45, 0.22 | 0.65, 0.23 |
| 200 | 135 | a | 12.1 / 11.7 | 0.36 | -0.188 | 23.5 / 83 | 38, 0.21 | 0.68, 0.20 |
| 200 | 135 | b | 12.4 / 12.0 | 0.41 | -0.188 | 22.1 / 96 | 44, 0.22 | 0.66, 0.21 |
| 200 | 240 | a | 15.8 / 15.1 | 0.66 | -0.190 | 64.5 / 229 | 42, 0.17 | 0.65, 0.19 |
| 200 | 240 | b | 17.1 / 16.1 | 0.97 | -0.204 | 33.9 / 90 | 45, 0.29 | 0.70, 0.16 |

The switched cells' beta drops by 0.19-0.20 from OFF to ON (1.20-1.30 to 1.017); a ray crossing a string of 24 cells accumulates up to 24 x k0 x 0.19 x 11.6 mm = about 370 deg of phase difference between the states, so the static-cell phase error is large (rms 19-65 deg) and G_A,ON changes by up to 1 dB.
The table row (slot at 90 deg to the diode axis) is used for every switched cell, but the hologram sets the slot axis per cell: the mean angle between slot and string is about 40-45 deg and only about 20 % of the switched cells are within 22.5 deg of 90 deg; with a theta-aware alpha (rows theta = 0 / 45 / 90, same q and W) the mean switched amplitude would be only about 0.65-0.70 of the value used here, and about 20 % of the switched cells would fall below half amplitude. **The results are therefore optimistic for B; this is the main limit of the model.**

## Limits
Cell values are periodic-lattice values (every cell has a diode to its +x neighbour, disc r 5.30, gap 1.00); string ends and neighbours of a different kind are not modelled; the edge-diode rows are not resonance-checked; no bias supply, no HFSS check; the rows' sheet resolution is 6 x 2 mesh cells. Strings of 10-24 cells need 8.6-21.9 V at 10 mA (all strings individually, more than one string cannot be in series unless adjacent). Diode count N_sw - n_lines: 124-128 for N 135 (within the 150 on hand), 221-227 for N 240 (above 150). Stopped here.
