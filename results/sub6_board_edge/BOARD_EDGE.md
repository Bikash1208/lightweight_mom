# SUB6 board model with EDGE-DIODE switched strings (stage: sub6-board-edge, + EDGE_TO_EDGE amendment)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED. NEED_ANGLE_TABLE (see below): every gate on this page is provisional, because the 90-deg table row is applied to cells whose hologram angle to the string is not in {0, 45, 90} deg.**
Code `scripts/research/sub6_board_edge.py`, `sub6_board_edge2.py` (chords); they import `sub6_board_island2.py` and `sub6_board_model4.py` unchanged (V3 method as the island follow-up). Data `points_edge.csv` (192 rows = 48 designs x V0-V3), `points_edge_chords.csv`, `provenance_run_edge.json`, input `input_edge_diode.json/csv` (pushed commit 34c4a37 of compute/sub6-edge-diode-cell). Launch fraction 0.656.

## NEED_ANGLE_TABLE
The edge-diode table has the slot at theta = 0, 45 and 90 deg to the diode axis; only the 90-deg rows are used here. The hologram sets the slot axis of every cell, so the angle to its string takes all values: for the switched cells of the V3 designs the angle is within 7.5 deg of 0, 45 or 90 deg for only 20-41 % of the cells (mean over designs 34 %); within 7.5 deg of 90 deg only about 20 % (histograms in the table below). The remaining cells need an angle table (alpha, beta, Delta of the edge cell versus the slot angle) from the cell session; **nothing was interpolated or guessed, and the earlier "theta-aware amplitude factor" (0.65-0.70) is withdrawn: it was an interpolation between three rows.** Until the angle table exists the B-beam numbers are not valid for those cells.

## Free segments (10-24 cells, series strings along x or y)
STATIC = NOCUT rows (64); SWITCHED = edge-diode theta-90 rows filtered by me to alpha_OFF >= 0.9 and sqrt(alpha_ON/alpha_OFF) <= 0.5: 5 rows ['EDGE_th90_q0.7_W2.4', 'EDGE_th90_q0.7_W2.9', 'EDGE_th90_q0.7_W3.6', 'EDGE_th90_q0.8_W2.4', 'EDGE_th90_q0.8_W2.9'] (beta_OFF 1.20-1.30, beta_ON 1.017). Edge rows carry no resonance proxy (resonance_flip = null: not computed, not "no resonance").
MEETS_ALL: 1 of 192 rows. | variant | designs | MEETS_ALL | mean margin | best margin |
|---|---|---|---|---|
| V0 | 48 | 0 | -4.52 | -0.05 |
| V1 | 48 | 0 | -3.64 | -0.19 |
| V2 | 48 | 0 | -4.03 | -0.45 |
| V3 | 48 | 1 | -3.22 | 0.33 |

Best V3 design per (side, N, placement):
| side | N_sw | place | best V3 config (taper / excl / w) | OFF G_A / G_B / bal | OFF AR_A / AR_B (-3 dB) | ON G_A / B supp / AR_A | SLL off / on | eta off / on | lines (x) len | string V min-max | diodes | gate (margin) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 220 | 135 | a | 0.5 / no / 0 | 12.9 / 13.4 / 0.53 | 1.82 / 1.89 | 11.4 / 17.6 / 1.75 | -11.5 / -12.9 | 0.20 / 0.06 | 8 (7) 11-21 | 9.5-19.0 V | 127 | MISSES_SLL (-1.48) |
| 220 | 135 | b | 0.5 / yes / 0 | 13.6 / 12.6 / 0.95 | 1.70 / 1.75 | 12.1 / 19.9 / 2.12 | -12.8 / -17.5 | 0.17 / 0.06 | 7 (4) 10-24 | 8.6-21.9 V | 128 | MISSES_SLL (-0.17) |
| 220 | 240 | a | 0.5 / no / 0.75 | 16.7 / 16.2 / 0.52 | 1.79 / 1.82 | 16.8 / 14.7 / 1.85 | -14.7 / -13.3 | 0.36 / 0.17 | 19 (4) 10-20 | 8.6-18.1 V | 221 | MEETS_ALL (0.33) |
| 220 | 240 | b | 0.5 / yes / 0.5 | 15.9 / 16.0 / 0.05 | 1.72 / 1.81 | 16.0 / 17.1 / 2.37 | -9.1 / -12.6 | 0.34 / 0.18 | 13 (9) 10-24 | 8.6-21.9 V | 227 | MISSES_SLL (-3.87) |
| 200 | 135 | a | 1 / yes / 0 | 12.2 / 11.8 / 0.42 | 1.72 / 2.41 | 11.7 / 15.8 / 2.18 | -11.3 / -15.9 | 0.20 / 0.07 | 11 (6) 10-15 | 8.6-13.3 V | 124 | MISSES_SLL (-1.67) |
| 200 | 135 | b | 0.5 / yes / 0 | 13.3 / 11.8 / 1.46 | 2.09 / 2.32 | 12.0 / 18.6 / 2.61 | -9.9 / -10.3 | 0.20 / 0.08 | 7 (4) 10-24 | 8.6-21.9 V | 128 | MISSES_SLL (-3.09) |
| 200 | 240 | a | 1 / no / 0.25 | 14.5 / 16.6 / 2.14 | 1.79 / 1.84 | 15.1 / 15.3 / 2.16 | -10.7 / -8.9 | 0.45 / 0.20 | 19 (2) 10-21 | 8.6-19.0 V | 221 | MISSES_SLL (-4.11) |
| 200 | 240 | b | 1 / no / 0 | 13.9 / 17.7 / 3.79 | 2.01 / 2.03 | 16.1 / 10.3 / 2.06 | -14.3 / -11.0 | 0.41 / 0.25 | 13 (2) 10-24 | 8.6-21.9 V | 227 | MISSES_BALANCE+SUPPRESSION+SLL (-2.00) |

### Beta-shift effect on beam A and slot-to-string angles
| side | N_sw | place | G_A ON, beta shift off / on | cost (dB) | sw d_beta (ON-OFF) | static dphi rms / max (deg) | slot-to-string angle: mean (deg), fraction within 7.5 deg of 90 / of 0, 45 or 90 | hist 15-deg bins (0-90) |
|---|---|---|---|---|---|---|---|---|
| 220 | 135 | a | 11.9 / 11.4 | 0.46 | -0.188 | 31.6 / 123 | 40, 0.07 / 0.31 | 27/28/26/24/11/19 |
| 220 | 135 | b | 12.3 / 12.1 | 0.26 | -0.188 | 19.4 / 69 | 42, 0.06 / 0.31 | 31/18/25/19/18/24 |
| 220 | 240 | a | 16.1 / 16.8 | -0.71 | -0.188 | 37.8 / 260 | 45, 0.07 / 0.31 | 39/38/44/34/42/43 |
| 220 | 240 | b | 16.1 / 16.0 | 0.02 | -0.188 | 33.4 / 124 | 45, 0.08 / 0.35 | 34/45/34/53/44/30 |
| 200 | 135 | a | 12.1 / 11.7 | 0.36 | -0.188 | 23.5 / 83 | 38, 0.07 / 0.36 | 38/28/16/18/17/18 |
| 200 | 135 | b | 12.4 / 12.0 | 0.41 | -0.188 | 22.1 / 96 | 44, 0.07 / 0.31 | 19/31/22/18/26/19 |
| 200 | 240 | a | 15.8 / 15.1 | 0.66 | -0.190 | 64.5 / 229 | 42, 0.04 / 0.29 | 37/53/45/34/42/29 |
| 200 | 240 | b | 17.1 / 16.1 | 0.97 | -0.204 | 33.9 / 90 | 45, 0.09 / 0.31 | 47/40/28/36/44/45 |

The switched cells' beta falls by 0.19-0.20 from OFF to ON (1.20-1.30 to 1.017); a ray over a 24-cell string accumulates about 370 deg of phase difference between the states; static-cell phase error rms 19-65 deg; cost to G_A,ON 0.3-1.0 dB (one design -0.7 dB).

## EDGE_TO_EDGE: every switched line is a full chord of the hexagon (both string ends at the rim)
Chords = complete lattice rows (x) or columns (y) from rim to rim; the row j = 0 and column i = 0 are interrupted by the feed clearance and are not chords; with ring exclusion on, no chord within 2 rows/columns of the feed. Chords are chosen greedily on the ideal-cell projection objective until the cell count is close to the target (135 / 240); x and y chords may not share a cell (a crossing is not a series string); y-only and x+y gave the same set whenever the greedy picks only y chords. Same V3 evaluation.
| side | target N | chords | actual N | lines (x) | chord len min-max | longest string (diodes, V) | total diodes | best gate (margin) | free segments at the same N: best gate (margin) |
|---|---|---|---|---|---|---|---|---|---|
| 200 | 135 | x only | 126 | 6 (6) | 5-29 | 28 diodes, 26.6 V | 120 | MISSES_SLL (-3.16) | MISSES_BALANCE+SLL (-3.08, place a) |
| 200 | 135 | y only | 137 | 5 (0) | 19-31 | 30 diodes, 28.5 V | 132 | MISSES_AR+SUPPRESSION+SLL (-4.98) | MISSES_SLL (-2.57, place a) |
| 200 | 135 | x+y | 128 | 6 (2) | 5-31 | 30 diodes, 28.5 V | 122 | MISSES_AR+SLL (-6.67) | MISSES_SLL (-3.32, place b) |
| 200 | 240 | x only | 223 | 11 (11) | 5-29 | 28 diodes, 26.6 V | 212 | MISSES_SLL (-4.40) | MISSES_SLL (-3.65, place b) |
| 200 | 240 | y only | 232 | 8 (0) | 19-33 | 32 diodes, 30.4 V | 224 | MISSES_BALANCE+SLL (-7.60) | MISSES_SLL (-0.33, place a) |
| 200 | 240 | x+y | 232 | 8 (0) | 19-33 | 32 diodes, 30.4 V | 224 | MISSES_BALANCE+SLL (-7.60) | MISSES_SLL (-0.33, place a) |
| 220 | 135 | x only | 129 | 7 (7) | 3-33 | 32 diodes, 30.4 V | 122 | MISSES_SLL (-2.74) | MISSES_SLL (-2.54, place b) |
| 220 | 135 | y only | 128 | 4 (0) | 19-37 | 36 diodes, 34.2 V | 124 | MISSES_SLL (-9.07) | MISSES_SLL (-5.56, place a) |
| 220 | 135 | x+y | 128 | 4 (0) | 19-37 | 36 diodes, 34.2 V | 124 | MISSES_SLL (-9.07) | MISSES_SLL (-5.56, place a) |
| 220 | 240 | x only | 229 | 11 (11) | 3-33 | 32 diodes, 30.4 V | 218 | MISSES_SLL (-2.98) | MEETS_ALL (1.18, place a) |
| 220 | 240 | y only | 223 | 7 (0) | 19-37 | 36 diodes, 34.2 V | 216 | MISSES_SLL (-5.39) | MISSES_AR+SLL (-1.40, place b) |
| 220 | 240 | x+y | 223 | 7 (0) | 19-37 | 36 diodes, 34.2 V | 216 | MISSES_SLL (-5.39) | MISSES_AR+SLL (-1.40, place b) |

Chord designs MEETS_ALL: 0 of 72; free segments at the same N: 8 of 108. **Longest strings: 28 diodes = 26.6 V (x chord, 200 mm), 32 diodes = 30.4 V (y chord, 200 mm), up to 36 diodes = 34.2 V (y chord, 220 mm)**; a full chord is thus always too long for a 24-cell string, the required voltage is reported regardless. Chord total diodes 120-224 (150 on hand: chords at N about 130 fit, N about 225 does not).

## Limits
NEED_ANGLE_TABLE (above); cell values are periodic-lattice values (every cell has a diode and a 1.00 mm gap, disc r 5.30), string ends and neighbours are not modelled; no HFSS check; no bias supply modelled; the chord selection and the row filter are mine. Stopped here.
