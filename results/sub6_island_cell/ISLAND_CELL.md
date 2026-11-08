# SUB6 island cell scan (stage: sub6-island-cell)

**MODEL_ONLY / NOT_HFSS_VALIDATED.** Slotted disc (anchor r 5.435, p 11.6, RT5880 3.175, 5.8 GHz) + floating metal island in the slot + ONE BAR50 sheet (1.00 x 0.35 mm) across the 1.00 mm clearance between the island's +y long edge and the slot's +y edge. Code `scripts/research/sub6_island_cell.py` (+ `sub6_slotdisc_lib.py`, `sub6_slotdisc_grid.py` unchanged); data `rows.csv` (48 rows), `island_rows.json` (full states, references), `ckpt/` (MoM nodes), `provenance_run.json`.

## Geometry conventions (mine, flagged)
* **Mesh 64 on p = 11.6 mm: one mesh cell = 0.18125 mm** (the mask is sampled at cell centres, so every edge is snapped to a multiple of it). Slot half-length and half-width = floor(dim/2/0.18125 + 1/2) cells; island half-length = round(Li/2/0.18125); clearance 0.3 mm -> 2 cells = **0.3625 mm**, 0.6 mm -> 3 cells = **0.544 mm**; the 1.00 mm +y gap -> 5 cells = **0.906 mm** (same convention as the 1.00 mm rim notch of the cell scan, sheet length 1.00 mm). The realised dimensions are in `rows.csv`. A 0.3 mm gap is only 1.7 mesh cells: resolution-limited, treat the clearance dependence as qualitative.
* The swept clearance applies to the island's -y long edge (the "other side"); the clearance to the two short slot ends follows from Li/L ((L - Li)/2, 1.1-2.7 mm), not from the 0.3/0.6 sweep.
* ISLAND_ABSENT = the closed slotted disc (NOCUT), one reference per (W, L/2r) (8 references); W = 4.4 with L/2r = 0.9 skipped (slot corner leaves 0.07 mm of rim, below the 0.2 mm rule of the cell scan).
* States: OFF_BAR50 ((0.15 pF || 5 k) + 0.6 nH), ON_BAR50 (3.0 ohm + 0.6 nH), as in the cell scan. rho = Delta_ON / Delta_OFF (complex). Gate OK: alpha_ON or alpha_OFF >= 1 Np/m AND |rho - 1| >= 1.

## Result
Gate OK: **28 of 48** rows. Leakage classes: {'OK_ON_ONLY': 6, 'LOW': 16, 'OK_BOTH': 3, 'OK_OFF_ONLY': 23}. Resonance proxy flagged: 32 rows.

Top 12 by max(alpha_OFF, alpha_ON):
| id | real slot W / island Li / clear -y (mm) | alpha OFF / ON | beta OFF / ON | Delta OFF / ON (mS) | |A_perp| OFF / ON | rho | |rho-1| | class | resonance | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| ISL_W2.9_q0.8_Li0.3_c0.6 | 2.900 / 2.538 / 0.544 | 2.33 / 3.02 | 1.0994 / 1.1493 | -0.06588-6.76j / -0.4957-12.8j | 1.44e-14 / 1.04e-12 | 1.89-0.0549j | 0.891 | OK_BOTH | FLIP:ON_BAR50:y | NOT_OK |
| ISL_W3.6_q0.8_Li0.3_c0.3 | 3.625 / 2.538 / 0.362 | 3 / 0.12 | 1.1377 / 1.0098 | -0.2406-11.3j / -3.983+23.3j | 8.49e-13 / 3.2e-13 | -2.05-0.396j | 3.08 | OK_OFF_ONLY | FLIP:ON_BAR50:y | OK |
| ISL_W3.6_q0.8_Li0.5_c0.6 | 3.625 / 4.350 / 0.544 | 2.98 / 0.32 | 1.1407 / 1.0320 | -0.302-11.6j / -0.7689+6.17j | 1.61e-12 / 5.03e-14 | -0.529-0.0799j | 1.53 | OK_OFF_ONLY | FLIP:OFF_BAR50:y | OK |
| ISL_W3.6_q0.8_Li0.5_c0.3 | 3.625 / 4.350 / 0.362 | 2.96 / 0.251 | 1.1546 / 1.0391 | -0.4473-13.3j / -0.4679+3.9j | 6.01e-13 / 2.47e-14 | -0.291-0.0449j | 1.29 | OK_OFF_ONLY | FLIP:OFF_BAR50:y | OK |
| ISL_W3.6_q0.8_Li0.7_c0.6 | 3.625 / 6.162 / 0.544 | 2.95 / 0.112 | 1.1445 / 1.0470 | -0.396-12j / -0.3377+1.92j | 4.57e-12 / 2.42e-14 | -0.158-0.0333j | 1.16 | OK_OFF_ONLY | FLIP:OFF_BAR50:y,ON_BAR50:y | OK |
| ISL_W3.6_q0.8_Li0.3_c0.6 | 3.625 / 2.538 / 0.544 | 2.94 / 0.00924 | 1.1287 / 1.0015 | -0.1742-10.2j / -149+86.1j | 4.48e-12 / 9.42e-12 | -8.18-14.7j | 17.4 | OK_OFF_ONLY | FLIP:ON_BAR50:y | OK |
| ISL_W3.6_q0.8_Li0.7_c0.3 | 3.625 / 6.162 / 0.362 | 2.91 / 0.0532 | 1.1590 / 1.0503 | -0.5851-13.8j / -0.2585+1.19j | 9.19e-13 / 3.74e-15 | -0.0853-0.0223j | 1.08 | OK_OFF_ONLY | FLIP:OFF_BAR50:y | OK |
| ISL_W2.9_q0.8_Li0.7_c0.3 | 2.900 / 6.162 / 0.362 | 2.7 / 0.155 | 1.1118 / 1.0439 | -0.1926-8.23j / -0.4237+2.52j | 1.2e-12 / 1.13e-14 | -0.305-0.0586j | 1.31 | OK_OFF_ONLY | none | OK |
| ISL_W2.9_q0.8_Li0.5_c0.3 | 2.900 / 4.350 / 0.362 | 2.66 / 0.284 | 1.1091 / 1.0218 | -0.1441-7.95j / -1.887+10.8j | 1.51e-12 / 1.72e-13 | -1.36-0.262j | 2.37 | OK_OFF_ONLY | none | OK |
| ISL_W2.9_q0.7_Li0.5_c0.3 | 2.900 / 3.625 / 0.362 | 0.68 / 2.65 | 1.0746 / 1.1213 | -0.02884-3.1j / -0.9887-8.98j | 2.22e-14 / 2.02e-12 | 2.9-0.292j | 1.92 | OK_ON_ONLY | FLIP:ON_BAR50:y | OK |
| ISL_W4.4_q0.8_Li0.3_c0.6 | 4.350 / 2.538 / 0.544 | 2.62 / 0.302 | 1.1852 / 1.0219 | -0.5592-17.1j / -1.093+11.1j | 4.98e-13 / 6.51e-14 | -0.648-0.0852j | 1.65 | OK_OFF_ONLY | FLIP:OFF_BAR50:y | OK |
| ISL_W2.9_q0.8_Li0.7_c0.6 | 2.900 / 6.162 / 0.544 | 2.54 / 0.274 | 1.1066 / 1.0365 | -0.1431-7.58j / -0.7594+4.57j | 1.45e-12 / 7.59e-14 | -0.601-0.112j | 1.6 | OK_OFF_ONLY | none | OK |

Alpha is the first-order tensor leakage at 25 deg (as the cell scan); beta is the self-consistent TM constant of the mean sheet; the A_perp of the island cell is the antisymmetric part of Y_yy at kt = beta (mS).

## Limits
Mesh 64 (not converged against 96 here), MoM only, no HFSS check; a 0.3 mm clearance is 2 mesh cells; the island is floating metal in a hole and its coupling to the BAR50 model (lumped sheet 1.00 x 0.35 mm) is the only DC/RF connection to the disc; no bias network is modelled. Stopped here: no board model, no HFSS package.
