# SUB6 board model, round 3: diode cap 150 (N_switched <= 135), same method as round 2 (stage: sub6-board-model-r3)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Same cells (`results/sub6_slotdisc/rows_RT5880.jsonl`, mesh-64 MoM, no HFSS check), same architecture and gates as `BOARD_MODEL.md` (round 1, untouched). Code `scripts/research/sub6_board_model2.py` (+ `sub6_board_model.py`, `sub6_board_req.py` unchanged); data `points_r2.csv` (108 points), `tolerance_C.csv`, `layout_<side>mm_N<n>_<a|b>_t0.8.csv` (taper 0.8 with ring exclusion; columns include `rotated90`), `provenance_run2.json`. Launch fraction 0.6263. 13 s for 108 points on 64 workers.

## What changed
1. **Full tensor in the far field.** Cell current J = |E| [ (Delta/2) e(2a - phi_R) + A_eff e(phi_R - 90 deg) ], A_eff = A_perp sin^2(theta') + A_par cos^2(theta'), theta' = phi_R - a. The Delta term is the engine's (amplitude sqrt(2 alpha W / rho)); the A term is added with the complex ratio A_eff/(Delta/2) of the cell state (OFF_BAR50 / ON_BAR50; NOCUT rows have A = 0, so only switched cells carry it; mean |ratio| OFF about 5, ON about 0.05, but the OFF amplitude itself is tiny). A 90 deg rotation of a slot flips the Delta term and swaps sin^2 <-> cos^2. **Assumption (mine): the table alpha is the Delta-term leakage and the A term is added on top.** The A term alone has wavevector beta (bound), so it matters only through the discrete sampling of the sites.
2. **Joint axis choice (joint_opt2 style).** Greedy single 90-degree rotations (psi + 90 deg, same in OFF and ON) over ALL live sites, static and switched, with the corrected RHCP index (`copol_component(...,'rhcp')`). Cost: maximise |R_A|^2 in OFF and ON and |R_B|^2 in OFF, penalise |R_B|^2 in ON; constraints: LHCP/RHCP <= -20 dB at beam A (both states) and beam B (OFF), co-pol floors 0.5 dB (A) / 1.0 dB (B) below the unrotated sums, B-in-ON no higher than max(unrotated, 17 dB below B OFF); at most 200 rotations. Rotated slots per point: 1-4 % of the cells (CP_STEP3 risk 3: mutual coupling of rotated cells not modelled).
3. **ON SLL.** Taper c in {0.5, 0.8, 1.0} (g = 1 - c (rho/rho_edge)^2 on the march target; round 1 used 0.5); switched cells excluded from lattice rings 1-2 around the feed (ring = max(|i|,|j|) <= 2) and from the outer ring (cells with a missing 8-neighbour) when "ring exclusion" is on. Both placements are restricted to the allowed sites (placement b = the stage-1b optimiser with the same swap rule).
4. **Sweep.** N_sw 110/120/135 at 180, 200 and 220 mm, a/b, taper x3, ring exclusion off/on = 108 points. The diode cap (150) allows at most 135 switched cells.

## Gates (balance <= 3 dB, AR <= 3 dB, B suppression >= 12 dB, SLL <= -13 dB; SLL = worst of both states)
**Peak-AR gate (the job's gate, `metrics`: AR at the RHCP peak within 6 deg of the beam): MEETS_ALL 11 of 108. With AR taken as the worst over the -3 dB RHCP contour of the beam (stricter, my addition): MEETS_ALL 0 of 108.** Misses are mostly SLL (ON state, -12.3 to -13.3 dB at the best configurations) and, on the stricter AR, ON AR over the -3 dB beam (up to 4.2 dB). The tolerance study (below) uses the point with the best stricter margin, which still misses SLL by a hair (see table), and its 5-row switched table misses SLL by about 1.6 dB. Best point by the stricter margin: 200 mm, N_sw 135, placement a, taper 1, ring exclusion no.

Best configuration per (side, N_sw), chosen by the stricter margin:

| side | N_sw (f) | best config pl / taper / ring-excl | OFF G_A / G_B / bal | OFF AR A / B (peak; -3 dB beam) | ON G_A / B supp | ON AR_A (peak; -3 dB beam) | SLL off / on | eta off / on | rotated slots | gate peak-AR / -3dB-AR (margin) |
|---|---|---|---|---|---|---|---|---|---|---|
| 220 | 110 (0.12) | b / 1 / yes | 8.1 / 8.3 / 0.17 | 1.74 / 1.58; 2.38 / 1.96 | 13.2 / 16.8 | 1.64; 1.99 | -19.5 / -12.3 | 0.04 / 0.16 | 11 (5 sw) | MISSES_SLL / MISSES_SLL (-0.73) |
| 220 | 120 (0.13) | b / 0.5 / yes | 8.9 / 9.0 / 0.18 | 1.74 / 1.69; 2.36 / 2.23 | 14.0 / 17.0 | 1.85; 3.01 | -18.4 / -12.5 | 0.05 / 0.18 | 10 (6 sw) | MISSES_SLL / MISSES_AR+SLL (-0.47) |
| 220 | 135 (0.14) | a / 1 / yes | 9.9 / 8.3 / 1.62 | 1.73 / 1.40; 2.09 / 2.45 | 14.4 / 18.6 | 1.89; 3.03 | -18.7 / -12.6 | 0.06 / 0.19 | 16 (11 sw) | MISSES_SLL / MISSES_AR+SLL (-0.41) |
| 200 | 110 (0.15) | b / 1 / yes | 8.2 / 6.7 / 1.43 | 1.74 / 1.58; 2.33 / 2.99 | 12.5 / 17.4 | 1.62; 4.21 | -15.4 / -13.3 | 0.05 / 0.17 | 11 (8 sw) | MEETS_ALL / MISSES_AR (-1.21) |
| 200 | 120 (0.16) | a / 1 / no | 9.0 / 9.1 / 0.07 | 1.71 / 1.20; 2.06 / 2.14 | 14.1 / 16.2 | 1.68; 3.04 | -14.2 / -12.8 | 0.07 / 0.21 | 9 (7 sw) | MISSES_SLL / MISSES_AR+SLL (-0.17) |
| 200 | 135 (0.18) | a / 1 / no | 9.8 / 9.8 / 0.04 | 1.72 / 1.25; 2.09 / 2.17 | 14.7 / 15.9 | 1.37; 2.22 | -14.6 / -13.0 | 0.08 / 0.24 | 15 (11 sw) | MISSES_SLL / MISSES_SLL (-0.00) |
| 180 | 110 (0.18) | b / 0.5 / yes | 9.0 / 9.2 / 0.23 | 1.64 / 1.25; 1.73 / 2.40 | 13.3 / 15.6 | 1.92; 3.19 | -12.6 / -12.6 | 0.08 / 0.21 | 18 (9 sw) | MISSES_SLL / MISSES_AR+SLL (-0.40) |
| 180 | 120 (0.19) | b / 0.5 / yes | 9.7 / 9.9 / 0.15 | 1.70 / 1.26; 1.80 / 2.47 | 13.9 / 15.5 | 1.87; 3.01 | -13.4 / -12.7 | 0.09 / 0.24 | 19 (11 sw) | MISSES_SLL / MISSES_AR+SLL (-0.26) |
| 180 | 135 (0.21) | b / 1 / no | 10.5 / 10.5 / 0.02 | 1.69 / 1.40; 2.00 / 2.70 | 14.7 / 16.3 | 1.94; 3.39 | -15.0 / -12.9 | 0.11 / 0.27 | 19 (12 sw) | MISSES_SLL / MISSES_AR+SLL (-0.39) |

Gain is referred to accepted power x model conversion eta, so absolute G_A are model numbers; ON G_A is above OFF because the switched cells add conversion in ON (eta_on > eta_off).

## What each lever did (means over the subset)
| subset | points | ON AR_A: no A term / A term, unrotated / after joint rotation | SLL_on: same three | SLL_off after | B supp ON: no A / unrotated / after | MEETS_ALL peak-AR / -3 dB-AR |
|---|---|---|---|---|---|---|
| all | 108 | 4.70 / 4.65 / 1.79 | -10.7 / -10.7 / -11.7 | -15.2 | 22.2 / 22.8 / 16.3 | 11 / 0 |
| taper 0.5 | 36 | 4.64 / 4.60 / 1.78 | -10.7 / -10.7 / -11.6 | -15.2 | 22.1 / 22.9 / 16.2 | 4 / 0 |
| taper 0.8 | 36 | 4.64 / 4.60 / 1.78 | -10.7 / -10.7 / -11.6 | -15.2 | 22.1 / 22.9 / 16.2 | 4 / 0 |
| taper 1.0 | 36 | 4.80 / 4.76 / 1.81 | -10.7 / -10.7 / -11.9 | -15.3 | 22.4 / 22.6 / 16.5 | 3 / 0 |
| ring exclusion off | 54 | 5.52 / 5.48 / 1.83 | -11.2 / -11.1 / -11.4 | -15.0 | 22.6 / 23.9 / 15.9 | 3 / 0 |
| ring exclusion on | 54 | 3.87 / 3.83 / 1.75 | -10.2 / -10.3 / -12.1 | -15.5 | 21.8 / 21.6 / 16.7 | 8 / 0 |
| 200 mm | 36 | 4.54 / 4.50 / 1.73 | -10.6 / -10.6 / -12.1 | -14.5 | 23.6 / 23.6 / 16.9 | 8 / 0 |
| 220 mm | 36 | 4.73 / 4.70 / 1.88 | -11.2 / -11.2 / -11.5 | -18.0 | 21.1 / 22.3 / 16.7 | 0 / 0 |
| placement a | 54 | 4.80 / 4.76 / 1.81 | -10.8 / -10.8 / -11.8 | -15.1 | 22.3 / 22.4 / 16.5 | 8 / 0 |
| placement b | 54 | 4.59 / 4.54 / 1.78 | -10.6 / -10.6 / -11.6 | -15.4 | 22.0 / 23.1 / 16.2 | 3 / 0 |

"no A term" = the round-1 current model; "A term, unrotated" = tensor with A but the designed axes; "after" = after the joint rotations. B suppression after rotation is held near the 17 dB cap (unrotated 19-29 dB): the rotations trade suppression for AR and SLL, the gate (12 dB) is still met.

## Tolerance of the diode capacitance at the best point (switched table restricted to the 5 rows of `results/sub6_slotdisc_ls_tol`)
| diode C (pF) | switched OFF alpha mean (Np/m) | OFF G_A / G_B | balance | B contrast OFF-ON (dB) | ON B suppression | OFF AR A / B (peak; -3 dB) | ON AR_A (peak; -3 dB) | SLL off / on | gate |
|---|---|---|---|---|---|---|---|---|---|
| 0.12 | 0.0205 | 10.4 / 9.1 | 1.25 | 9.4 | 14.7 | 1.78 / 1.50; 2.30 / 2.77 | 1.78; 3.05 | -15.8 / -11.4 | MISSES_SLL |
| 0.15 | 0.0043 | 9.8 / 9.8 | 0.02 | 10.1 | 14.7 | 1.72 / 1.38; 2.00 / 2.38 | 1.79; 3.07 | -14.5 / -11.4 | MISSES_SLL |
| 0.20 | 0.1020 | 10.2 / 9.4 | 0.80 | 9.7 | 14.7 | 1.65 / 0.93; 2.20 / 2.20 | 1.80; 3.09 | -13.7 / -11.4 | MISSES_SLL |

Source data `input_ls_tol.csv` (copied from github/compute/sub6-slotdisc-ls-tol 6d9c0ba). Design: the best point rebuilt with only the 5 ls_tol rows as switched cells (C = 0.15 pF nominal), then OFF alpha, beta and Delta of those rows replaced by the C = 0.12 / 0.20 values, sites, slot axes and rotations unchanged; A_perp kept at its nominal value (not in ls_tol). C acts on the OFF state only, so ON B suppression and ON AR do not change; the C effect is the OFF leakage of the switched cells (OFF alpha mean 0.021 at 0.12 pF and 0.102 at 0.20 pF against 0.004 nominal; at 0.20 pF one row exceeds the 0.1 Np/m filter), seen as OFF balance / OFF B and the B contrast OFF-ON.

## Limits
Cells are mesh-64 MoM, not HFSS-checked (2 of 80 rows have an HFSS package). The A-term model (see 1), the AR/XPD target at a single direction, and rotated slots next to unrotated ones without coupling are my modelling choices. Mixed gap / W neighbours: coupling not modelled. Stopped here: no cell choice, no board design, no HFSS package.
