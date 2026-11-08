# SUB6 board model, round 2: polarisation-aware synthesis (stage: sub6-board-model-r2)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Same cells (`results/sub6_slotdisc/rows_RT5880.jsonl`, mesh-64 MoM, no HFSS check), same architecture and gates as `BOARD_MODEL.md` (round 1, untouched). Code `scripts/research/sub6_board_model2.py` (+ `sub6_board_model.py`, `sub6_board_req.py` unchanged); data `points_r2.csv` (120 points), `tolerance_C.csv`, `layout_<side>mm_N<n>_<a|b>_t0.8.csv` (taper 0.8 with ring exclusion; columns include `rotated90`), `provenance_run2.json`. Launch fraction 0.6263. 15 s for 120 points on 64 workers.

## What changed
1. **Full tensor in the far field.** Cell current J = |E| [ (Delta/2) e(2a - phi_R) + A_eff e(phi_R - 90 deg) ], A_eff = A_perp sin^2(theta') + A_par cos^2(theta'), theta' = phi_R - a. The Delta term is the engine's (amplitude sqrt(2 alpha W / rho)); the A term is added with the complex ratio A_eff/(Delta/2) of the cell state (OFF_BAR50 / ON_BAR50; NOCUT rows have A = 0, so only switched cells carry it; mean |ratio| OFF about 5, ON about 0.05, but the OFF amplitude itself is tiny). A 90 deg rotation of a slot flips the Delta term and swaps sin^2 <-> cos^2. **Assumption (mine): the table alpha is the Delta-term leakage and the A term is added on top.** The A term alone has wavevector beta (bound), so it matters only through the discrete sampling of the sites.
2. **Joint axis choice (joint_opt2 style).** Greedy single 90-degree rotations (psi + 90 deg, same in OFF and ON) over ALL live sites, static and switched, with the corrected RHCP index (`copol_component(...,'rhcp')`). Cost: maximise |R_A|^2 in OFF and ON and |R_B|^2 in OFF, penalise |R_B|^2 in ON; constraints: LHCP/RHCP <= -20 dB at beam A (both states) and beam B (OFF), co-pol floors 0.5 dB (A) / 1.0 dB (B) below the unrotated sums, B-in-ON no higher than max(unrotated, 17 dB below B OFF); at most 200 rotations. Rotated slots per point: 1-7 % of the cells (CP_STEP3 risk 3: mutual coupling of rotated cells not modelled).
3. **ON SLL.** Taper c in {0.5, 0.8, 1.0} (g = 1 - c (rho/rho_edge)^2 on the march target; round 1 used 0.5); switched cells excluded from lattice rings 1-2 around the feed (ring = max(|i|,|j|) <= 2) and from the outer ring (cells with a missing 8-neighbour) when "ring exclusion" is on. Both placements are restricted to the allowed sites (placement b = the stage-1b optimiser with the same swap rule).
4. **Sweep.** N_sw 180/200/220/240/260 at 200 mm (936 cells) and 220 mm, a/b, taper x3, exclusion x2 = 120 points. 220 mm is not scaled by the cell count here (N_sw as given).

## Gates (balance <= 3 dB, AR <= 3 dB, B suppression >= 12 dB, SLL <= -13 dB; SLL = worst of both states)
**Peak-AR gate (the job's gate, `metrics`: AR at the RHCP peak within 6 deg of the beam): MEETS_ALL 87 of 120. With AR taken as the worst over the -3 dB RHCP contour of the beam (stricter, my addition): MEETS_ALL 34 of 120.** The joint rotation puts the cross-pol at the beam-A direction at its floor (-20 dB gives AR 1.74 dB), so peak-AR margins are all about 1.26 dB; the -3 dB-beam AR is 0.2-1.2 dB worse. Best point by the stricter margin: 220 mm, N_sw 240, placement b, taper 0.5, ring exclusion yes.

Best configuration per (side, N_sw), chosen by the stricter margin:

| side | N_sw (f) | best config pl / taper / ring-excl | OFF G_A / G_B / bal | OFF AR A / B (peak; -3 dB beam) | ON G_A / B supp | ON AR_A (peak; -3 dB beam) | SLL off / on | eta off / on | rotated slots | gate peak-AR / -3dB-AR (margin) |
|---|---|---|---|---|---|---|---|---|---|---|
| 220 | 180 (0.19) | b / 1 / yes | 11.8 / 12.0 / 0.23 | 1.74 / 1.69; 2.01 / 2.30 | 16.4 / 17.3 | 1.73; 2.48 | -19.4 / -14.0 | 0.11 / 0.28 | 24 (16 sw) | MEETS_ALL / MEETS_ALL (0.52) |
| 220 | 200 (0.21) | a / 1 / yes | 12.8 / 13.1 / 0.22 | 1.74 / 1.40; 1.78 / 2.38 | 17.3 / 17.8 | 1.68; 2.67 | -19.7 / -14.7 | 0.14 / 0.33 | 29 (23 sw) | MEETS_ALL / MEETS_ALL (0.33) |
| 220 | 220 (0.23) | b / 0.5 / yes | 14.3 / 14.6 / 0.23 | 1.61 / 1.65; 1.87 / 2.32 | 18.2 / 16.9 | 1.74; 2.26 | -17.6 / -14.1 | 0.18 / 0.37 | 23 (13 sw) | MEETS_ALL / MEETS_ALL (0.68) |
| 220 | 240 (0.26) | b / 0.5 / yes | 14.9 / 15.2 / 0.24 | 1.66 / 1.56; 1.92 / 2.30 | 18.7 / 17.1 | 1.67; 2.20 | -17.7 / -14.6 | 0.21 / 0.41 | 29 (16 sw) | MEETS_ALL / MEETS_ALL (0.70) |
| 220 | 260 (0.28) | b / 1 / yes | 14.6 / 14.8 / 0.19 | 1.74 / 1.65; 2.09 / 2.65 | 18.6 / 17.3 | 1.50; 1.90 | -17.3 / -14.9 | 0.22 / 0.43 | 37 (29 sw) | MEETS_ALL / MEETS_ALL (0.35) |
| 200 | 180 (0.24) | b / 1 / yes | 12.2 / 12.3 / 0.11 | 1.77 / 1.52; 2.22 / 2.91 | 16.3 / 16.9 | 1.70; 3.26 | -14.3 / -14.8 | 0.16 / 0.34 | 21 (15 sw) | MEETS_ALL / MISSES_AR (-0.26) |
| 200 | 200 (0.26) | a / 0.8 / no | 13.3 / 13.5 / 0.13 | 1.62 / 0.84; 2.10 / 2.10 | 17.4 / 17.1 | 1.69; 3.20 | -14.9 / -14.8 | 0.19 / 0.39 | 28 (20 sw) | MEETS_ALL / MISSES_AR (-0.20) |
| 200 | 220 (0.29) | a / 1 / no | 14.1 / 13.8 / 0.26 | 1.81 / 1.13; 2.12 / 2.75 | 17.6 / 16.6 | 1.89; 3.07 | -14.9 / -15.2 | 0.23 / 0.43 | 20 (15 sw) | MEETS_ALL / MISSES_AR (-0.07) |
| 200 | 240 (0.32) | a / 0.8 / no | 15.0 / 15.1 / 0.14 | 1.77 / 1.23; 2.33 / 2.76 | 18.4 / 16.2 | 1.82; 2.92 | -15.7 / -13.7 | 0.29 / 0.50 | 17 (14 sw) | MEETS_ALL / MEETS_ALL (0.08) |
| 200 | 260 (0.34) | b / 0.5 / no | 15.7 / 15.9 / 0.15 | 1.78 / 1.06; 2.24 / 2.57 | 18.9 / 17.0 | 1.78; 2.59 | -16.1 / -13.6 | 0.34 / 0.54 | 27 (20 sw) | MEETS_ALL / MEETS_ALL (0.41) |

Gain is referred to accepted power x model conversion eta, so absolute G_A are model numbers; ON G_A is above OFF because the switched cells add conversion in ON (eta_on > eta_off).

## What each lever did (means over the subset)
| subset | points | ON AR_A: no A term / A term, unrotated / after joint rotation | SLL_on: same three | SLL_off after | B supp ON: no A / unrotated / after | MEETS_ALL peak-AR / -3 dB-AR |
|---|---|---|---|---|---|---|
| all | 120 | 3.49 / 3.48 / 1.75 | -12.0 / -12.0 / -13.8 | -16.5 | 22.4 / 22.7 / 16.8 | 87 / 34 |
| taper 0.5 | 40 | 3.37 / 3.37 / 1.79 | -12.2 / -12.1 / -13.5 | -16.6 | 22.6 / 22.9 / 16.8 | 30 / 9 |
| taper 0.8 | 40 | 3.41 / 3.41 / 1.76 | -12.0 / -12.1 / -13.6 | -16.4 | 22.5 / 22.9 / 16.8 | 28 / 10 |
| taper 1.0 | 40 | 3.68 / 3.67 / 1.71 | -11.8 / -11.8 / -14.3 | -16.4 | 22.1 / 22.2 / 16.8 | 29 / 15 |
| ring exclusion off | 60 | 4.29 / 4.26 / 1.81 | -12.6 / -12.8 / -13.5 | -16.8 | 23.3 / 24.4 / 16.3 | 43 / 12 |
| ring exclusion on | 60 | 2.69 / 2.70 / 1.69 | -11.4 / -11.2 / -14.1 | -16.1 | 21.6 / 21.0 / 17.3 | 44 / 22 |
| 200 mm | 60 | 2.80 / 2.80 / 1.76 | -11.4 / -11.3 / -13.9 | -14.8 | 23.2 / 23.1 / 17.0 | 43 / 5 |
| 220 mm | 60 | 4.18 / 4.16 / 1.75 | -12.6 / -12.7 / -13.7 | -18.1 | 21.6 / 22.2 / 16.6 | 44 / 29 |
| placement a | 60 | 3.50 / 3.49 / 1.77 | -12.0 / -12.1 / -13.8 | -16.4 | 22.4 / 22.2 / 17.0 | 42 / 12 |
| placement b | 60 | 3.48 / 3.48 / 1.74 | -12.0 / -11.9 / -13.8 | -16.6 | 22.5 / 23.2 / 16.6 | 45 / 22 |

"no A term" = the round-1 current model; "A term, unrotated" = tensor with A but the designed axes; "after" = after the joint rotations. B suppression after rotation is held near the 17 dB cap (unrotated 19-29 dB): the rotations trade suppression for AR and SLL, the gate (12 dB) is still met.

## Tolerance of the diode capacitance at the best point (switched table restricted to the 5 rows of `results/sub6_slotdisc_ls_tol`)
| diode C (pF) | switched OFF alpha mean (Np/m) | OFF G_A / G_B | balance | B contrast OFF-ON (dB) | ON B suppression | OFF AR A / B (peak; -3 dB) | ON AR_A (peak; -3 dB) | SLL off / on | gate |
|---|---|---|---|---|---|---|---|---|---|
| 0.12 | 0.0138 | 15.6 / 14.4 | 1.17 | 12.8 | 17.1 | 1.83 / 1.43; 2.06 / 2.45 | 1.67; 2.20 | -18.3 / -14.6 | MEETS_ALL |
| 0.15 | 0.0031 | 14.9 / 15.2 | 0.24 | 13.6 | 17.1 | 1.66 / 1.56; 1.92 / 2.30 | 1.67; 2.20 | -17.7 / -14.6 | MEETS_ALL |
| 0.20 | 0.1178 | 15.4 / 14.5 | 0.94 | 12.9 | 17.1 | 1.58 / 1.52; 1.88 / 2.42 | 1.66; 2.20 | -16.8 / -14.6 | MEETS_ALL |

Source data `input_ls_tol.csv` (copied from github/compute/sub6-slotdisc-ls-tol 6d9c0ba). Design: the best point rebuilt with only the 5 ls_tol rows as switched cells (C = 0.15 pF nominal), then OFF alpha, beta and Delta of those rows replaced by the C = 0.12 / 0.20 values, sites, slot axes and rotations unchanged; A_perp kept at its nominal value (not in ls_tol). C acts on the OFF state only, so ON B suppression and ON AR do not change; the C effect is the OFF leakage of the switched cells (OFF alpha mean 0.014 at 0.12 pF and 0.118 at 0.20 pF against 0.003 nominal; at 0.20 pF one row exceeds the 0.1 Np/m filter), seen as OFF balance / OFF B and the B contrast OFF-ON.

## Limits
Cells are mesh-64 MoM, not HFSS-checked (2 of 80 rows have an HFSS package). The A-term model (see 1), the AR/XPD target at a single direction, and rotated slots next to unrotated ones without coupling are my modelling choices. Mixed gap / W neighbours: coupling not modelled. Stopped here: no cell choice, no board design, no HFSS package.
