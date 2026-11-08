# SUB6 board-island follow-up (stage: sub6-board-island, follow-up)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Code `scripts/research/sub6_board_island2.py` (+ `sub6_board_model4.py` unchanged); data `budget_best4.json`, `followup_points.csv` (72 designs x 4 variants), this file.

## 1. Why ON G_A < OFF G_A (best 200 / 220 mm points, N 135 and N 240)
| point (side N pl taper excl) | G_A OFF / ON | dG (ON-OFF) = dD + d10log(eta) | eta OFF / ON | A cone | B cone | other radiated | unconverted | est. sheet absorption (switched part) | G_A ON: beta shift off / on | switched dbeta, static dphi rms |
|---|---|---|---|---|---|---|---|---|---|---|
| 200 135 a 1.0 0 | 13.47 / 12.75 | -0.72 = +3.24 -3.96 | 0.196 / 0.079 | 0.088 / 0.067 | 0.064 / 0.0009 | 0.044 / 0.011 | 0.804 / 0.921 | 0.0099 (0.0077) / 0.0237 (0.0213) | 12.78 / 12.75 | -0.085, 13.6 deg |
| 200 240 b 0.5 1 | 16.70 / 17.32 | +0.62 = +2.73 -2.11 | 0.392 / 0.241 | 0.189 / 0.198 | 0.142 / 0.0033 | 0.061 / 0.040 | 0.608 / 0.759 | 0.0139 (0.0110) / 0.0281 (0.0250) | 17.41 / 17.32 | -0.070, 17.5 deg |
| 220 135 b 0.5 1 | 13.71 / 12.41 | -1.29 = +3.02 -4.31 | 0.157 / 0.058 | 0.069 / 0.048 | 0.054 / 0.0006 | 0.034 / 0.009 | 0.843 / 0.942 | 0.0087 (0.0066) / 0.0165 (0.0143) | 12.45 / 12.41 | -0.070, 9.6 deg |
| 220 240 b 0.8 0 | 17.17 / 17.04 | -0.13 = +2.67 -2.80 | 0.344 / 0.181 | 0.155 / 0.147 | 0.125 / 0.0029 | 0.064 / 0.031 | 0.656 / 0.819 | 0.0142 (0.0114) / 0.0356 (0.0324) | 17.16 / 17.04 | -0.087, 20.6 deg |

G = D_A + 10 log10(launch x eta): the ON deficit is the **conversion efficiency**, not the beta shift. In ON the island cells stop radiating, the static cells keep their (alpha-fixed) amplitudes, so eta falls from 0.16-0.39 to 0.06-0.24 and 76-94 % of the accepted power reaches the aperture edge unconverted; the directivity rises by +2.7 to +3.2 dB (one beam instead of two), which only partly compensates. At N 240 the two terms nearly cancel (dG +0.6 / -0.1 dB). **The switched-cell beta shift (-0.07 to -0.09; static-cell phase error rms 10-21 deg) costs only 0.03-0.12 dB of G_A,ON** (columns "beta shift off / on"). Radiated shares come from the far-field power split (cone <= 15 deg around each beam) times the bookkeeping eta, so they are model-level (eta from the alpha bookkeeping, shares from the pattern integral). Residual B in ON: 0.0006-0.0033 of the accepted power (0.2-5 % of the radiated power). Absorption is an ESTIMATE from Re(Y_mean) of each row (alpha_abs = Re Y / (d Im Y_tot/d beta), integrated along the rays with the leakage-only ray power); it is 1-3.6 % of the accepted power (ON diodes: about 1.4-3.2 %), is not subtracted from the gains, and the diode share in ON is about 2x the OFF value.

## 2. Beam-A design phase (averaged / optimised beta)
Since the beta shift is not the cause, the gain change is small: designing the static axes with the phase P_d = P_OFF + w (P_ON - P_OFF) and choosing w* from {0, .25, .5, .75, 1} to maximise min(G_A,OFF, G_A,ON) changes min(G_A) by +0.06 dB on average (max +0.13 dB over the 72 designs); w = 0.5 alone gives +0.03 dB and loses MEETS_ALL points (SLL). Not worth adopting. The real lever for G_A,ON is eta_on (more static conversion in ON, or fewer unused island cells), not the phase.

## 3. AR in OFF (rotation optimisation at several directions)
Rotation optimisation of ALL sites (90-degree flips, round-2 method, same costs) with the cross-pol constraint (LHCP/RHCP <= -20 dB) imposed at seven directions per beam (centre + 6 on a 5-degree ring) for A (OFF, ON) and B (OFF); candidates w in {0, w*}, ring weight in {0.1, 0.3, 1}, plus V2 itself; the design with the best gate margin is kept (V3).
| variant | points | MEETS_ALL (-3 dB AR) | mean OFF AR_A / AR_B (-3 dB beam) | mean ON AR_A | mean min(G_A) | mean B supp ON | mean SLL off / on | mean rotated |
|---|---|---|---|---|---|---|---|---|
| V0 round 'island' (w = 0, single-direction rotations) | 72 | 12 | 3.24 / 3.01 | 2.18 | 13.28 | 23.3 | -14.5 / -17.3 | 7 |
| V1 w = 0.5 (averaged beta) | 72 | 9 | 3.28 / 2.98 | 2.17 | 13.31 | 23.1 | -14.5 / -17.7 | 7 |
| V2 w* (best min G_A) | 72 | 13 | 3.30 / 2.97 | 2.16 | 13.33 | 22.5 | -14.6 / -17.8 | 7 |
| V3 = best of V2 and multi-direction rotations | 72 | 65 | 2.12 / 2.19 | 2.06 | 13.36 | 23.6 | -14.7 / -17.5 | 13 |

| side | N_sw | best V3 config (pl / taper / excl / w) | OFF G_A / G_B / bal | OFF AR_A / AR_B (-3 dB) | ON G_A / B supp / AR_A | SLL off / on | eta off / on | rotated | V3 gate (margin) | V0 best margin | MEETS_ALL V0 -> V3 (of 12) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 220 | 120 | a / 0.5 / yes / 0 | 12.8 / 12.0 / 0.73 | 1.93 / 1.87 | 10.4 / 15.8 / 2.05 | -14.3 / -17.6 | 0.14 / 0.04 | 15 | MEETS_ALL (0.95) | 0.23 | 3 -> 12 |
| 220 | 135 | a / 0.5 / yes / 0 | 13.6 / 12.9 / 0.66 | 2.46 / 2.39 | 11.5 / 15.9 / 2.33 | -13.7 / -18.0 | 0.16 / 0.05 | 7 | MEETS_ALL (0.54) | 0.39 | 3 -> 5 |
| 220 | 240 | a / 0.5 / yes / 0 | 16.7 / 16.5 / 0.21 | 1.85 / 1.81 | 16.2 / 15.5 / 2.06 | -15.9 / -17.8 | 0.32 / 0.15 | 20 | MEETS_ALL (0.94) | -0.17 | 0 -> 12 |
| 200 | 120 | a / 0.5 / yes / 0.75 | 12.8 / 11.6 / 1.18 | 2.02 / 2.07 | 11.9 / 25.8 / 1.99 | -14.4 / -16.6 | 0.17 / 0.06 | 10 | MEETS_ALL (0.93) | 0.06 | 4 -> 12 |
| 200 | 135 | b / 0.5 / yes / 0 | 13.7 / 11.8 / 1.90 | 1.67 / 1.93 | 11.8 / 18.1 / 1.81 | -15.3 / -18.8 | 0.18 / 0.06 | 14 | MEETS_ALL (1.07) | 0.17 | 2 -> 12 |
| 200 | 240 | a / 0.5 / no / 0 | 17.3 / 16.3 / 1.01 | 1.87 / 1.75 | 17.4 / 32.0 / 1.88 | -16.9 / -17.2 | 0.41 / 0.24 | 18 | MEETS_ALL (1.12) | -0.12 | 0 -> 12 |

The round-island rotation fixed the cross-pol only at the beam centre, leaving the -3 dB-beam AR of beam B at 2.2-3.1 dB; the multi-direction constraint brings OFF AR_A / AR_B to 1.7-2.5 dB at the cost of SLL (the optimiser does not see SLL, so the selection by gate margin keeps only designs that still meet -13 dB) and of 2x more rotated slots (mean 13, max 22; the cap is 200).

## Limits
MODEL_ONLY; mesh-64 island cells (switched table 9 rows); the selection of the best of several optimiser runs by the gate margin is a design search on the evaluation metric; rotated slots near each other with no coupling model; ON residual-B cancellation relies on the mean rho_d (-0.33); G_A numbers are referred to accepted power (launch 0.656 with the standard SMA pin of sub6-launcher). Stopped here.
