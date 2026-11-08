# SUB6 board (island cells) with HFSS-calibrated switched cells: DRY RUN (stage: sub6-board-island-hfsscal)

**MODEL_ONLY / DRY RUN. ON = MoM placeholder (no HFSS ON data yet), OFF calibrated from ONE phase per direction. Do not use for decisions.** Code `scripts/research/sub6_board_hfsscal.py` (imports `sub6_board_island2.py` and `sub6_board_model4.py` unchanged); data `hfss_data.json` (the user's HFSS numbers, quoted from the chat, not a pushed commit), `calibration.json`, `points_hfsscal.csv` (288 rows = 72 designs x V0-V3).

## Method
* Cell: the single island row `ISL_W2.9_q0.8_Li0.5_c0.6` (the row HFSS was run for; `base` and `cal` both use only this row, so they differ only by the calibration; the 9-row board-island numbers are not directly comparable).
* HFSS (ISL_vias2: vias + bias stubs, FR4 0.8, OFF, 85 deg): along 5.963 GHz Q 250, across 5.566 GHz Q 189 (surface-wave modes); not used: via resonance 4.54 / 4.70 GHz (Q about 125), air mode 6.023 GHz.
* Per state and direction the MoM dispersion f_M(phase) is corrected, f_cal = f_M + a (+ b (phi - phi_ref) once two phases exist); beta/k0 at 5.8 GHz from f_cal; v_g from its slope; alpha_x = omega / (2 Q v_g) = f (d beta/d f) / (2 Q). The MoM leakage alpha of the row is kept (MoM-only); alpha_x is absorbed, not radiated.
* On the board: beta(theta') = beta_along cos^2 theta' + beta_across sin^2 theta' with theta' = ray angle - slot axis of the hologram (3 iterations, axes before the 90-degree flips; the flips are not fed back); alpha_x likewise. ON: MoM beta and alpha, alpha_x = 0.
| direction | HFSS points (phase, f, Q) | MoM f at that phase | shift a (GHz) | beta/k0 at 5.8 GHz: MoM -> calibrated | d phi / d f (deg/GHz) | alpha_x = omega/(2 Q v_g) (Np/m) |
|---|---|---|---|---|---|---|
| along | 85 deg, 5.963 GHz, Q 250 | 5.784 | +0.179 | 1.0556 -> 1.0168 | 17.3 | 0.301 |
| across | 85 deg, 5.566 GHz, Q 189 | 5.519 | +0.047 | 1.1671 -> 1.1461 | 36.8 | 0.849 |

## Result (V3, 200 mm, 12 designs per N: a/b x taper 0.5/0.8/1.0 x ring exclusion off/on)
| N_sw | mode | MEETS_ALL of 12 | best design (pl / taper / excl) | gate (margin) | OFF G_A / G_B / bal | OFF AR_A / AR_B | ON G_A / B supp / AR_A | SLL off / on | eta off / on | absorbed OFF | switched beta OFF / ON (mean) | alpha_x (Np/m) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 120 | base (MoM, same single row) | 4 | a / 0.5 / 0 | MEETS_ALL (0.04) | 12.0 / 11.0 / 1.03 | 2.22 / 2.20 | 9.8 / 30.5 / 1.74 | -13.0 / -19.1 | 0.150 / 0.036 | 0.000 | 1.104 / 1.005 | 0.00 |
| 120 | cal (ON = MoM placeholder) | 0 | a / 0.5 / 0 | MISSES_SLL (-0.13) | 12.0 / 10.8 / 1.17 | 2.07 / 2.50 | 9.7 / 28.1 / 1.75 | -12.9 / -17.8 | 0.148 / 0.036 | 0.027 | 1.085 / 1.005 | 0.59 |
| 135 | base (MoM, same single row) | 3 | a / 1 / 1 | MEETS_ALL (0.97) | 12.3 / 11.6 / 0.73 | 1.50 / 2.03 | 10.2 / 28.5 / 1.82 | -14.0 / -18.4 | 0.170 / 0.043 | 0.000 | 1.104 / 1.005 | 0.00 |
| 135 | cal (ON = MoM placeholder) | 3 | a / 1 / 1 | MEETS_ALL (0.76) | 12.2 / 11.7 / 0.48 | 1.86 / 2.24 | 10.3 / 26.0 / 1.84 | -14.2 / -16.9 | 0.167 / 0.043 | 0.030 | 1.084 / 1.005 | 0.58 |
| 240 | base (MoM, same single row) | 11 | a / 0.5 / 0 | MEETS_ALL (1.13) | 16.0 / 15.5 / 0.56 | 1.79 / 1.86 | 15.6 / 27.0 / 1.87 | -15.0 / -19.8 | 0.345 / 0.150 | 0.000 | 1.104 / 1.005 | 0.00 |
| 240 | cal (ON = MoM placeholder) | 12 | a / 0.5 / 1 | MEETS_ALL (1.12) | 15.9 / 15.1 / 0.80 | 1.88 / 1.75 | 15.3 / 27.4 / 1.81 | -15.2 / -16.3 | 0.330 / 0.149 | 0.046 | 1.084 / 1.005 | 0.58 |

Paired change (cal - base):
| N_sw | paired designs | mean change OFF G_A (dB) | mean change ON G_A (dB) | mean change of the gate margin (dB) | mean eta_off ratio cal/base | designs MEETS_ALL base -> cal |
|---|---|---|---|---|---|---|
| 120 | 12 | +0.03 | -0.03 | -0.15 | 0.984 | 4 -> 0 |
| 135 | 12 | -0.03 | +0.03 | -0.14 | 0.983 | 3 -> 3 |
| 240 | 12 | -0.05 | -0.03 | +0.05 | 0.974 | 11 -> 12 |

## Reading
The OFF calibration lowers the switched-cell beta (mean 1.104 -> 1.085) and adds 0.58 Np/m of absorbed loss (2.7-4.6 % of the launched power); it changes the gains by <= 0.05 dB and the gate margin by -0.15 ... +0.05 dB on average, but N 120 loses all 4 passing designs (SLL misses by 0.1-0.2 dB), N 135 keeps 3, N 240 goes 11 -> 12 of 12. The ON-state numbers are the MoM ones and the ON calibration (beta shift ON-OFF, ON loss) will matter more: the OFF-ON beta change of the switched cells moves from 0.10 (MoM) to 0.08 now but is untested.

## Limits
One HFSS point per direction (shift only, MoM slope assumed for v_g and for the phase-to-frequency map); Q includes every loss of the HFSS cell (diode R, FR4, copper, the bottom radiation boundary through the anti-pads) and is treated as extra absorption of the switched cells only; the static cells are not calibrated; the theta' weighting (cos^2) is a model choice; the flips are not fed back into beta; HFSS v2 has the 2.0 x 0.3 mm bias stubs. Waiting for: OFF 95 deg, ON 80 / 85 / 95 deg.
