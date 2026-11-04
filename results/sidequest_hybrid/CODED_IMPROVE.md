# Coded HMS K = 2: efficiency items I1-I7 at equal or lower diode / line count — MODEL_ONLY
Chosen cell: crossed split disc, rim 0.200 mm, two diodes per cell (x, y), measured MA4AGP907. Board: hexagon of side/circumradius 150 mm (58.5e3 mm^2, Dmax = 4 pi A / lambda^2 = 29.6 dBi), centre probe, states A RHCP (25, 0), B LHCP (25, 180).
Code: `scripts/research/sidequest_coded_improve.py` (model), `sidequest_coded_lib.py` (MoM sub-library), `sidequest_coded_improve_run.py` (I1-I6, alone and cumulative), `sidequest_coded_improve_i7.py` (I7), `sidequest_launcher_opt.py` (I4). Data `coded_lib.json`, `coded_improve*.json`, `launcher_opt.json`. Tests `tests/test_sidequest_coded_improve.py`.

## Baseline B0 (rebased)
The 18.2 dBi of the reference used the mesh-64 cell (alpha 0.89 Np/m). With the mesh-converged estimate (alpha_beam 0.70 Np/m, alpha_abs 0.095 Np/m with half the diodes ON, beta/k0 1.07; `CODED_RIM.md`) and the same least-squares design (diodes <= 766, lambda tuned) the baseline is **17.22 / 17.22 dBi, SLL -13.7 dB, AR 0.84 dB, XPD 26.3 dB, aperture efficiency 5.8 %, SW->beam 16.4 %, edge power 82.5 %, 73.2 % static, 766 diodes, 1 line + inverter per axis**. Ohmic loss is now counted only for ON diodes (static cells are metal): ON weight 0.33.
`B0L` = same but with the static cells limited to the MoM library (open slit m = +0.59, shorted slit m = -1.46) instead of the idealised +-1: 17.46 dBi (the library's short is a bigger swing than the diode ON state).

## Model (what each item can and cannot do)
Per element value m (m = +1 diode OFF, -1 diode ON of the baseline cell; library Y_xx mapped as m = (Y - Y_ref)/H_base, Y_ref = mean of baseline OFF/ON, H_base its half swing 1.62 mS). Cell current J = sum_e m_e (pi/2) amp (u_e . rho_hat) u_e e^(-j beta rho); alpha_i = A0 sum_e (m_e u_e.rho_hat)^2, power bookkeeping along 12 azimuth sectors (6.4 mm radial bins); per-cell least squares to the complex target g amp e^(-j eps_k) over discrete option sets (static levels / switched pairs), + lambda per diode, lambda re-tuned so that diodes <= 766; g (target scale) scanned {1, 1.5, 2, 2.7}, best gain with SLL <= -12.7 dB.
**Library (MoM, mesh 64; 23 variants):** switched (diode OFF/ON; slit half-length 2.6/2.8 mm x diode position 0.8/1.2/1.6/2.0 mm): Y_xx OFF 7.6-10.0j, ON 5.9-6.4j mS, ON conductance 0.05-0.15 mS (position 0.8: swing 2.37 vs 2.00 in m units and 36 % lower ON conductance than the baseline position 1.6); static: no slit 5.18j, shorted slit 5.2-5.5j, open slit 5.3-8.8j mS.

## Results: each item alone on B0
| item | delta vs B0 |
|---|---|
| I1 analog static cells (library-limited range -1.66..+0.59) | +0.23 / +0.23 dB, SLL +0.1 dB, AR -0.14 dB, ap_eff +0.32 %pt |
| I2 pre-phased switched cells (ls 2.6/2.8, pos 1.6) | -0.01 / -0.01 dB, SLL +0.5 dB, AR +0.02 dB, ap_eff -0.02 %pt |
| I3 radial taper (uniform illumination) | -2.39 / -2.39 dB, SLL +0.8 dB, AR +0.03 dB, ap_eff -2.46 %pt |
| I4 launcher (0.80, S11 NOT met) | +1.06 / +1.06 dB, SLL +0.0 dB, AR +0.00 dB, ap_eff +1.60 %pt |
| I5 diode position + loss-aware assignment | +0.56 / +0.56 dB, SLL -0.2 dB, AR +0.15 dB, ap_eff +0.79 %pt |
| I6 edge rings | -0.16 / -0.16 dB, SLL +0.4 dB, AR +0.02 dB, ap_eff -0.21 %pt |

## Results: cumulative (adopted: I1, I2, I4, I5; I3 and I6 lose gain and are not adopted)
| step | delta vs previous step |
|---|---|
| +I1 | +0.23 / +0.23 dB, SLL +0.1 dB, AR -0.14 dB, ap_eff +0.32 %pt |
| +I2 | +0.02 / +0.02 dB, SLL +0.1 dB, AR -0.03 dB, ap_eff +0.02 %pt |
| +I3 | +0.00 / +0.00 dB, SLL +0.0 dB, AR +0.00 dB, ap_eff +0.00 %pt |
| +I4 | +1.06 / +1.06 dB, SLL +0.0 dB, AR +0.00 dB, ap_eff +1.70 %pt |
| +I5 | +0.85 / +0.85 dB, SLL -1.1 dB, AR +0.20 dB, ap_eff +1.70 %pt |
| +I6 | -0.25 / -0.25 dB, SLL -0.5 dB, AR -0.07 dB, ap_eff -0.53 %pt |
| cumulative (adopted I1+I2+I4+I5) vs B0 | +2.16 / +2.16 dB, SLL -0.9 dB, AR +0.03 dB, ap_eff +3.74 %pt |

## Full table
| run | gain A/B [dBi] | SLL [dB] | AR [dB] | ap. eff. [%] | edge [%] | static [%] | SW->beam [%] | diodes | lines |
|---|---|---|---|---|---|---|---|---|---|
| B0 | [17.22, 17.22] | [-13.7, -13.7] | [0.84, 0.84] | 5.81 | 82.5 | 73.2 | 16.4 | 766 | 1 |
| B0L | [17.46, 17.46] | [-13.7, -13.7] | [0.77, 0.77] | 6.14 | 80.8 | 73.5 | 18.1 | 758 | 1 |
| I1_alone | [17.45, 17.45] | [-13.6, -13.6] | [0.7, 0.71] | 6.13 | 82.7 | 73.7 | 16.3 | 752 | 1 |
| I2_alone | [17.21, 17.21] | [-13.2, -13.2] | [0.86, 0.86] | 5.79 | 82.8 | 73.4 | 16.2 | 762 | 1 |
| I3_alone | [14.83, 14.83] | [-12.9, -12.9] | [0.87, 0.87] | 3.35 | 90.5 | 73.4 | 8.6 | 762 | 1 |
| I4_alone | [18.28, 18.28] | [-13.7, -13.7] | [0.84, 0.84] | 7.41 | 82.5 | 73.2 | 16.4 | 766 | 1 |
| I5_alone | [17.78, 17.78] | [-13.9, -13.8] | [0.99, 0.99] | 6.60 | 80.9 | 73.3 | 18.3 | 764 | 1 |
| I6_alone | [17.06, 17.06] | [-13.3, -13.3] | [0.86, 0.86] | 5.60 | 82.5 | 73.4 | 16.4 | 762 | 1 |
| cum_I1 | [17.45, 17.45] | [-13.6, -13.6] | [0.7, 0.71] | 6.13 | 82.7 | 73.7 | 16.3 | 752 | 1 |
| cum_I2 | [17.47, 17.47] | [-13.5, -13.5] | [0.67, 0.67] | 6.15 | 82.9 | 73.3 | 16.2 | 764 | 1 |
| cum_I3 | [17.47, 17.47] | [-13.5, -13.5] | [0.67, 0.67] | 6.15 | 82.9 | 73.3 | 16.2 | 764 | 1 |
| cum_I4 | [18.53, 18.53] | [-13.5, -13.5] | [0.67, 0.67] | 7.85 | 82.9 | 73.3 | 16.2 | 764 | 1 |
| cum_I5 | [19.38, 19.38] | [-14.6, -14.6] | [0.87, 0.87] | 9.55 | 79.1 | 73.5 | 20.2 | 759 | 1 |
| cum_I6 | [19.13, 19.13] | [-15.1, -15.1] | [0.8, 0.8] | 9.02 | 79.5 | 73.4 | 19.8 | 760 | 1 |

## Item notes
* **I1 analog static cells.** Static levels 9 points over the library range; gain +0.23 dB vs B0 but ~0 vs the library-limited B0L (17.46); AR 0.70 dB, XPD 28 dB (finer amplitude shaping); static fraction 73.7 %, diodes 752. Static cells may not be *added* beyond this: the diode penalty already puts every cell with a cheap static compromise to static.
* **I2 pre-phased switched cells.** Sub-library of ls 2.6/2.8 at position 1.6: asymmetric pairs, e.g. (m_OFF, m_ON) = (0.0, -1.2). Effective phase resolution (equivalent uniform step, rms error x sqrt(12)): 138 -> 127 deg: no useful gain; a reactive sheet radiates with the fixed position-dependent phase, so the per-cell pair offers amplitude, not phase, freedom. Gain +0.02 dB.
* **I3 radial taper.** Uniform illumination with edge power 15 % needs a mean alpha of about 7 Np/m; the cell supplies 0.7-1.0 (the solved C is 5-20 times too small in every sector, `taper_info_*` in the JSON), so the taper is saturated everywhere and a graded depth only removes conversion: -2.4 dB alone, -4 dB cumulative. **Not adopted; edge power stays 79-83 %.** The target ("edge <= 15 %") cannot be met with this cell; it needs a cell with alpha_beam of order 5-10 Np/m.
* **I4 launcher.** Repo PinLaunch model, probe top 2-10 mm x top-hat radius 0-3 mm (`launcher_opt.json`). Launch fraction >= 0.8 only for a 10 mm pin with a 2-3 mm top hat (0.80 / 0.82). Its probe impedance is far from 50 ohm: 167 - 286j (r 2 mm) / 83 - 220j ohm (r 3 mm); bare S11 -1.3 dB, with a lossless reactive match about -12 dB (R = 83 ohm), which still misses -15 dB. The bare 6 mm probe: 0.627, S11 -15 dB (HFSS calibrated). Best S11-compliant alternative 3 mm pin + 1 mm disc: 0.609, -14.2 dB: no gain. **I4 gain (+1.06 dB, 0.627 -> 0.80) is shown but flagged: S11 not met, needs an HFSS launcher check (pin height 10 mm = 0.35 lambda, anti-pad, match).**
* **I5 diode loss / position.** Position 0.8 mm (and 1.2) in the sub-library raises the swing and cuts the ON conductance; with the ON-weighted assignment penalty mu: ON weight 0.33 -> 0.20, SW->beam 16 -> 20 %, +0.85 dB cumulative, SLL -1.1 dB. This is the only item that moves conversion.
* **I6 edge rings.** The outer two rings (12.8 mm) as static, doubled target: -0.25 dB (gain), SLL -0.5 dB. They cannot radiate the 80 % edge power in 12.8 mm (alpha_max about 1.9 Np/m -> 5 % of the incident power). **Not adopted.** A resistive rim would terminate the power (+SLL) but not raise gain; not computed.
* **I7 joint.** eps_B, g and lambda re-optimised on the adopted chain (24 configs): the optimum is eps_B = 0, g = 2.7: 19.38 dBi, no gain over the chain (`coded_improve_i7.json`). Static fraction does not rise at fixed gain: 73.5 %, 759 diodes. Trading diodes for gain: 597 diodes (79.1 % static) 19.05 dBi but SLL -12.7 dB; 499 diodes (82.6 %) 18.77 dBi, SLL -10.7 dB (`coded_improve_fewer_diodes.json`).

## Final (adopted I1 + I2 + I4 + I5, I7 settings)
**19.38 / 19.38 dBi (18.32 without the flagged I4), SLL -14.6 dB, AR 0.87 dB, XPD 26.0 dB, aperture efficiency 9.55 % (8.0 % without I4), SW->beam 20.2 %, edge power 79.1 %, 73.5 % static, 759 diodes, 1 control line + inverter per axis.** Constraint (diodes <= 766, lines <= 1 per axis): OK.
The efficiency gap is not closed by I1-I7: 90 % of the surface-wave power is never radiated or absorbed (edge 79 %). The next lever is the cell: alpha_beam of 5-10 Np/m (larger ΔY, e.g. a loaded or multi-gap cell) and a lossy/radiating edge.

## HFSS checklist (before any claim)
Rim 0.2 mm cell, 4 states; position-0.8 mm / 1.2 mm / 2.0 mm cell variants (ON/OFF Y, loss); the 10 mm / 3 mm top-hat launcher (launch fraction, S11, match); the sub-library interpolation (reduced-swing cells); first-order leakage and scalar element factor on the board.
