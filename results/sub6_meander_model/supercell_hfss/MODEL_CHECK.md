# Frequency-specific leakage model check (MODEL_ONLY)

The 40 mm PEC-lid unit-cell TM dispersion was inverted to sheet susceptance, then passed through `leakage_alpha_first_order_tensor` with the lid removed. The earlier anchor “expected range” of **1.158–1.480 Np/m was calculated at 5.8 GHz**, from two possible y-branch pairings. It was not calculated at the measured supercell frequency 5.682 GHz; the earlier report's claim that the 5.682 GHz value was within that expected range was unsupported.

| Anchor HFSS f (GHz) | Model alpha at 30 deg (Np/m) | HFSS alpha (Np/m) | Model/HFSS | Unit-cell x/y bracketed? |
|---:|---:|---:|---:|---|
| 5.594682 | 2.01884 | 1.1119 | 1.816 | yes/no |
| 5.681668 | 1.67573 | 1.4530 | 1.153 | yes/no |
| 5.752876 | 1.28460 | 1.8311 | 0.702 | yes/yes |

The model trends downward with frequency while HFSS trends upward. It is a screening approximation, not a calibrated frequency predictor. The y branch is extrapolated at 5.595 and 5.682 GHz from 85/95-degree unit-cell endpoints. The supercell extraction uses the matched negative-fold per-cell phases 88.125/86.25/84.375 deg and Q_base 1939.74/2013.42/2087.10. The 45-degree baseline phase is extrapolated 0.625 deg below the 85-degree endpoint. The group velocity is 176.165 Mm/s from the 15-to-45-degree supercell slope using **92.8 mm**, not 11.6 mm. These checks do not explain the opposite frequency trends. The supercell branch identity and the frequency-dependent sheet approximation remain uncertainties.

| 5.8-GHz beta-pair branch | beta_x/k0 | beta_y/k0 | alpha at 25/30/35/40 deg (Np/m) |
|---|---:|---:|---|
| Anchor, y main | 1.13153 | 1.07006 | 1.043 / 1.077 / 1.103 / 1.109 |
| Anchor, y alternate | 1.13153 | 1.06501 | 1.376 / 1.416 / 1.442 / 1.442 |
| B2p2 D3p0 | 1.07889 | 1.11278 | 0.336 / 0.346 / 0.353 / 0.353 |
| B2p2 D3p6 x low | 1.11204 | 1.16519 | 0.262 / 0.275 / 0.292 / 0.311 |
| B2p2 D3p6 x high | 1.05700 | 1.16519 | 3.032 / 3.152 / 3.255 / 3.308 |

These are MODEL_ONLY values at 5.8 GHz. Mode coupling and which branch a pin launches are not established by the beta-pair chain. The two D3p6 x branches give substantially different leakage predictions, so both require HFSS tracking.

Full inputs, outputs, and angle points: `model_check.csv` and `model_check.json`. Reproduction script: `scripts/research/sub6_supercell_model_check.py`.
