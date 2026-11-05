# K2-default board after HFSS state 1: leakage / edge / launcher options - MODEL_ONLY (2026-11-05)
HFSS state 1 (150 mm RT5880, 18 min adaptive, 1.28 M tets, 53 GB): beam (24, 2) deg, 9.1 dBi (3.5 % aperture efficiency); 13 % of the power at theta < 30, 32 % at 60-90, 35 % back; strongest sidelobe -0.8 dB. Model (design (i)): 11.19 dBi, SW->beam 13.1 %, 82.5 % of the launched SW reaches the edge, launch 0.627 (so 37 % direct, 52 % edge, 8 % beam of the input).
Budget (input = 1): direct 0.373 + edge 0.627 x 0.825 = 0.89 non-beam, against HFSS 0.87 outside theta < 30. The model's SW->beam fraction is right; the leakage is too weak.
## Calibration to the ONE HFSS result (flagged, single point)
gain offset -2.09 dB; ring (60-90) 36 % and back 39 % of the non-beam power; the non-beam pattern has peaking factor 7.6 (set by the -0.8 dB sidelobe). Every predicted gain below = model - 2.09 dB; SLL = max(aperture SLL of the model, calibrated non-beam level). Edge/back/SLL are therefore indicative, the gain trend is the model's.
Method: `scripts/research/sidequest_leakage_study.py` (design (i), g 2, mu 30, eps_B 180, lam re-tuned to <= 514 diodes unless stated), data `k2_leakage_study.json`, `k2_leakage_table.json`.
## Results (gain = calibrated min-state; ring/back/beam30 = % of input power; SW at edge = engine, % of launched SW)
| variant | diodes | gain [dBi] | SW at edge | ring 60-90 | back | theta<30 | SLL [dB] |
|---|---|---|---|---|---|---|---|
| baseline (HFSS: 9.1 / 32 / 35 / 13 / -0.8) | 514 | 9.1 | 82 % | 32 | 35 | 13 | -0.8 |
| (1) A0 x1.8 (1.5 Np/m) | 508 | 11.3 | 74 % | 30 | 33 | 21 | -3.3 |
| A0 x3 (2.5) | 510 | 13.2 | 63 % | 28 | 30 | 32 | -5.6 |
| A0 x4.7 (4.0) | 506 | 14.6 | 49 % | 25 | 27 | 45 | -7.5 |
| A0 x7 (6.0) | 508 | 15.7 | 37 % | 22 | 24 | 57 | -9.1 |
| A0 x10.7 (9.0) | 514 | 16.4 | 23 % | 19 | 20 | 70 | -10.5 |
| **A0 x14 (12.0): edge <= 15 %** | 508 | 16.6 | 15 % | 17 | 18 | 76 | -10.8 |
| A0 x18 (15.0) | 514 | 16.4 | 9 % | 16 | 17 | 80 | -11.2 |
| (2) absorber rim 12.8 mm, 10 dB | 470 | 8.5 | 15 % | 17 | 18 | 13 | -3.0 |
| absorber rim 20 dB | 442 | 8.1 | 4 % | 14 | 16 | 13 | -3.3 |
| (3) launch 0.8 | 514 | 10.2 | 82 % | 31 | 34 | 17 | -2.0 |
| launch 0.9 | 514 | 10.7 | 82 % | 30 | 33 | 19 | -2.6 |
| combined: A0 x4.7 + rim 10 dB + launch 0.8 | 506 | 15.4 | 8 % | 10 | 11 | 58 | -11.5 |
## Findings
* Leakage is the lever: A0 is the cell's leakage rate for the 50/50 pattern (0.844 Np/m now). The 15 % edge target on the 150 mm board needs A0 about 12 Np/m (x14, modulation depth about x3.8 if A0 goes with the square of the X swing); the gain saturates at about 16.5 dBi (SW->beam 77 %); the diode count stays at about 510 (lam re-tuned).
* An absorber rim alone does not help the gain (it removes the edge SW and some outer cells: 8-9 dBi) and cannot go below the direct-probe floor: the 37 % direct share stays. It only helps the ring/back share and the sidelobe, and only together with stronger leakage.
* Launcher: 0.627 -> 0.8 / 0.9 gives +1.1 / +1.6 dB; the direct share (1 - launch) is the floor of back radiation and sidelobes (about 35 % of input at the baseline). A better launcher (choked / cavity-backed probe, or the side-fed Yagi of Zhao et al.) is worth about +1.5 dB and removes about 1/3 of the back power; it does not fix the leakage.
* Stack: the board cell (ls 2.8 OFF/ON) on RT5880 has X (lid) 121 -> 100 ohm, a swing of only 21 ohm (`k2_stack_contrast.json`). The same disc cell on RO3006 / RO3010 is NOT usable: the disc resonates there (Y of the wrong sign, ON Y 70-99 mS, no X root), so a stack route needs a re-scaled cell. Indicative only: the 10.5 GHz leaf cell on RO3006 with the HFSS-confirmed MA4AGP907 model had X 109 -> 17 ohm (swing 92 ohm, mesh-aligned but the leaf is the reconstructed shape), i.e. x4.4 in swing, x19 in A0 if A0 ~ (dX)^2: the order that is needed. RO3010 0.635 / 1.27: NA (the leaf run crashed on the Newton step; stopped per instruction). Nothing here is validated; the RT5880 cell stays the HFSS-confirmed one.
* Honest limit: the model has no edge-diffraction or probe-pattern physics. The ring / back split and the non-beam peaking come from one HFSS result; gain trend and SW-edge numbers are the model's.
