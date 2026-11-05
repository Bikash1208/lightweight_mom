# SLL fixes from the board itself (no rim, no extra parts) - HFSS far-field analysis, cone test, ranking, package r 2.30 - MODEL_ONLY (2026-11-06)
Package built (model reaches SLL <= -13 dB in both states): `hfss_exports/sidequest_board_RO3006_r2p30/` -> `C:\hmsa_sync\data\sidequest\board_RO3006_r2p30\`. **Read section 5 first: the model underestimates the HFSS cone by 10-17 dB, so the model SLL is not a prediction of the HFSS SLL.** Literature: `docs/briefs/SIDEQUEST_SLL_LIT.md`.

## 1. HFSS far fields of board_K2_default (`sidequest_ffd_analysis.py`, `hfss_ffd_analysis.json`, `hfss_ring_levels.json`; ffd 91 theta x 181 phi, 2 deg)
* State 1: D_peak 9.47 dBi at (24, 2); state 2: 8.27 dBi at (24, 176). Power in theta < 30 / 30-60 / 60-90 / back: 13.5 / 15.9 / 36.6 / 34.1 % (S1), 10.2 / 13.8 / 42.3 / 33.8 % (S2); ring 60-90 and back are uniform in azimuth (2.6-3.9 % per 30 deg sector), D about +1.5 dBi and -1.7 dBi on average: this is the axisymmetric diffuse floor (direct probe + edge, 71-76 %), NOT the strongest sidelobes.
* **97 % of the radiated power is common to both states** ((E1+E2)/2 against (E1-E2)/2). Strongest sidelobes (S1): (20,24) -0.7 dB, (26,-34) -1.1, (12,140) -1.9, (26,-146) -1.9, (26,94) -2.4 ...; (S2): (20,154) -0.1, (26,-148) -1.2, (12,40) -1.2 ...
* **Azimuthal MEDIAN level relative to the peak: theta 12: -4.6 / -4.2 dB; theta 24: -5.7 / -6.3 dB; theta 16-20: -8..-9 dB; theta 30: -10 dB (S1 / S2)** = a conical ring at theta 24 and a stronger one at theta 12, as you said. The earlier "state-independent floor 8.3 dBi" assumption is withdrawn (the floor is about 1.6 dBi); the non-beam floor model now uses a peaking factor 2.0 instead of 7.45.

## 2. Cone test in the model (`sidequest_cone_test.py`, `k2_cone_test.json`; baseline RT5880 board)
* The static cells ARE the model's cone: azimuthal median relative to the peak at theta 24: full -16.1 dB, **static-only -15.6 dB**, switched-only -21.5 dB; theta 12: full -21.6, static-only -22.8, switched-only -29.8 (S1; S2 the same within 1 dB). The static field peaks at (25, 180) in both states (-5.0 / -5.2 dB): it is the state-common compromise of the two targets.
* **But the level is 10-17 dB too low against HFSS** (HFSS median -5.7 / -6.3 dB at theta 24 and -4.6 / -4.2 dB at theta 12 against the model -16 / -22). Scaling the statics by x1.33 / x2 (decoupled depth) raises the model ring by only 1.8 / 3 dB. The size of the cone (and the theta 12 ring) is therefore NOT explained by the modelled static compromise alone; the mechanism that makes it 10-17 dB stronger in HFSS (cell mutual coupling and spatial dispersion not in the library, higher-order harmonics of the 7-level hologram, radial standing waves, feed) is unidentified.

## 3. Fixes, ranked on the model (RO3006, decoupled depth: statics at full MoM depth, switched pairs at the calibrated contrast; <= 514 diodes, one line; r 2.30 baseline F0 = decoupled, design_null, g 2, mu 70: SLL total -8.4 / -8.9 dB, gain 14.5 / 13.7 dBi, iso 12.7 / 12.6)
| rank | method | source DOI | dSLL S1 / S2 [dB] | dgain S1 / S2 [dB] | fab impact | modelled |
|---|---|---|---|---|---|---|
| 1 | **Neutral static cells where the two states disagree** (statics at the level of least radiation; penalty kappa x disagreement x m^2 in the joint optimisation) | none (own proposal) | -3.4 / -4.4 (kappa 2) | +0.9 / -1.0 | none (cell choice) | y |
| 2 | Inverse taper (target amplitude rising towards the rim, 1 + 0.5 (rho/rho_edge)^2) | 10.1038/s41598-020-74674-w | -0.7 / -0.3 | 0 / +0.1 | none | y |
| 3 | Edge-ring boost x1.5 (outer 12.8 mm radiates the leftover wave into the beam) | 10.1038/s41598-020-74674-w, 10.1109/TAP.2021.3090801 | -0.2 / -1.4 | +0.2 / +0.2 | none | y |
| 4 | Ring-penalty joint optimisation (theta 8-32 cone outside 12 deg of the own beam, quartic, coordinate descent over the cell options) | 10.1109/TAP.2015.2399937, 10.1038/s41467-017-00164-9 | 0 / +0.8 | 0 / -0.1 | none | y |
| 5 | More diodes (676 / 892 instead of 514) | - | -0.5 / -0.8 (at ring CD) | +0.6 / +0.6 | +160 / +380 diodes | y |
| 6 | Larger leakage: r 2.30 instead of 2.22 (A0 11.3 / 3.3 uncorrected-decoupled, edge 8 % / 60 %) | 10.1038/s41598-020-74674-w | -0.7 / +0.9 (design_null, g 2) | +3.4 / +1.2 | cell radius | y |
| 7 | Decoupled static depth (against the coupled calibration of all levels) | user proposal | +2.5 worse | +0.8 | none | y |
| 8 | Launcher with less direct radiation (pad, array launcher, cavity feed) | 10.1109/TAP.2021.3090801, 10.1038/ncomms10360 | about 0 in the model (floor not binding, -10..-16 dB) | +1.3 (launch 0.96) | pad / feed structure | y (floor only); the cone is not a feed effect in the model |
| 9 | Via fence / shorted wall edge, corrugated or choke edges, extended ground | no journal source verified | NA | NA | ground plane | n |
Combination (best in the scan, 514 diodes): **r 2.30, kappa 2, g 2.0, mu 30, inverse taper 0.5, ring x1.5: SLL total -13.6 / -13.2 dB** (aperture -13.6 / -13.2, floor -16.3 / -13.8), gain 15.7 / 13.2 dBi, isolation 22.7 / 17.4 dB, SW to the edge 21 %, SW to the beam 73 %, A0 16.9 Np/m, beta/k0 1.404. Neighbouring points of the same scan give -12.0 .. -13.4 dB in the worse state (one point of about 30; the best-of-many selection is worth about 1 dB). At r 2.26 the best is -11.0 / -14.1 dB (A0 8.5, beta/k0 1.28). At r 2.22 (A0 3.3, beta/k0 1.21) the best with neutral statics is -8.8 / -12.2 dB.
The mechanism of fix 1: kappa 2 sets 1033 of the 1093 static cells to the neutral level (open slit 0.929, m about 0), so the state-common field is mostly removed and the beams are made by the 257 switched cells (514 diodes) plus the few non-neutral statics: a SPARSE hologram (19 % of the cells radiate the beam).

## 4. Diode count against the SLL (model, ring-penalty CD, decoupled, g 2 mu 70; the search stops where the cost stops improving)
| diodes | r 2.22 SLL S1 / S2 | r 2.26 | r 2.30 |
|---|---|---|---|
| 514 | -9.3 / -9.2 | -9.6 / -8.8 | -8.2 / -8.1 |
| 700 | -9.1 / -10.8 | -10.4 / -10.2 | -8.9 / -9.5 |
| 900 | -9.5 / -12.2 | -10.5 / -11.5 | -9.2 / -9.7 |
| all cells allowed (design settles at 892-938) | -9.5 / -12.2 | -10.2 / -11.4 | -9.2 / -9.7 |
More diodes alone do not reach -13 dB; the fix that does is the neutral static choice (514 diodes). No design in this study needs more than 514 diodes.

## 5. Limits and what it would take
* **The model SLL cannot be trusted at the -13 dB level**: for the HFSS-measured board (RT5880) the model gives -11.7 / -11.8 dB against the HFSS strongest sidelobes -0.7 / -0.1 dB; the HFSS ring (theta 12 / 24 median -4.6 / -5.7 dB) is 10-17 dB above the model ring. A model result of -13.6 / -13.2 dB therefore means "the modelled hologram has no modelled sidelobe above -13 dB", not an HFSS SLL.
* Sparse neutral-static design: its benefit rests on the statics being really neutral (open slit 0.929 at m = +0.05 in the MoM normal-incidence library; its neutrality at the HFSS level is untested and the OFF-branch error found in the cell check applies to open-slit cells: a 10-15 % level error of the many neutral cells is a cone of its own); the state imbalance is 2.5 dB (15.7 / 13.2 dBi); beta/k0 1.404 (bound-mode and spatial-dispersion risk, K2_POINTING.md; r 2.30 is on the A0 plateau, pointing breaks up at beta/k0 >= 1.40 under tolerance); r 2.30 has not been HFSS-checked as a cell.
* What it would take to trust a number: three HFSS runs on the existing machinery: (a) board with the probe + pad and NO cells (direct radiation and any theta 12 ring from the feed / edge); (b) board with only the static cells (switched cells replaced by neutral statics): it shows whether the cone is the state-common field; (c) the r 2.30 cell eigenmode check (OFF / ON at 90 deg) to tie the library. With (b) the model's cone amplitude can be calibrated (factor 10-17 dB) and the fixes re-ranked against it.
