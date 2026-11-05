# Fewer diodes: ideas D1-D6 side by side, K = 2 and K = 4 (CP set B) — MODEL_ONLY, nothing adopted

> **PENDING CELL VALIDATION (2026-11-05):** HFSS (air box, pec metal) shows the crossed two-diode rim-disc cell switches ONLY with both diodes ON (C/C 10.670, R/C 10.665, C/R 10.690, R/R 11.006 GHz at 90 deg); the model's per-axis independence is falsified. Every number in this document that depends on independent per-axis switching of this cell (alpha per axis, polarisation switching, board gain, diode counts) is on hold. See `results/sidequest_hybrid/CELL_INDEPENDENCE.md`.

Same RT5880 cell (rim disc, diodes at 0.8 mm), same engine and option sets (I1+I2+I5; K = 4 with 3 line classes per axis), `scripts/research/sidequest_diode_reduction*.py`, data `diode_reduction.json` / `_all.csv` / `_tables.md` (every row), figure `diode_reduction_pareto.png`, rotated cell `rotated_board_chi0.json`, `rotated_board_chi45.json`, tests `tests/test_sidequest_diode_reduction.py`.
**Baselines (reproduced):** K = 4: 1120 diodes, 3 + 3 lines (6 + 6 without inverters), 19.22 dBi mean (19.16-19.26), SLL -15.7 dB, XPD 22.4, ap. eff. 9.2 %. K = 2: 743 diodes, 1 + 1 lines (2 + 2), 20.24 dBi, SLL -16.1 dB, XPD 28.9, ap. eff. 11.6 %.
**Feasible** = SLL <= -12.7 dB and mean gain within 2 dB of the baseline and diodes <= baseline. Rows with SLL above -12.7 dB are shown faint in the figure and never reported as designs.

## What each idea does (alone, K = 2 / K = 4; best feasible row at the stated diode count)
| idea | implementation | K = 2 | K = 4 | verdict |
|---|---|---|---|---|
| base lambda sweep | diode penalty lambda | 536 diodes: 19.76 dBi, SLL -12.6 (fails); 743: 20.24 | 777: 18.23, SLL -13.5; 1120: 19.22 | the reference front |
| **D1** power-ranked sparsification | switched elements ranked by amp^2 x code disagreement, lowest forced static | 662: 19.64, SLL -13.2; below that SLL fails | 790: 17.35, SLL -12.3 (fails) | **worse than the lambda front**: the lambda penalty already ranks by marginal LS error, which weights by power |
| **D2** superpixels 1x2x / 1x2y / 2x2 | a group of cells shares one option (one diode pair) only where cheaper | **1x2y: 677 diodes, 20.37 dBi, SLL -15.9** (above the baseline); 528: 20.03, -16.0; 2x2: 535: 19.82, -13.4 | 1x2y: 709: 18.72, SLL -16.3; 2x2: 680: 18.34, -15.4 | **best single idea** |
| **D3** one-axis switching | no cell with both elements switched ('one'); strict x only | 'one': 660: 20.09, SLL -14.8; strict-x: SLL -3..-6 (fails everywhere, y static 100 %) | 'one': 945: 18.81, -14.4 | 'one' small gain; strict fails |
| **D5** state overlap | per-state +-2 deg pointing and free reference phases, coordinate descent | 569: 20.07, -13.2 (base at 536: 19.76, -12.6) | 762: 18.30, -13.3 (base 777: 18.23, -13.5) | +0.3-0.4 dB at equal diodes (K = 2), ~0 at K = 4; offsets found K = 2 (dphi -2/+2, dtheta -2/-2), K = 4 (dphi 0/-2/0/+2, dtheta 0/0/-2/0) |
| **D6** core + static rim | switched cells only inside radius R | R = 115 mm: 593: 19.67, -13.0; smaller R fails SLL | R = 115 mm: 966: 18.51, -15.0; R = 95: 747: 17.47, -12.4 (fails) | worse than the lambda front: the static rim cannot follow the code (SLL) |
| **D4** one-diode rotated disc | 8 axis angles, 1 diode/cell (fix fcc6691); alpha 0.50 Np/m (axis-aligned, exact) / 1.04 (45 deg, staircased) | chi0: 15.6 dBi at 675 diodes; chi45: **489 diodes: 18.12 dBi, SLL -14.2** (2.1 dB below baseline) | chi0: 13.0 at 706 (3 lines); chi45: 15.8 at 695, 16.3 at 962 | not competitive: the disc leaks 2-4x less than the crossed rim cell |
Lines: unchanged by every idea except D4 (K = 4: 3 per axis stays 3; D4 has one diode type: 3 total, 6 without inverters). D2 sampling (Nyquist = hologram phase advance per group pitch > pi): 1x2x: 28.7 % of the amplitude-weighted cells in the x-beam states (A, C), 0 in B, D; 1x2y the reverse; 2x2: ~28 % in all states. Grating lobes cannot appear (group pitch 12.8 mm = 0.45 lambda_0, lambda_0/p = 2.2 > 1); the aliasing shows up as raised SLL (K = 2 2x2 at 359 diodes: -9.9 dB) not as separate beams. D3: the fraction of cells that already keep a static y in the 'one' design is 75.7 % (K = 2, 660 diodes) and 77.0 % (K = 4); forcing 100 % costs the pattern.

## Best combinations (wave 1 + wave 2, D5 parameters carried)
Pareto front (SLL <= -12.7 dB), mean gain per state:
| K | diodes | gain per state [dBi] | SLL | XPD | ap. eff. | lines x / y (no inv.) | design |
|---|---|---|---|---|---|---|---|
| 2 | 743 (baseline) | 20.24 | -16.1 | 28.9 | 11.6 % | 1 / 1 (2 / 2) | |
| 2 | 625 | 20.46 | -15.9 | 27.6 | 12.2 | 1 / 1 | D2 1x2y + D5, lambda 0.32 |
| 2 | 509 | 20.16 | -15.8 | 28.8 | 11.4 | 1 / 1 | D2 1x2y + D5, lambda 0.43 |
| **2** | **484** | **20.09** | -14.2 | 28.6 | 11.2 | 1 / 1 | **D2 1x2y + D3 + D5, lambda 0.40 (best <= 500)** |
| 2 | 410 | 19.81 | -13.3 | 29.1 | 10.5 | 1 / 1 | D2 1x2y + D5, lambda 0.53 (smallest SLL-ok design) |
| 4 | 1120 (baseline) | 19.22 | -15.7 | 22.4 | 9.2 | 3 / 3 (6 / 6) | |
| 4 | 844-856 | 19.2 | -16.5 / -15.9 | 22 | 9.1 | 3 / 3 | D2 1x2x/1x2y + D5, lambda 0.67 (same gain, -24 % diodes) |
| 4 | 665 | 18.6 | -15.1 | 20.7 | 8.0 | 3 / 3 | D2 1x2y + D5, lambda 0.89 |
| 4 | 573-611 | 18.3-18.5 | -13.3..-13.8 | 21.9-22.5 | 7.4-7.7 | 3 / 3 | D2 1x2x + D3 + D5, lambda 0.83-0.89 |
| **4** | **515** | **18.0** (18.05 / 17.87 / 18.24 / 17.78) | **-12.8** | 19.8 | 6.9 | 3 / 3 (6 / 6) | **D2 1x2y + D5, lambda 1.11 (smallest SLL-ok design; 1.2 dB below baseline)** |
**Marked designs** (SLL <= -12.7 dB, loss <= 2 dB against the baseline, no more diodes than the baseline):
* **<= 360 diodes: none, at either K.** The smallest SLL-ok designs are 410 diodes (K = 2) and 515 (K = 4); below that the sidelobes exceed -12.7 dB (K = 2: 359 diodes -9.9 dB; K = 4: 447 diodes -12.6 dB; 299 diodes -7.1 dB).
* **<= 500 diodes: K = 2: 484 diodes, 20.09 dBi (-0.15 dB), SLL -14.2 dB; K = 4: none** (515 diodes, 18.0 dBi, SLL -12.8 dB is the nearest, 15 diodes over).
Read against the RIS estimate (150-360 diodes): at K = 2 the board needs 410-484 diodes for ~20 dBi, at K = 4 about 515-665 diodes for 18-18.6 dBi; the 300-500 diode / 18+ dBi window is reached at K = 2 (and at K = 4 only at its upper edge, 515 diodes, 18.0 dBi). The lever is code-sharing (D2 plus D5), not geometry (D6) or heuristic ranking (D1). 

## Caveats
Model only (first-order leakage, scalar element factor); the library shapes for the analog static levels and the second geometry are normal-incidence MoM rows; D2 shares one option over a group (the larger switched element is assumed to reach the same Y as the single cell: not solved); D5 offsets move the beam by up to 2 deg from nominal (pointing error reported in the CSV); D4 numbers depend on the rotated-disc alpha (0.50-1.04 Np/m, factor 2 between axis-aligned and staircased 45 deg at mesh 64); line counts are per axis; the wave-2 lambda grid is coarse near the SLL threshold.
