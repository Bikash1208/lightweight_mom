# K = 4 coded board on RT5880 (model only, round 13) — MODEL_ONLY
Cell: rim disc, diodes at 0.8 mm (A0 = 1.84 Np/m, AA0 0.094, background loss 0.15, beta 1.079 k0), launch 0.627, hexagon side 150 mm, centre probe, theta 25 deg, phi 0/90/180/270 deg; the same engine as round 11/12 generalised to K states (`sidequest_coded_k4.py`, `..._run.py`; data `coded_k4_{K2,K4a,K4b}.json`, `coded_k1_dedicated.json`; tests `tests/test_sidequest_coded_k4.py`). Options: static levels (I1: 9 analog levels), switched cells = K-bit ON-patterns ("words", both complements) of a diode cell, I2 = two geometries (ls 2.8 / 2.6 at position 0.8 mm), I5 = ON-weighted assignment penalty. The library shapes are the round-11 normal-incidence MoM rows renormalised to the chosen cell (flagged). Lines per axis = distinct ON-patterns of the switched elements of that axis (without inverters) / complement classes (with inverters). No diode cap for K = 4 (K = 2 reference: 743 diodes = 371.5 per state).
**Normalisation note:** with the chosen pos-0.8 cell as the unit cell, I1/I2/I5 no longer add gain at K = 2 (20.56 -> 20.24 dBi, SLL -14.2 -> -16.1 dB, AR 0.83 -> 0.62); they reduce the diodes at K = 4 (below). Earlier K = 2 numbers (21.5 dBi) mixed a pos-1.6 library with the pos-0.8 alpha; the consistent K = 2 value is 20.2-20.6 dBi.

## Results (diodes tuned by lambda; I1I2I5 at lambda 0.555 as K = 2 bit design)
Polarisation set **B = RHCP/LHCP/RHCP/LHCP** (phi 0/90/180/270 deg), epsilon (0, 0, 0, 180) deg:
| design | gain [dBi] | SLL [dB] | XPD [dB] | ap. eff. [%] | edge [%] | static [%] | diodes (per state) | lines x / y (with inv.; without) |
|---|---|---|---|---|---|---|---|---|
| bit (lambda 0.555) | 20.8 / 20.7 / 20.8 / 20.7 | -14.1..-15.9 | 25-26 | 13.1 | 58 | 36.1 | 1827 (457) | 7 / 7 (14 / 14) |
| **I1+I2+I5** | **19.8 / 19.8 / 19.7 / 19.9** | **-15.2..-15.9** | 23-24 | **10.5** | 63 | **55.5** | **1274 (318)** | 7 / 7 (14 / 14) |
| K = 2 reference (I1I2I5) | 20.2 / 20.2 | -16.1 | 28.9 | 11.6 | 63 | 74.0 | 743 (371.5) | 1 / 1 (2 / 2) |
Static against SLL (K = 4b, I1I2I5 option set, lambda scan): 
| lambda | static % | diodes | gain [dBi] | SLL [dB] |
|---|---|---|---|---|
| 0.2 | 29.3 | 2021 | 20.85 | -14.1 |
| 0.4 | 44.4 | 1589 | 20.4 | -15.2 |
| 0.555 | 55.5 | 1274 | 19.8 | -15.2 |
| **0.8** | **69.4** | **876** | **18.7** | **-14.5** |
| 1.2 | 82.1 | 513 | 17.1 | -9.7 (fails -12.7) |
| 2.0 | 93.0 | 200 | 14.9 | -5.3 |
At -12.7 dB SLL the static fraction reaches about 69-75 % for K = 4 (gain 18.7 dBi, 876 diodes = 219 per state, below the K = 2 per-state count) at a cost of 1.1-2.1 dB relative to 1274-1589 diodes. Lines stay 7 per axis (14 without inverters), the 2^K - 2 limit halved; they do not shrink with static fraction.

Polarisation set **A = RHCP/LHCP/LP-H/LP-V**: states A and B (CP) are fine (19.5 dBi, SLL -11.3..-12.5), **states C and D (LP) fail**: gain 16.6 / 17.5 dBi, **SLL -3.8 / -5.2 dB** (bit design: 17.4 / 18.3 dBi, SLL -4.2 / -5.0). The crossed cell has diagonal Delta Y only, so the radiated current J = Delta Y rho_hat cannot follow a linear polarisation over the whole aperture (a dedicated static single-state LP board reaches XPD only 5 dB); LP states need the xy-coupled (rotated) element, i.e. the rotated-disc cell of rounds 8-9. Static 67.7 %, 924 diodes, lines 7 / 4 (14 / 4). With a diode in every element (full 1-bit) the LP states are usable (17.7 / 18.4 dBi, SLL -14.5 / -14.2 dB, XPD 13-17 dB) at 2860 diodes.

## Comparison (same cell, model; CP set B)
| board | gain per state [dBi] | diodes | feeds / RF switch | lines |
|---|---|---|---|---|
| coded, 55.5 % static | 19.8 | 1274 | 1 / none | 7 + 7 (with inverters) |
| coded, 69.4 % static (SLL -14.5) | 18.7 | 876 | 1 / none | 7 + 7 |
| full 1-bit (diode in every element) | 20.9 | 2860 | 1 / none | 8 + 8 (16 + 16 without inverters) |
| dedicated single-state board (ceiling, no switching) | 21.0 | 0 | 1 | 0 |
| 4 feeds + SP4T, quadrant per state (estimate) | **14.2** = 21.0 - 6.0 (quarter area) - 0.8 (SP4T loss, assumed) | SP4T only (~4 junctions, assumed) | 4 / SP4T | 2-3 logic |
The segmented number is an estimate (quarter aperture with the same conversion fraction; feed spacing, gate isolation and the 0.8 dB SP4T loss are assumptions). The coded board is about 5.6 dB above it and about 1.2 dB below the full 1-bit board with 55 % fewer diodes.

## Risks
Model-level only; library shapes from the normal-incidence rows; LP states unsupported by this cell; K = 4 line count 7 per axis is the design rule, not an optimum (not minimised by a word-merging step here).
