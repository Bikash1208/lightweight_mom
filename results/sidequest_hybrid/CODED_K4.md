# K = 4 coded board on RT5880 (model only, round 13) — MODEL_ONLY

> **PENDING CELL VALIDATION (2026-11-05):** HFSS (air box, pec metal) shows the crossed two-diode rim-disc cell switches ONLY with both diodes ON (C/C 10.670, R/C 10.665, C/R 10.690, R/R 11.006 GHz at 90 deg); the model's per-axis independence is falsified. Every number in this document that depends on independent per-axis switching of this cell (alpha per axis, polarisation switching, board gain, diode counts) is on hold. See `results/sidequest_hybrid/CELL_INDEPENDENCE.md`.

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

# Round 14 — minimum line count and a real quadrant (4 probes + SP4T) model

## A. Minimum number of control lines, K = 4 (CP set B)
A line drives one complement class of ON-patterns (with an inverter); the 14 non-constant 4-bit ON-patterns form **7 classes**. Design restricted to L allowed classes on both axes (`sidequest_k4_lines.py`; I1+I2+I5 options, lambda 0.555, eps (0, 0, 0, 180) deg; exhaustive over class subsets for L <= 3, greedy backward elimination for L = 4-6; score = mean gain - 0.5 dB per dB of SLL above -12.7 dB):
| L (per axis) | lines without inverters | mean gain [dBi] | worst SLL [dB] | static [%] | diodes | best classes |
|---|---|---|---|---|---|---|
| 7 | 14 | 19.78 | -15.2 | 55.5 | 1274 | all |
| 6 | 12 | 19.66 | -14.9 | 56.4 | 1248 | |
| 5 | 10 | 19.54 | -15.2 | 57.8 | 1206 | |
| 4 | 8 | 19.38 | -16.0 | 59.3 | 1164 | 0011, 0101, 0110, 0111 |
| **3** | **6** | **19.22** | **-15.7** | **60.8** | **1120** | **0011, 0101, 0110** (the three two-of-four patterns) |
| 2 | 4 | 17.74 | **-7.1 (fails)** | 72.2 | 796 | 0011, 0101 |
| 1 | 2 | 15.49 | -0.1 (fails) | 85.4 | 417 | 0110 |
**Minimum: 3 lines per axis with inverters (6 without), 0.56 dB below the 7-line design**, SLL -15.7 dB, XPD 22-25 dB, 1120 diodes (280 per state). L = 2 fails the -12.7 dB SLL limit for every class pair and every epsilon in {0, pi}^3 (168 designs, best -7.3 dB, `coded_k4_lines_eps.json`). If x and y diodes of the same class share a physical line the total is 3 lines + inverters (6 if the axes are wired separately); the K = 2 board needs 1 per axis. Subsets for L = 4-6 are greedy (not proven optimal); the search is over the same class set on both axes and a fixed lambda (static fraction adapts).

## B. Segmented comparator: 4 probes, quarter-aperture holograms
Model (`sidequest_quadrant_board.py`, `quadrant_board.json`): state k uses the 90-degree wedge of the hexagon centred on phi_k (about 380 cells, 14.3-15.0e3 mm^2), its own probe on the wedge axis at 0.75 of the centroid distance (62 mm from the board centre; 0.3-1.3 scanned, best 0.75-1.0), omnidirectional launch (launch fraction 0.627), only the wedge's cells radiate (**ideal isolation of the other wedges: an upper bound**), static cells only (0 diodes; 1-bit and 9-level options, best chosen), the same cell, alpha, loss and power bookkeeping as the coded board with the radial geometry measured from the probe. Gain = geometric gain - IL(SP4T) (parameter).
| state | pol | geometric gain [dBi] | SLL [dB] | XPD [dB] | AR [dB] | ap. eff. (wedge) [%] | SW->beam [%] | edge [%] |
|---|---|---|---|---|---|---|---|---|
| A (phi 0) | RHCP | 12.86 | -14.4 | 21.8 | 1.42 | 8.4 | 20 | 76 |
| B (phi 90) | LHCP | 12.48 | -15.4 | 25.3 | 0.94 | 8.0 | 20 | 77 |
| C (phi 180) | RHCP | 12.86 | -14.4 | 21.8 | 1.42 | 8.4 | 20 | 76 |
| D (phi 270) | LHCP | 12.48 | -15.4 | 25.3 | 0.94 | 8.0 | 20 | 77 |
Gain against SP4T insertion loss:
| IL [dB] | segmented A/B [dBi] | coded 7 lines [dBi] | coded 3 lines [dBi] | coded - segmented (7 lines) |
|---|---|---|---|---|
| 0 | 12.9 / 12.5 | 19.8 | 19.2 | +6.9 / +7.3 |
| 0.8 | 12.1 / 11.7 | 19.8 | 19.2 | +7.7 / +8.1 |
| 1.5 | 11.4 / 11.0 | 19.8 | 19.2 | +8.4 / +8.8 |
| 2.0 | 10.9 / 10.5 | 19.8 | 19.2 | +8.9 / +9.3 |
The coded board is **6.9-9.3 dB** above the segmented board; the earlier estimate (14.2 dBi = dedicated 21.0 - 6.0 - 0.8) was optimistic by about 2 dB, because a quarter-aperture board converts only 20 % of the probe power (edge 76 %) against 36 % for the full aperture, with the same cell alpha. The segmented board needs 0 diodes in the aperture (plus the SP4T and four probes). Caveats: omnidirectional probe at one position per wedge (a directional launcher or gate ring, not modelled, would raise it), ideal isolation, no coupling between probes.
