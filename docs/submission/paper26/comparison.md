# 26 GHz 1-bit coded 5-port board: this model vs Mishra et al. (TVT 2026)

Cell pair read from the paper's Fig. 2 (cell 1 = bare substrate, cell 2 = square patch with annular slot), not from its text.

## Cell pair, 2.57 mm period

| f (GHz) | empty phase (deg) | slotted phase (deg) | separation (deg) | empty |S11| (dB) | slotted |S11| (dB) | |S11| difference (dB) | beta empty (rad/m) | beta slotted (rad/m) |
|---|---|---|---|---|---|---|---|---|
| 25 | 56.38 | -130.16 | 173.45 | -0.009 | -0.009 | 0.000 | 580.16 | 666.53 |
| 26 | 48.96 | -135.16 | 175.89 | -0.010 | -0.008 | 0.002 | 607.92 | 702.29 |
| 27 | 41.30 | -139.48 | 179.22 | -0.012 | -0.007 | 0.004 | 636.11 | 738.13 |

## Surface-wave beta modes at 26 GHz

- Paper mode, beta = k0 sqrt((eps_r+1)/2): 689.27 rad/m.
- Model mode, beta of the mixed 1-bit surface (equal-area mean of the two sheet admittances): 662.07 rad/m.
- Aperture-edge phase difference between the modes at the 28.8 mm edge of a 5 lambda_0 sub-aperture: 44.9 deg.

## Board

| | Paper | This model, paper mode | This model, model mode |
|---|---|---|---|
| Board | 20 lambda_0 x 10 lambda_0 | 230.6 x 115.3 mm | 230.6 x 115.3 mm |
| Cells | 92 x 46 | 90 x 45 | 90 x 45 |

## Beams: pointing and LHCP/RHCP (per requested beam)

The paper gives the target directions, not measured beam angles in its text (measured patterns are in its Fig. 5 only). Levels are relative to each port's own peak; the paper's gains are realised dBic, ours are lossless directivity (dBi).

| port | requested (theta, phi) | realised (theta, phi) | pointing error (deg) | LHCP level (dB rel. port peak) | RHCP suppression (dB) | directivity (dBi, lossless) | paper realised LHCP gain (dBic) | mode |
|---|---|---|---|---|---|---|---|---|
| 1 | (15, 0) | (14.0, 0.0) | 1.0 | 0.00 | 26.7 | 20.28 | 16.8 | paper |
| 2 | (30, 90) | (29.0, 88.8) | 1.2 | 0.00 | 16.3 | 19.53 | 16.5 | paper |
| 3 | (45, 180) | (43.0, 182.5) | 2.6 | 0.00 | 7.9 | 18.36 | 15.5 | paper |
| 4 | (60, 270) | (55.0, 271.2) | 5.1 | 0.00 | 3.6 | 18.18 | 14.8 | paper |
| 5 | (15, 0) | (16.0, 335.3) | 6.6 | -0.31 | 3.3 | 13.63 | 15.8 | paper |
| 5 | (30, 90) | (28.0, 103.6) | 6.9 | 0.00 | 1.2 | 13.94 | 15.8 | paper |
| 5 | (45, 180) | (47.0, 172.6) | 5.7 | -0.16 | 1.5 | 13.79 | 15.8 | paper |
| 5 | (60, 270) | (54.0, 281.1) | 11.1 | -0.67 | -0.3 | 13.28 | 15.8 | paper |
| 1 | (15, 0) | (15.0, 0.0) | 0.0 | 0.00 | 23.8 | 21.23 | 16.8 | model |
| 2 | (30, 90) | (29.0, 88.8) | 1.2 | 0.00 | 11.2 | 20.60 | 16.5 | model |
| 3 | (45, 180) | (44.0, 177.5) | 2.0 | 0.00 | 7.9 | 19.01 | 15.5 | model |
| 4 | (60, 270) | (55.0, 271.2) | 5.1 | 0.00 | 4.8 | 18.17 | 14.8 | model |
| 5 | (15, 0) | (15.0, 0.0) | 0.0 | -0.38 | 40.0 | 15.96 | 15.8 | model |
| 5 | (30, 90) | (30.0, 88.8) | 0.6 | -0.00 | 12.5 | 16.34 | 15.8 | model |
| 5 | (45, 180) | (43.0, 182.5) | 2.6 | 0.00 | 5.1 | 16.34 | 15.8 | model |
| 5 | (60, 270) | (57.0, 271.2) | 3.2 | -0.94 | 3.2 | 15.40 | 15.8 | model |

RHCP suppression is shown as the RHCP level below LHCP at the beam, so a positive number means RHCP is that many dB down. Paper: 10-15 dB.

## Other paper targets

| Quantity | Paper | This model |
|---|---|---|
| Isolation between ports | < -30 dB (measured) | not computed: first-order edge-gap flags only; check in HFSS with the STEP |
| |S11| | < -10 over 25-27 GHz | not comparable: a reflectarray cell reflects almost all power (|S11| about 0 dB by construction); the paper's |S11| is a feed-match quantity |

## Port 5 equalisation

Beam levels (dB relative to port 5's peak) per iteration:

- paper mode: [0.0, -0.51, -1.76, -4.11]; [0.0, -1.05, -0.91, -0.05]; [-0.31, 0.0, -0.16, -0.67]
- model mode: [0.0, -1.64, -4.12, -6.05]; [-1.66, -1.82, -1.03, 0.0]; [-0.38, -0.0, 0.0, -0.94]

