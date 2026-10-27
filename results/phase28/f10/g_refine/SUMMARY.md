# G screening and package-inductance sensitivity (10.0 GHz, mesh 64, MODEL_ONLY)

States 0 and 1 stand for all four. 'dX design' = anisotropy at the item-2 design loading (beta/k0 1.06); 'alpha' = leakage at the cell's OWN loading. Mesh 64 only (G's 64 -> 96 change is about 1 %); no variant was confirmed at mesh 96. Variant B3 (hub 1.0 + neck gap 0.6) is geometrically invalid (zero-width diagonal neck), A-variant gap_t 0.4 likewise: not evaluated.

| name | state | b_bar [mS] | delta [mS] | dX own [ohm] | dX design [ohm] | alpha [Np/m] |
|---|---|---|---|---|---|---|
| A1_neck_width_0p6 | 0 | 4.47 | 2.308 | 44.3 | 89.3 | 0.404 |
| A1_neck_width_0p6 | 1 | 4.46 | 2.294 | 44.0 | 88.7 | 0.397 |
| A2_diag_outer_3p0 | 0 | 4.96 | 2.493 | 52.7 | 98.2 | 0.629 |
| A2_diag_outer_3p0 | 1 | 6.37 | 3.895 | 121.4 | 186.2 | 3.419 |
| A3_hub_half_1p0 | 0 | 4.79 | 2.701 | 56.0 | 108.7 | 0.668 |
| A3_hub_half_1p0 | 1 | 4.90 | 2.807 | 59.7 | 114.3 | 0.768 |
| A4_hub_half_1p4 | 0 | 4.53 | 2.245 | 43.5 | 86.4 | 0.397 |
| A4_hub_half_1p4 | 1 | 4.23 | 1.942 | 35.3 | 73.0 | 0.250 |
| A5_neck_gap_0p2 | 0 | 4.58 | 2.272 | 44.4 | 87.7 | 0.417 |
| A5_neck_gap_0p2 | 1 | 4.29 | 1.985 | 36.5 | 74.8 | 0.270 |
| B1_hub1p0_diag3p0 | 0 | 5.14 | 2.737 | 60.6 | 110.6 | 0.843 |
| B1_hub1p0_diag3p0 | 1 | 7.40 | 4.988 | 241.6 | 310.5 | 8.507 |
| B2_hub1p0_diag2p6 | 0 | 4.55 | 2.681 | 53.2 | 107.7 | 0.572 |
| B2_hub1p0_diag2p6 | 1 | 3.62 | 1.752 | 28.6 | 64.9 | 0.145 |
| G_baseline | 0 | 4.63 | 2.469 | 49.1 | 97.0 | 0.508 |
| G_baseline | 1 | 4.46 | 2.293 | 43.9 | 88.6 | 0.396 |
| LS_0p1nH | 0 | 5.02 | 2.848 | 62.0 | 116.5 | 0.848 |
| LS_0p1nH | 1 | 4.72 | 2.553 | 51.9 | 101.1 | 0.574 |
| LS_0p2nH | 0 | 5.53 | 3.354 | 83.5 | 146.6 | 1.592 |
| LS_0p2nH | 1 | 5.05 | 2.872 | 63.0 | 117.8 | 0.880 |
| LS_0p3nH | 0 | 6.24 | 4.055 | 125.5 | 199.9 | 3.473 |
| LS_0p3nH | 1 | 5.46 | 3.275 | 79.9 | 141.5 | 1.457 |

LS_* rows: ASSUMED package inductance in both diode states (sensitivity, not device data). For G the inductance RAISES anisotropy and leakage (unlike cell J).

Hexagon size for conversion: `hexagon_for_conversion_50_70.json` (script `phase28_hexagon_for_conversion.py`).
