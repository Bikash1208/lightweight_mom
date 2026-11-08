# CELL_SCAN: slotted disc + rim-notch diode cell, RT5880 3.175 mm, p = 11.6 mm (server 1)

**MODEL_ONLY / NOT_HFSS_VALIDATED.** 5.8 GHz, mesh 64, F = 2n = 128, Phase 26 slotted circular patch + 1.00 mm rim notch at the +u slot end bridged by a lumped sheet (1.00 across x 0.35 along the rim, centred in the rim). `SlotDiscNotch` (duck-typed, `scripts/research/sub6_slotdisc_lib.py`; no frozen code touched) equals `SlottedCircularPatch` without a notch (`tests/test_sub6_slotdisc_geometry.py`). Job `jobs/sub6-slotdisc-cell-scan.md` with Amendments 1-3. RO3006 (server 2 / AEM-LAB1) is NOT in this file; TMM6 cancelled.

## Gates
| gate | state |
|---|---|
| ANCHOR_REPRODUCED | **ANCHOR_REPRODUCED**: scaled hex150 majority cell (r 5.45, 8.72 x 2.91 slot, p 11.63, closed, RT5880 3.175) α(25°) 1.64 Np/m vs 3.14 (10.54 GHz) x 5.8/10.54 = 1.73, ratio 0.947. β/k0 1.087 vs 1.068, Y_par/Y_perp 5.31/10.4 mS vs 5.38/10.7 mS (`anchor.json`) |
| LEAKAGE (per row, 1 Np/m at 25°) | OK_BOTH **0**, OK_ON_ONLY **45**, OK_OFF_ONLY **1**, LOW **34** of 80 |
| MESH_64_96 | see the end of this file (3 best rows by max(α_ON, α_OFF)) |
| RESONANCE (5.0-6.6 GHz) | **PROXY** (sign flip of Im Y at normal incidence at 5.0 vs 6.6 GHz): 16 rows flagged, 14 of them at L/2r = 0.9 |

## Findings
- **The beam is carried by the ON state.** α_OFF has median 0.03 Np/m (max 1.74 in one row); α_ON has median 1.30, max 4.05 (25 rows ≥ 2 Np/m, 45 rows ≥ 1). In OFF_BAR50 the anisotropy Δ is small (|Δ_OFF| ≤ 6 mS, often < 1 mS) because the empty notch (Δ_OPEN +5.7 … +14.4 mS over the grid) and the closed-slot cell (Δ_NOCUT negative, e.g. −5.0 mS in the anchor-nearest row) have opposite sign and the OFF load (−j161 Ω + 6.7 Ω) sits in between: OFF is a near-isotropic sheet.
- **ρ = Δ_ON/Δ_OFF is ill-conditioned** (Δ_OFF ≈ 0): |ρ| 5.8-53 in the 12 best rows (1.3-295 over all 80), not the −1.4…−3.5 of the cell scan. Against the open notch, ρ = Δ_ON/Δ_OPEN is −1.1…−4: 12 rows lie in [−3.5, −1.3] with α_ON 1.8-4.1 (table 2). The OFF_OPEN state is a different circuit from OFF_BAR50 (no diode), so this ρ is NOT the diode ρ; reported for comparison only.
- α_ON by L/2r: 0.53 (0.6), 1.30 (0.7), 2.25 (0.8), 1.64 (0.9); by gap and by W it is flat (1.4-1.5): the slot length (rim) sets the leakage, the gap and W do not.
- **A (antisymmetric part)**: zero along the slot (A_par = 0 by the mirror about the slot axis); across the slot A_perp is non-zero because the notch breaks the other mirror: OFF +0.9 … +1.3 mS, ON −1.7 … −0.4 mS in the 12 best rows; ON−OFF change of A_perp −3.0 … −1.3 mS there (−14.9 … +3.2 over all 80). A_perp(ON) is 3-10 % of |Δ_ON| in these rows.
- NOCUT α: 0.0-4.1 over the grid (the closed slot; the anchor row repeats the anchor).

## Table 1: best leakage (top 12 of 80 by max(α_ON, α_OFF)); all 80 rows in `rows_summary.csv`
| row (gap_LoverR2_W) | α_OFF | α_ON | ON/OFF | Δ_OFF mS | Δ_ON mS | ρ=Δ_ON/Δ_OFF | ρ=Δ_ON/Δ_OPEN | A_perp OFF→ON mS | β/k0 OFF/ON | leakage | res |
|---|---|---|---|---|---|---|---|---|---|---|---|
| g1_q0.9_W2.9 | 0.014 | 4.05 | 293 | -0.61 | -13.47 | 21.4∠12° | -2.26 | 0.95→-0.51 | 1.050/1.140 | OK_ON_ONLY |  |
| g1_q0.9_W2.4 | 0.002 | 4.05 | 1956 | -0.25 | -10.78 | 37.9∠28° | -1.78 | 0.91→-0.44 | 1.049/1.118 | OK_ON_ONLY |  |
| g0.8_q0.8_W3.6 | 0.003 | 3.24 | 1049 | -0.26 | -9.66 | 32.6∠27° | -1.43 | 0.98→-0.51 | 1.054/1.119 | OK_ON_ONLY |  |
| g0.73_q0.8_W3.6 | 0.001 | 3.16 | 3998 | 0.13 | -10.35 | 52.0∠129° | -1.50 | 1.06→-0.59 | 1.054/1.126 | OK_ON_ONLY |  |
| g0.8_q0.9_W1.6 | 0.007 | 3.07 | 456 | 0.42 | -8.34 | 18.8∠160° | -1.29 | 0.88→-0.37 | 1.050/1.107 | OK_ON_ONLY |  |
| g0.6_q0.9_W1.6 | 0.027 | 2.96 | 108 | 0.81 | -10.56 | 12.8∠167° | -1.54 | 1.01→-0.54 | 1.052/1.132 | OK_ON_ONLY |  |
| g0.73_q0.9_W1.6 | 0.178 | 2.91 | 16 | 3.11 | -17.77 | 5.8∠170° | -2.68 | 1.29→-1.72 | 1.040/1.184 | OK_ON_ONLY | PROXY |
| g0.6_q0.8_W3.6 | 0.016 | 2.91 | 177 | -0.53 | -10.45 | 18.8∠14° | -1.47 | 0.99→-0.54 | 1.060/1.131 | OK_ON_ONLY |  |
| g0.73_q0.8_W2.9 | 0.011 | 2.91 | 265 | 0.53 | -8.30 | 15.3∠163° | -1.17 | 1.02→-0.52 | 1.051/1.109 | OK_ON_ONLY |  |
| g0.8_q0.8_W2.9 | 0.001 | 2.90 | 2672 | 0.16 | -7.77 | 36.9∠139° | -1.11 | 0.95→-0.46 | 1.052/1.103 | OK_ON_ONLY |  |
| g1_q0.9_W3.6 | 0.100 | 2.84 | 28 | -1.48 | -23.04 | 15.5∠4° | -4.01 | 1.03→-0.78 | 1.054/1.211 | OK_ON_ONLY | PROXY |
| g0.6_q0.8_W2.9 | 0.000 | 2.71 | 9904 | -0.07 | -8.38 | 52.6∠61° | -1.14 | 0.95→-0.48 | 1.056/1.113 | OK_ON_ONLY |  |

## Table 2: rows with ρ(ON/OPEN) in [−3.5, −1.3] (top 12 by α_ON)
| row (gap_LoverR2_W) | α_OFF | α_ON | ON/OFF | Δ_OFF mS | Δ_ON mS | ρ=Δ_ON/Δ_OFF | ρ=Δ_ON/Δ_OPEN | A_perp OFF→ON mS | β/k0 OFF/ON | leakage | res |
|---|---|---|---|---|---|---|---|---|---|---|---|
| g1_q0.9_W2.9 | 0.014 | 4.05 | 293 | -0.61 | -13.47 | 21.4∠12° | -2.26 | 0.95→-0.51 | 1.050/1.140 | OK_ON_ONLY |  |
| g1_q0.9_W2.4 | 0.002 | 4.05 | 1956 | -0.25 | -10.78 | 37.9∠28° | -1.78 | 0.91→-0.44 | 1.049/1.118 | OK_ON_ONLY |  |
| g0.8_q0.8_W3.6 | 0.003 | 3.24 | 1049 | -0.26 | -9.66 | 32.6∠27° | -1.43 | 0.98→-0.51 | 1.054/1.119 | OK_ON_ONLY |  |
| g0.73_q0.8_W3.6 | 0.001 | 3.16 | 3998 | 0.13 | -10.35 | 52.0∠129° | -1.50 | 1.06→-0.59 | 1.054/1.126 | OK_ON_ONLY |  |
| g0.6_q0.9_W1.6 | 0.027 | 2.96 | 108 | 0.81 | -10.56 | 12.8∠167° | -1.54 | 1.01→-0.54 | 1.052/1.132 | OK_ON_ONLY |  |
| g0.73_q0.9_W1.6 | 0.178 | 2.91 | 16 | 3.11 | -17.77 | 5.8∠170° | -2.68 | 1.29→-1.72 | 1.040/1.184 | OK_ON_ONLY | PROXY |
| g0.6_q0.8_W3.6 | 0.016 | 2.91 | 177 | -0.53 | -10.45 | 18.8∠14° | -1.47 | 0.99→-0.54 | 1.060/1.131 | OK_ON_ONLY |  |
| g0.8_q0.9_W2.4 | 0.042 | 2.65 | 63 | -0.92 | -21.78 | 23.3∠8° | -3.51 | 1.04→-0.78 | 1.056/1.208 | OK_ON_ONLY | PROXY |
| g0.4_q0.8_W2.4 | 0.003 | 2.41 | 708 | 0.23 | -13.85 | 47.1∠140° | -1.77 | 1.17→-0.79 | 1.062/1.169 | OK_ON_ONLY |  |
| g0.4_q0.8_W2.9 | 0.000 | 2.28 | 7457 | -0.07 | -16.53 | 81.7∠67° | -2.14 | 1.21→-0.89 | 1.064/1.189 | OK_ON_ONLY |  |
| g0.4_q0.9_W1.6 | 0.001 | 2.19 | 1525 | 0.15 | -18.21 | 76.8∠126° | -2.54 | 1.09→-0.74 | 1.062/1.200 | OK_ON_ONLY |  |
| g0.4_q0.8_W3.6 | 0.036 | 1.80 | 50 | -0.72 | -23.65 | 31.6∠13° | -3.12 | 1.28→-1.17 | 1.069/1.235 | OK_ON_ONLY | PROXY |

## Flags and assumptions
- Notch 1.00 mm = 5.5 mesh cells at mesh 64 (p 11.6, dx 0.181 mm); load/notch edges are off the grid, so the effective width is quantised: the mesh 64/96 check below is the only test of this.
- Rim rule r − L/2 ≥ 0.35 mm kept (the diode 0.35 mm fits exactly at L/2r = 0.9: rim 0.54-0.56 mm for p 11.6). The corner-rim skip rule removed no RT5880 row.
- OFF_BAR50 = (0.15 pF ‖ 5 kΩ) + 0.6 nH, ON_BAR50 = 3.0 Ω + 0.6 nH; OFF_BAR50_HI and ON_BAR50_MAX not run (later, 5 best rows).
- α uses only Im(Y) (as Phase 26), mean-sheet β self-consistent with the Bloch tensor at kt = β, linear Y(kt) between 2 nodes per row (kt/k0 1.06 and a second near 1.07-1.11).
- LEAKAGE gate as corrected by the user: OK_BOTH | OK_ON_ONLY | OK_OFF_ONLY | LOW at 1 Np/m, 25°.
- Not HFSS-validated; the 10.5 GHz cell scan HFSS comparison problem (single-diode switching) applies here too.

## MESH_64_96 (MODEL_ONLY / NOT_HFSS_VALIDATED)
3 best rows by max(alpha_ON, alpha_OFF); mesh 64/F128 vs mesh 96/F192; data `mesh_check_RT5880.json`. The job text gives no numeric threshold; the gate state below is my reading (flagged).

| row | state | alpha 64 | alpha 96 | d% | beta/k0 64 | beta/k0 96 | d% | abs(Delta) mS 64 | 96 | d% |
|---|---|---|---|---|---|---|---|---|---|---|
| g1_q0.9_W2.9 | ON | 4.05 | 4.05 | 0.0 | 1.1401 | 1.1360 | -0.36 | 13.5 | 12.9 | -4.1 |
| g1_q0.9_W2.9 | OFF | 0.0138 | 0.00306 | -77.9 | 1.0502 | 1.0490 | -0.11 | 0.629 | 0.333 | -47.0 |
| g1_q0.9_W2.4 | ON | 4.05 | 3.88 | -4.1 | 1.1184 | 1.1103 | -0.73 | 10.8 | 9.73 | -9.9 |
| g1_q0.9_W2.4 | OFF | 0.00207 | 0.000619 | -70.1 | 1.0488 | 1.0474 | -0.13 | 0.285 | 0.195 | -31.5 |
| g0.8_q0.8_W3.6 | ON | 3.24 | 3.17 | -2.3 | 1.1185 | 1.1139 | -0.41 | 9.67 | 9.12 | -5.7 |
| g0.8_q0.8_W3.6 | OFF | 0.00309 | 0.00325 | +5.2 | 1.0541 | 1.0541 | 0.00 | 0.297 | 0.301 | +1.5 |

MESH_64_96 = CONVERGED_ON_STATE (alpha_ON within 4.1%, beta within 0.73%, |Delta_ON| within 9.9%); OFF_NOT_RESOLVED (alpha_OFF ~0.001-0.014 Np/m, |Delta_OFF| < 0.7 mS, relative change up to 78% but absolute change < 0.011 Np/m). Conclusion unchanged: beam carried by ON, OFF nearly isotropic. Mesh 96 lowers |Delta_ON| 4-10%, so rho values carry that uncertainty. Notch quantisation (5.5 cells at 64) not separately isolated.
