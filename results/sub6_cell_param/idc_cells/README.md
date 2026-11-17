# Split-disc / interdigital (IDC) cell probes, CORRECTED tool (`bar90_tune_work/tune_idc2.py`, local-only; see CORRECTION.md)
HFSS eigenmode, periodic cell, C40B floor (PEC lid 40 mm), RT5880 3.175 mm + same-material bias layer, pitch 11.6 mm, disc r 5.6, slot 0.4 mm, 0.3 mm minimum feature (checked). Rows `MN_R5p6_W0p4_N<teeth>_T1p2_D<depth>_E_OFF_<direction>_<phase>.json`; "across" = propagation along x (along the meander axis), "along" = along y; beta/k0 = (phase/80.8)*(5.8/f). Floor MN_MINF 3.0 GHz. Build report records substrate eps_r 2.2; log prints DISC copper area (88.2 mm2 for N2 = two halves; 94.0 mm2 for the straight split).
## Edge-to-edge 2-tooth IDC (N2, D3.0, slot rim to rim, two equal halves, gap centreline 25.8 mm)
Lowest bound branch (mode 1; along 88 also has a second bound mode at 6.079 GHz, beta 1.039; mode 1 continues the slope 0.052 GHz/deg like the static cell 0.055): across 88 deg 5.820 GHz (1.085, Q 2439), 100 deg 6.276 GHz (1.144), slope 0.038 GHz/deg; along 80 deg 5.478 GHz (1.048), 88 deg 5.890 GHz (1.073), slope 0.052. Interpolated to 5.8 GHz: across about 1.083, along about 1.0675; mean 1.075, delta beta about 0.015. Gate: branch rising PASS, no flat branch below 6.5 GHz PASS, mean beta/k0 1.05-1.15 PASS, delta beta >= 0.04 FAIL (0.015; static slot cell 0.022, closed 2-tooth meander 0.025 with opposite sign). Hypothesis: splitting the disc lowered beta for propagation across the split (gap series capacitance reduces the loading), flipping the sign of delta beta; finger count/length is then the contrast knob. Not yet tested.
## Straight split (N0, gap 11.2 mm, two halves 47.0 mm2 each, DISC area 94.0 mm2 verified)
Reruns done one at a time, all converged (first-run SOLVER_FAILUREs were transient). Lowest bound branch: across 88 deg 5.877 GHz (1.075, Q 2929), 100 deg 6.476 GHz (1.108), slope 0.050; along 80 deg 5.522 GHz (1.040), 88 deg 6.035 GHz (1.047), slope 0.064. Interpolated to 5.8 GHz: across about 1.070, along about 1.044; mean 1.057, delta beta about 0.026. Gate: rising PASS, no flat branch PASS, mean PASS, delta beta >= 0.04 FAIL (0.026, static level).
## Reading
Delta beta: straight split 0.026, edge-to-edge N2 0.015, static slot cell 0.022, closed N2 meander 0.025 (opposite sign). More interdigital finger capacitance lowered the contrast, so fingers are NOT a contrast knob in the direction hoped; the split-disc family at 11.6 mm pitch gives static-level contrast (0.015-0.026). Not yet tested: slot width (gap capacitance) as a knob (0.3 mm minimum), larger or smaller pitch.
## Gap centreline length (one convention: leg and horizontal centrelines plus end runs to the rim)
straight split 11.2 mm; edge-to-edge N1 D2.4 18.1 mm; edge-to-edge N2 D3.0 25.8 mm; closed 2-tooth about 21.8 mm (in band); closed N3 D2.2 25.4 mm (flat near 5.9 GHz); closed N2 D4.0 26.8 mm (flat near 5.8 GHz); closed N3 D3.0 31 mm (flat near 5.3 GHz).
## Fabrication check (0.3 mm minimum feature)
All cells pass: halves at least 0.4 mm apart, 0.3 mm morphological opening keeps >= 99.9 % of each half, cell-to-cell gap 0.4 mm.

## 14.5 mm pitch test (MN_R7p1_W0p4_N2_T1p2_D3p0_E_P14p5, RT5880 eps 2.2, DISC area 146.8 mm2): FAIL
Stored beta_over_k0, below_light_line and side fields in these four rows use the 11.6 mm factor (80.8 deg) and are wrong; true beta/k0 = stored x 0.800 (101 deg per k0 at 14.5 mm).
| Mode | Direction | f at two phases (GHz) | slope GHz/deg | beta/k0 at 5.8 GHz |
|---|---|---|---|---|
| m1 (Q ~1.6-1.9e3) | along 100/110 | 5.23 -> 5.32 | 0.009 flat | - |
| m1 | across 105/115 | 5.38 -> 5.57 | 0.019 flat | - |
| m2 | along | 5.49 -> 5.93 | 0.044 (expected) | ~1.059 |
| m2 | across | 5.83 -> 6.02 | 0.019 (hybridised) | ~1.02 |
Gate failed: flat branch 5.2-5.6 GHz in band, mean beta/k0 ~1.04 (<1.05). The along-vs-across difference is not reported as contrast (two points, hybridised mode). Q>1e5 modes are box/lid modes. Hypothesis (unverified): the split-half resonance moves into band as the disc grows (1/r from r 5.6: ~6.8 GHz -> 5.4 GHz at r 7.1); radius ceiling ~5.9 mm, pitch <= ~12.1 mm. Physical modes at 6.1-7.4 GHz exist in the 11.6 mm N2_E rows but are not identified, so no numeric ceiling is claimed.

## Rounded-rectangle outline + 2-tooth edge-to-edge split (MN_R5p6_W0p4_N2_T1p2_D3p0_E_O0p75; x semi-axis 5.6, y 4.2, n=8; copper 81.6 mm2, RT5880 eps 2.2): PROVISIONAL
Per-mode numbers are beta/k0 at each mode's own frequency; beta at 5.8 GHz comes only from interpolating an identified branch.
| Direction | Branch through 5.8 GHz | slope GHz/deg | beta/k0 at 5.8 (interpolated) |
|---|---|---|---|
| along y (phase on y) | 5.537 GHz @80 (Q 8.7e3) -> 6.053 @88 (Q 6.8e3) | 0.065 clean | ~1.041 |
| across x (phase on x) | 5.628 @88 (Q 1.8e3) -> 6.026 @100 (Q 1.55e3) | 0.033, hybridised | ~1.153 |
If the x branch holds: difference ~0.11, mean ~1.10 (round N2_E: x 1.083, y 1.0675). Gate clause "no flat branch below 6.5 GHz" FAILS as written: low-Q/flat modes along y 5.19-5.33 and 6.33-6.35 GHz, across x 6.22-6.43 GHz (second bound branch just above the x branch). Round N2_E 11.6 mm rows have no Q~1.9e3 mode at 5.2-5.3 GHz (lowest physical modes 5.48 and 5.82), so these low-Q modes are new with the outline (also seen at 14.5 mm/r7.1); the 1/r hypothesis in the 14.5 mm section is NOT refuted. Pending: across 65/90/95 solves to confirm the x branch.
Other agent (MODEL_ONLY, inversion of HFSS beta, theta 30): N2_E round alpha 0.0997 Np/m (conversion 2.95 % at R150), static ST 0.138 Np/m (4.06 %); 50 % needs 2.31 Np/m. alpha does not scale as dBeta^2 (ratio 0.72, not 0.46); earlier "9/4/12 %" estimates withdrawn.
