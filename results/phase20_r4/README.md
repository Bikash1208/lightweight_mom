# Phase 20 R4 — CP-2 state library: metrics and caveats

MODEL_ONLY, NOT_HFSS_VALIDATED. The ~10 deg model uncertainty (R3) is **not** included in any number here. Library: `cp2_library_v1.npz` (+ manifest with sha256),
CP-2 octagon cell (period 15 mm, MADP-000907-14020, RT5880 1.575 mm), four one-hot states.

## How it was built (and what differs from the approved plan)
* Grid f 8–12 GHz step 0.25, theta 0–45 deg step 7.5, phi 0–330 deg step 30 (1241 points x 4 states) at mesh 60, **Floquet cutoff F = 2n** (see `FLOQUET_CONVERGENCE.md`).
  The approved grid was finer; the reduction was a cost decision made by the assistant and accepted afterwards by the user.
* Mesh correction: Richardson (p = 2) from mesh 60/90/120 at f every 0.5 GHz, theta in {0, 30, 45} deg (phi-averaged, interpolated in f and theta, added at every phi). Stored `uncertainty` = |p=1 vs p=2 extrapolation| + phi spread.
* Domain: f 8–12 GHz, theta <= 45 deg. theta > 30 deg is `EXTRAPOLATED_BEYOND_R3_TARGETS`. Points where a higher Floquet order propagates (e.g. theta = 45, phi = 0 above 11.71 GHz) are refused by the interpolator / registry
  (Wood guard band 3 %); outside the grid the answer is `NOT_AVAILABLE_FROM_CURRENT_DATA`.

## Interpolation error (item 4)
Held-out, seeded, 60 random off-grid points x 4 states (`validation_a.json`; library raw table vs direct mesh-60 / F = 120 solves):
theta <= 30 deg, r_RR: phase error median 0.09 deg, **95th percentile 1.81 deg (gate 3)**, max 8.7 deg (near a resonance); level error median 0.009 dB, **95th percentile 0.166 dB (gate 0.2)**, max 0.28 dB.
Result **GATE_PASS** (the dB margin is thin). theta > 30 deg (report only): phase 95th percentile 5.7 deg, level 1.8 dB, worst 178 deg / 11 dB (a resonance edge).

**Corrected table vs direct mesh-120 / F = 240 solves** (`validation_b.json`, 12 seeded points x 4 states): the mesh correction removes most of the mesh error but not all of it.
theta <= 30 deg: r_RR phase error median 1.5 deg, **95th percentile 12.1 deg** (max 15.8), level error median 0.29 dB, 95th percentile 0.96 dB (max 1.9 dB); theta > 30 deg: 95th percentile 15.2 deg / 2.0 dB.
The raw mesh-60 table is off by up to 11 deg / 1.4 dB at the same points, so the correction helps, but the residual is largest at f ~ 10.3–10.8 GHz (12–16 deg) next to the sharp resonance
that the 0.5 GHz sampling of the correction map cannot follow, and has a ~0.3–0.6 dB level offset at 8.5–9 GHz. This is **not** the brief's gate (the gate is the interpolation error above, which passes); it is the honest
residual of the mesh correction, comparable to the ~10 deg model uncertainty. The stored `uncertainty` (median 0.09, 95th percentile 0.23 in |C| units) was not checked against these residuals point by point. A finer correction grid (f step 0.25 GHz at 10–11 GHz) would reduce it; not done (cost).

## Target metrics at 10 GHz, worst over phi (item 5; `library_metrics.json`, corrected table)
R3 targets: lattice deviation <= 15 deg, step error <= 20 deg, min amplitude >= -1.5 dB, spread <= 1 dB, opposite-handed <= -15 dB.

| theta | lattice dev | step err | min amp | spread | opposite-handed |
|---|---|---|---|---|---|
| 0 | 10.7 | 18.2 | -1.16 | 0.24 | -15.2 |
| 7.5 | 11.2 | 18.9 | -1.15 | 0.21 | -15.3 |
| 15 | 12.4 | 20.5 | -1.15 | 0.23 | -14.6 |
| 22.5 | 14.4 | 22.8 | -1.19 | 0.27 | -12.4 |
| 30 | 17.0 | 26.0 | -1.35 | 0.38 | -10.0 |
| 37.5 (extrap.) | 19.9 | 29.4 | -1.74 | 0.67 | -7.9 |
| 45 (extrap.) | 23.5 | 34.0 | -2.33 | 1.03 | -5.9 |

* **Usable band.** All five targets are met only at the 10.0 GHz grid point for theta = 0 and 7.5 deg (the band is narrower than the 0.25 GHz grid, i.e. not resolved); nowhere for theta >= 15 deg. Dropping the opposite-handed target changes nothing (the ladder/amplitude limit it).
* **Opposite-handed term.** Converged (F = 2n) it is -15.2 dB at theta = 0, better than the R3 minimum-cutoff value (-13.6 dB): the R3 miss was partly numerical. It worsens quickly with angle (-10 dB at 30 deg), and it carries no Pancharatnam–Berry phase, so at array level it goes to a specular lobe, not the steered beam (R5 to verify at array level).
* Away from 10 GHz the cell is much worse: at 9.5 GHz the amplitude is about -3 dB, spread 1.4–1.7 dB, opposite-handed -5 dB (tolerance table below); the cell is a narrow-band design.

## LP interim — 1-bit codes of CP-2 (item 5, 10 GHz; `lp_codes_10GHz` in `library_metrics.json`)
co-LP code = states 1 vs 3, cross-LP code = states 2 vs 4. Phase difference between the two states (wanted 180 deg), minimum amplitude of the two, unwanted term re wanted, theta = 0 / 30 deg, phi = 0:
* TE co: 165° / 166°, -2.3 / -2.6 dB, unwanted -18.9 / -18.2 dB. TE cross: 177° / 177°, -1.0 / -1.3 dB, unwanted -11.4 / -8.3 dB.
* TM co: 176° / 169°, -2.2 / -2.3 dB, unwanted -16.8 / -16.6 dB. TM cross: 177° / 171°, -1.0 / -1.3 dB, unwanted -19.7 / -10.3 dB.
So the 1-bit codes give 155–177 deg phase difference (the co-LP code is the weaker, about 15 deg short at TE, with about -2.3 dB amplitude). 2-bit LP remains deferred (needs a multilayer / slot-coupled cell). Full phi = 0/90/180/270 data in the JSON.

## Tolerance corners (`tolerance_corners.json`; mesh 60, F = 120, nominal library correction added)
C_T +/-20 %, L_s 0 / 0.1 / 0.2 nH (ASSUMED_SENSITIVITY_NOT_DEVICE_DATA), f 9.5 / 10 / 10.5 GHz, theta 0 and 30 deg (phi 0, 45, 90). Selected rows (dev / step / opposite-handed):
* 10 GHz, theta 0: nominal 10.7 / 18.2 / -15.2; C_T -20 % 10.1 / 17.4 / -15.3; +20 % 12.0 / 20.8 / -14.2; L_s 0.1 nH 11.0 / 19.3 / -13.6; 0.2 nH 11.2 / 20.5 / -10.6. Phase ladder survives all corners; the opposite-handed level degrades with series inductance.
* 10.5 GHz, theta 0: nominal is good (7.1 / 14.2 / -18.6) but **any** corner breaks it (C_T -20 %: step 28 deg, amplitude spread 3.0 dB; L_s 0.1 nH: step 36 deg): a sharp resonance near 10.5 GHz, so the upper band edge is not robust.
* 9.5 GHz: amplitude -3.0 dB, spread 1.4 dB for all corners (the band's lower side is not usable).
* theta = 30 deg, 10 GHz: ladder within 17 deg lattice dev for all corners, opposite-handed -6.5 … -11 dB.

## Plots (`plots/`)
`state_response_vs_f.png` (phase / amplitude vs f per state), `ladder_deviation_vs_f.png`, `maps_theta_phi_10GHz.png` (opposite-handed vs theta, phi).

## Caveats
1: the Floquet-cutoff finding applies to **all** earlier RIS numbers computed at the minimum cutoff (R1–R3, and the committed `mom_reference.csv` files); see the Floquet-cutoff audit. 2: library is MODEL_ONLY until the user's HFSS runs exist; `src/ris/hfss_calibration.py`
compares and proposes offsets but applies a calibration only with explicit user approval. 3: isolated-cell values are reliable for regular gradients, not random codings at 15 mm pitch (R3 item 4); R5 must verify codings at array level.
