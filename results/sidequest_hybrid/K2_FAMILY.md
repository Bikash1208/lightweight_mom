# Stronger-leaking cell family + launcher, re-sized stacks - MODEL_ONLY (2026-11-05); no HFSS package built
Scripts: `sidequest_a0_direct.py` (A0 computed directly), `sidequest_family_scan.py` / `_eval.py` / `_design.py`, `sidequest_launcher.py`. Data: `k2_family_scan*.json`, `k2_family_lib.json`, `k2_family_summary.json`, `k2_launcher.json`, `k2_family_boards.json`, `k2_family_table.json`.

## 1. A0 computed directly (no A0 ~ dX^2 assumption)
First-order leakage of the -1 space harmonic of a +-dY/2 square sheet modulation on the real stack: alpha = P_rad / (2 S_sw), P_rad = 1/2 Re(Y_up)|J|^2/|Ybar + Y_up + Y_down(kx)|^2 at kx = k0 sin 25 deg, J = (2/pi)(dY/2)E0, S_sw = Poynting flow of the TM0 mode (grounded slab + sheet). No fitted constant. Check on the HFSS-confirmed cell (dY 3.13 mS, beta/k0 1.0731, RT5880 1.575): **0.836 Np/m against the HFSS back-out 0.844**. The MoM dY of the baseline cell is 2.63 mS (16 % below HFSS), so the MoM dY of every family is scaled by 3.13/2.63 (assumption, A0 x1.42).
Why not dX^2: A0 depends on |Ytot,-1| = |Ybar + Y_up + Y_down|: it grows when the mean loading cancels the slab's inductive Y_down at the radiating harmonic (leaky resonance) and falls with a thicker slab (higher beta, more stored SW power).

## 2. Cell families (10.5 GHz, 6.4 lattice, disc radius r, slit to r - 0.2, OFF = open slit, ON = 2.44 ohm R; MoM mesh 64)
| family | dY [mS] | X OFF / ON / dX [ohm] (lid) | beta/k0 | A0 [Np/m] (x baseline) | static levels m (solid ... strap) |
|---|---|---|---|---|---|
| baseline RT5880 1.575, r 3.0 (check) | 2.63 | 121 / 100 / 21 | 1.06 | 0.77 (x0.9; HFSS 0.844) | -1.67 ... +0.05 |
| RT5880 1.575, r 3.15 (gap 0.1 mm) | 8.31 | 199 / 120 / 79 | 1.13 | 2.86 (x3.4) | -1.44 ... -0.46, +1 (OFF) |
| RT5880 3.175, r 3.15 | 6.80 | 293 / 241 / 51 | 1.30 | 0.57 (x0.7) | -1.47 ... -0.43, +1 |
| **RO3006 1.27, r 2.2** | 3.03 | 293 / 185 / **108** | 1.20 | **4.77 (x5.6)** | -1.63 ... +0.21, +1 |
| RO3006 1.27, r 2.25 | 3.23 | 318 / 196 / 122 | 1.23 | 6.67 (x7.9) | -1.64 ... +0.21, +1 |

Result: **RT5880 does not reach x4.7** (best x3.4 at r 3.15, which is already a 0.1 mm gap between discs; 3.175 mm is worse, x0.7). **RO3006 1.27 reaches it (r 2.2 -> x5.6)**, with a dX swing x5 of the baseline. Static cells follow the same normalised levels (m = (Im Y - mid)/half-swing of the pair), so they carry the same A0; the range on RO3006 is -1.63 .. +0.21 (baseline -1.67 .. +0.05). A0 per state is not separable in the engine (alpha = A0 <m^2> over the two states): NA.
Sensitivity: A0 of RO3006 rises x4 from r 2.2 to 2.3 (MoM, uncorrected 3.4 -> 13.5) and beta/k0 from 1.20 to 1.36: a resonance, tolerance risk (r +-0.05 mm); r 2.45-2.5 gives beta/k0 about 1.9-2.0 (too tightly bound). dY includes the ON conductance; no HFSS check of any new cell.

## 3. Launcher (spectral vertical-pin model, grounded slab, thin pin)
Checks: lambda/4 monopole 39.6 ohm vs 36.5; **baseline probe launch 0.611 against the HFSS-based 0.627, R_in 42 ohm**; all launch values below are scaled by 0.627/0.611.
Direct radiation of the pin is intrinsic to a z-pin in a bare slab (a pin 1 mm above the surface gives the same split), so the lever is the loading around the probe:
| stack | bare around the probe (cells from 6.5 mm, as now) | cells loaded up to the probe (pad, upper bound) | R_in (t 6 mm, loaded) |
|---|---|---|---|
| RT5880 1.575 | 0.627 (direct 37 %) | 0.87 (r 3.0 loading) / 0.99 (r 3.15) | 56-57 ohm |
| RT5880 3.175 | 0.86 | 0.97 | NA |
| RO3006 1.27 | 0.71 | 0.96-0.98 | 57 ohm |

launch >= 0.8 therefore needs a loaded launch pad (cell-type annulus around the pin from small radius; no extra part, no diodes). The finite-radius pad is NOT computed (uniform loading = upper bound), so the real value lies between the bare and the loaded column. S11 from R_in alone (reactance not modelled): 57 ohm -> -23 dB, 36 ohm (bare RO3006) -> -16 dB; X is NA, so S11 <= -15 dB is not confirmed.

## 4. Boards (design (i), LP-V, g 2, mu 30, lam re-tuned to <= 514 diodes)
Gain / ring / back / SLL / isolation are calibrated to the one HFSS state-1 point: gain offset -1.78 dB, ring 36 % and back 39 % of the non-beam power, non-beam peaking 7.6, isolation floor factor 0.63. Edge = share of the launched SW arriving at the edge; ring = % of input at theta 60-90; back = % of input.
| board | launch | diodes | gain 1/2 [dBi] | SW edge | ring | back | SLL [dB] | isolation 1/2 [dB] |
|---|---|---|---|---|---|---|---|---|
| baseline (HFSS 9.1 dBi, ~10 dB) | 0.63 | 508 | 9.1 / 9.2 | 84 % | 32 | 35 | -0.8 | 10.0 / 10.3 (model-only 15.2 / 15.8) |
| RT5880 1.575 r 3.15, bare | 0.63 | 510 | 13.7 / 12.8 | 55 % | 26 | 28 | -5.5 | 13.2 / 12.5 |
| RT5880 1.575 r 3.15, pad | 0.99 | 510 | 15.7 / 14.8 | 55 % | 20 | 22 | -8.6 | 14.2 / 13.5 |
| RT5880 3.175 r 3.15, bare / pad | 0.86 / 0.97 | 514 | 9.3 / 7.9 ; 9.8 / 8.4 | 86 % | 31 | 34 | +0.3 / -0.2 | 10.2 / 9.3 ; 10.6 / 9.8 |
| **RO3006 r 2.2, bare** | 0.71 | 512 | 15.1 / 16.1 | 44 % | 21 | 23 | -8.6 | 13.8 / 14.6 |
| **RO3006 r 2.2, pad (launch >= 0.8)** | 0.96 | 512 | **16.5 / 17.4** | 44 % | 16 | 18 | -11.1 | 14.3 / 15.1 |
| RO3006 r 2.2, pad + rim 10 dB (extra part) | 0.96 | 502 | 16.1 / 16.9 | 7 % | 4 | 4 | -11.6 | 15.3 / 15.5 |
| RO3006 r 2.25, pad | 0.98 | 504 | 17.3 / 18.1 | 32 % | 12 | 13 | -13.2 | 13.9 / 15.2 |
| RO3006 r 2.25, pad + rim 10 dB | 0.98 | 508 | 17.0 / 17.8 | 5 % | 3 | 3 | -13.2 | 13.7 / 16.2 |

Edge power: RO3006 r 2.2 gives 44 % (the <= 15 % target is not met without the rim; with the rim 7 %); r 2.25 gives 32 %. The rim is an "extra part" (absorber, 12.8 mm wide, 10 dB) and costs about 0.4 dB.

## 5. Isolation gap (HFSS ~10 dB against the model's 15 dB)
The model's isolation counts only the hologram's own leakage of the other beam (static compromise cells; -4.1 dBi, 15 dB below the beam). HFSS also contains the non-hologram radiation: the direct probe (37 % of the input) and the edge diffraction (52 %) fill every direction at roughly 0.63 of the isotropic level, about -2.4 dBi at theta 25 deg, the same level as the other beam. With that floor (one factor, fitted so the baseline gives 10 dB) the model's 15 dB becomes 10 dB: isolation = gain - 10 log(other_beam_leak + NB x 0.63), NB = (1 - launch) + launch x edge. The new boards therefore gain isolation only through gain and a smaller NB (14-15.5 dB), not through the mu term. HFSS state 2 is not available yet: the factor is a one-point fit (state 1 only).

## 6. Caveats
* Everything is MODEL_ONLY: the new cells are not HFSS-checked (MoM normal-incidence extraction, dY corrected by the baseline ratio); the launcher is a thin-pin spectral model checked against one HFSS number; the pad is an upper bound; the ON diode loss (AA0) and the background loss are the baseline values (background scaled by sqrt(tan d) for RO3006).
* RO3006 needs a new HFSS probe calibration (the 0.627 and -15 dB numbers are RT5880).
* Leaf track closed; failure recorded in `results/sidequest_leaf/LEAF_CELL.md` section 1b.
