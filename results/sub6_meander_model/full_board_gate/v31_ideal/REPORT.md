# Full-board v3.1 premise test — A-only-IDEAL

## STATUS

COMPLETE. One new HFSS 2026 R1 solve at 5.68 GHz converged. This is an idealized, non-buildable upper bound. It suppresses the LHCP B peak by 10.59 dB but also reduces the RHCP A peak by 0.71 dB relative to the saved BOTH reference.

## GATES

The source CSV has 516 sites. Exactly **234** slotted B cells satisfying `beam=B`, `x<=0`, `rho>=29 mm` were replaced by plain r=5.45 mm discs. All 16 gate discs remain plain and ungrounded; there are **zero vias**. A row-by-row manifest comparison to BOTH found the same site coordinates, beam labels, gate labels, rotations, and no-via state outside the 234 intended replacements. STEP literal audit: zero inexact values. AEDT read back 518 PEC metal objects, the same materials, port, radiation boundary, absorber, and 5.68 GHz setup settings. The independent mesh differs because the copper geometry changed.

HFSS's first adaptive pass used **11.0 GB** maximum memory/process (298,089 solved tetrahedra; 46 granted cores). The solve converged in 6 passes, final max delta-S 0.00200, 52.7 GB peak memory/process, 22m45s solve wall time.

## DATA

| State | RHCP A peak | LHCP B peak | A drop vs BOTH | B drop vs BOTH | RHCP SLL | LHCP SLL | RHCP wide | LHCP wide | S11 | Accepted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BOTH reference | 9.03 dBi (32°, 0°) | 9.00 dBi (36°, 180°) | — | — | -8.43 dB | -9.52 dB | -4.85 dBi | -3.34 dBi | -16.76 dB | 97.89% |
| A-only-IDEAL | 8.32 dBi (32°, 0°) | -1.59 dBi (44°, 176°) | **0.71 dB** | **10.59 dB** | -8.22 dB | -0.52 dB | -4.11 dBi | -3.26 dBi | -17.22 dB | 98.10% |

At the *fixed design directions*, A (30°,0°) drops 0.56 dB and B (30°,180°) drops 13.98 dB. Peak drops use each hand's global above-horizon peak, as requested. SLL uses the same first-null definition as the v3.1 report; wide-angle is maximum gain at 60°<=theta<=90°.

## FLAGS

The suppressed LHCP pattern has no clear designed B beam; its -0.52 dB SLL is relative to a residual peak at (44°,176°) and does not represent a conventional beam sidelobe. HFSS gain is normalized to accepted power, while S11 and accepted power are reported separately. The solver read-back/profile indicates direct DCS despite the iterative setting request inherited from BOTH. No other state was built or solved.

## FILES

`metrics.json` holds the unrounded comparison; `geometry.json` records all 516 site choices and 234 replacements; `solve.prof`, `solve_status.json`, and `A_ONLY_IDEAL.ffd` preserve solver and far-field evidence. The STEP and AEDT project remain in private local `hfss_work/full_board_gate_v31_ideal/` and are not pushed.
