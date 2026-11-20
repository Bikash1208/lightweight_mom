# Full-board gate validation v3.1 — HFSS 2026 R1

## STATUS

COMPLETE. The full 516-site board was built and solved as BOTH, A-only, and B-only, one at a time. Both one-beam gate states **FAIL** the >=10 dB suppression and no-survivor-loss requirement. All three projects validated and reached normal adaptive convergence.

## GATES

- **First-pass resource gate: PASS.** HFSS reported 11.7 GB maximum memory/process for the first adaptive pass (323,608 solved tetrahedra; 1,873,383 initial mesh tetrahedra; 46 granted cores). The host has 256 GB physical RAM and ample free RAM. Measured continuation maxima were 44.8 GB (BOTH), 69.1 GB (A-only), and 51.3 GB (B-only).
- **BOTH beams: PASS.** RHCP peaks at (32°, 0°), LHCP at (36°, 180°); no rotation correction was applied. These are respectively 2° and 6° from the specified theta=30° targets.
- **A-only: FAIL.** LHCP peak drops only **0.15 dB**, while the surviving RHCP peak changes **-1.10 dB** against BOTH.
- **B-only: FAIL.** RHCP peak drops only **0.80 dB**, while the surviving LHCP peak changes **-0.11 dB** against BOTH.

## DATA

| State | RHCP peak dBi (theta, phi) | LHCP peak dBi (theta, phi) | RHCP SLL dB | LHCP SLL dB | RHCP wide dBi | LHCP wide dBi | S11 dB | Accepted power |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BOTH | 9.03 (32, 0) | 9.00 (36, 180) | -8.43 | -9.52 | -4.85 | -3.34 | -16.76 | 97.89% |
| A-only | 7.94 (32, 0) | 8.85 (36, 180) | -6.57 | -7.71 | -4.13 | -2.67 | -17.09 | 98.04% |
| B-only | 8.23 (32, 0) | 8.89 (36, 180) | -7.45 | -8.13 | -3.28 | -4.02 | -14.73 | 96.63% |

Gain is HFSS `dB(GainRHCP)` / `dB(GainLHCP)`; accepted power is `1-|S11|²`. SLL is the highest same-hand lobe above the horizon outside the direction-dependent first null (fallback: 2× mean HPBW). Wide-angle level is the maximum same-hand gain over 60°<=theta<=90°. Far-field grid: theta/phi 2° at 5.68 GHz. The untouched 30-mm dither assignment is present in the source CSV and all 516 site rotations.

| State | Adaptive passes | Maximum HFSS memory | Solve wall time | Requested/granted cores | Final max delta-S |
|---|---:|---:|---:|---:|---:|
| BOTH | 6 | 44.8 GB | 950 s continuation + 378 s first pass | 46/46 | 0.00248 |
| A-only | 6 | 69.1 GB | 1547 s | 46/46 | 0.00329 |
| B-only | 6 | 51.3 GB | 1395 s | 46/46 | 0.00192 |

The project uses 516 sites: 500 slotted anchor discs and 16 plain gate discs. A-only grounds eight B gate discs through 0.3-mm-square PEC vias; B-only grounds eight A gate discs. The centre feed cell is excluded. The STEP audit found zero inexact literals; AEDT read back 518 PEC metal objects (516 discs, ground, pin), plus RT5880 substrate, RT5880 bias layer, AIR, PORT, and one absorber. One lumped port, one radiation boundary, and one infinite sphere validated.

## FLAGS

- The physical model inherits the R100 testbed feed (pin radius 0.635 mm, ground opening radius 2.05 mm, port at z=-3.175 mm), RT5880 top slab 3.175 mm, RT5880 bias layer 0.254 mm, and idealized absorber. The 30-mm absorber occupies radius 156–186 mm, starting beyond all outer-cell copper; this radius is an implementation choice because the supplied site cutoff is at centre radius 150 mm.
- The CSV's 30-mm dither is represented by per-site beam assignments and copper rotations. No separate radial alpha profile was imposed; the 500 anchor cells all have the same copper geometry.
- An iterative solver setting was requested for the first pass. AEDT accepted the property update but read back `Direct Solver` on reopen, and its profile says `Type: DCS`. Results are from the direct solver. DDM was not used.
- The independent Jones model's -12.9/-14.3 dB suppression checks are **MODEL_ONLY** and disagree with these HFSS gate results. These 16 grounded plain-disc gates do not switch the beam-generating cells off in this full-board electromagnetic model.
- Obsolete layouts were excluded from the commit and no results from them are reported. Automatic approval review rejected recursive removal of their private local work folders; the obsolete folders remain outside version control and were not used here. All their own running PIDs were stopped before v3.1.

## FILES

`metrics_summary.json` has every raw peak, drop, SLL, wide-angle and accepted-power value. `*.ffd`, `*.prof`, `*_solve_status.json`, geometry manifests and the exact input CSV preserve the HFSS evidence. The build/solve/extract/analysis scripts are under `scripts/research/full_board_gate_v31_*.py`; AEDT project files and STEP intermediates remain in the private `hfss_work/full_board_gate_v31/` folder and are not pushed.
