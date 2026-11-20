# Full-board v3.1 ceiling test — HFSS 2026 R1

## STATUS

COMPLETE. The B-ALL-OFF-IDEAL state was solved at 5.68 GHz. Its global residual LHCP peak falls only **10.91 dB** below BOTH, so the requested conditional ALL-PLAIN state was also solved. Both new solves validated and converged; no other state was built or solved.

## GATES

All 258 B-labelled sites are plain discs in B-ALL-OFF-IDEAL, including eight B gate sites already plain in BOTH. Thus **250 slotted discs changed** relative to BOTH. All 16 gate discs remain plain and ungrounded; no vias exist. ALL-PLAIN has 516 plain discs, of which 500 replace slots relative to BOTH. Manifest readback found identical site coordinates, beam/gate labels, rotations and zero-via state everywhere else; STEP literal audits found zero inexact literals. Both projects used the same 5.68 GHz setup and automatic mesh settings as BOTH. Their adaptive meshes differ because their copper geometry differs.

| New state | Adaptive passes | Final max delta-S | Peak HFSS memory/process | Solve wall time | Cores granted |
|---|---:|---:|---:|---:|---:|
| B-ALL-OFF-IDEAL | 6 | 0.00136 | 42.5 GB | 1220 s | 46 |
| ALL-PLAIN | 6 | 0.00288 | 63.0 GB | 1214 s | 46 |

## DATA

Gain at the **same two directions** in dBi, read directly from each converged HFSS project:

| State | RHCP (32°,0°) | LHCP (32°,0°) | RHCP (36°,180°) | LHCP (36°,180°) |
|---|---:|---:|---:|---:|
| BOTH reference | +9.03 | -19.25 | -2.33 | +9.00 |
| A-only-IDEAL, 234 slot removals | +8.32 | -26.94 | -1.74 | -2.80 |
| B-ALL-OFF-IDEAL, 250 new slot removals / 258 B sites plain | +8.60 | -14.50 | -3.43 | -3.48 |
| ALL-PLAIN, 500 slot removals | -18.10 | -18.11 | -11.30 | -11.60 |

| State | Global LHCP residual peak (dBi, direction) | LHCP peak drop vs BOTH | Global RHCP peak (dBi, direction) | RHCP peak drop vs BOTH | SLL RHCP / LHCP | S11 | Accepted |
|---|---:|---:|---:|---:|---:|---:|---:|
| A-only-IDEAL, 234 | -1.59 (44°,176°) | 10.59 dB | +8.32 (32°,0°) | 0.71 dB | -8.22 / -0.52 dB | -17.22 dB | 98.10% |
| B-ALL-OFF-IDEAL, 258 B plain | **-1.91 (44°,144°)** | **10.91 dB** | +8.60 (32°,0°) | 0.43 dB | -9.15 / -0.91 dB | -17.59 dB | 98.26% |
| ALL-PLAIN, 516 | -3.78 (48°,84°) | 12.77 dB | -3.87 (48°,278°) | 12.91 dB | -0.07 / -0.03 dB | -17.84 dB | 98.36% |

The B-ALL-OFF global LHCP peak improves only **0.32 dB** over the 234-cell ideal state. At the original B direction (36°,180°), the all-plain LHCP floor is **-11.60 dBi**, 8.11 dB below B-ALL-OFF's -3.48 dBi there. That difference is not a source-by-source decomposition; A-cell scattering, feed radiation and coherent interference coexist in B-ALL-OFF.

## FLAGS

The suppressed-hand residuals have no well-defined intended B main beam. Their SLL values are computed from their global residual peaks with the same first-null algorithm as the v3.1 report and are descriptive only. HFSS `GainRHCP` and `GainLHCP` are relative to accepted power; S11 and accepted fraction are separately reported. The solver profiles indicate direct DCS despite the iterative setting request. No symmetry plane, frequency sweep, feed change or geometry change outside the specified plain-disc replacements was used.

## FILES

`metrics.json` has unrounded peaks, SLL and comparisons; `gain_at_reference_directions.json` contains direct HFSS gain queries; `*_geometry.json`, `*_solve_status.json`, `*.prof`, and `*.ffd` preserve build and solver evidence. Private AEDT projects and STEP intermediates remain local and are not pushed.
