# B-ALL-SHORT — full-board HFSS 2026 R1

## STATUS

COMPLETE. One new v3.1 board state was solved at 5.68 GHz. The requested three-part PASS gate **FAILS**: RHCP passes, while both LHCP limits fail.

## GATES

Exactly **250** slotted B cells received a central top-copper strip, 0.5 mm wide along the slot axis and spanning its full 2.91 mm width. Each 8.72 mm closed slot became two 4.11 × 2.91 mm closed slots; the strip rotated with the cell. The other eight B-labelled cells were already plain gate discs. All **16 gate discs** remained plain and ungrounded; zero vias. The per-site manifest matches BOTH for all coordinates, beam/gate labels, rotations, and non-target states. The added strip's 250 measured CAD body volumes agree with the two-slot area to within 0.002 mm³; STEP audit found zero inexact literals. AEDT read back 518 PEC metal objects and the same substrate, bias layer, absorber, feed, radiation boundary and 5.68 GHz setup settings. Automatic mesh settings were the same as BOTH; the adaptive mesh changed with the copper geometry.

The solve validated and converged in six adaptive passes (final max delta-S 0.00257), 46 cores, 82.9 GB peak HFSS memory/process, 1735 s wall time.

| PASS condition | HFSS result | Outcome |
|---|---:|---|
| Global LHCP peak ≤ −1.00 dBi | **−0.28 dBi** at (36°,178°) | FAIL |
| LHCP at (36°,180°) ≤ −1.00 dBi | **−0.41 dBi** | FAIL |
| RHCP at (32°,0°) ≥ +8.32 dBi | **+8.47 dBi** | PASS |

## DATA

| State | LHCP at (36°,180°) | Global LHCP peak | RHCP at (32°,0°) | SLL RHCP / LHCP | S11 | Accepted power |
|---|---:|---:|---:|---:|---:|---:|
| BOTH reference | +9.00 dBi | +9.00 dBi at (36°,180°) | +9.03 dBi | −8.43 / −9.52 dB | −16.76 dB | 97.89% |
| B-ALL-OFF-IDEAL | −3.48 dBi | −1.91 dBi at (44°,144°) | +8.60 dBi | −9.15 / −0.91 dB | −17.59 dB | 98.26% |
| **B-ALL-SHORT** | **−0.41 dBi** | **−0.28 dBi at (36°,178°)** | **+8.47 dBi** | **−9.26 / −1.94 dB** | **−17.73 dB** | **98.31%** |

SLL uses the previously fixed first-null main-lobe exclusion, rather than the 10° quick-look value. B-ALL-SHORT retains a strong LHCP residual near the original B direction; its global LHCP peak is 1.64 dB above the B-ALL-OFF-IDEAL residual.

## FLAGS

No vias, impedance sheets, feed changes, absorber changes, symmetry planes, frequency sweep or additional state solves were used. The two smaller closed slots were made by cutting them from each target disc, leaving a connected central 0.5 mm copper strip. The post-processed LHCP gain at (36°,180°) was computed from the stored 2° FFD grid, calibrated to HFSS's reported LHCP peak; this method exactly reproduces direct HFSS readbacks at the reference grid points in prior v3.1 states. The suppressed-hand SLL is relative to its residual global peak.

## FILES

`metrics.json` contains unrounded values and PASS booleans; `geometry.json` lists all site choices and 250 strips; `solve.prof`, `solve_status.json` and `B_ALL_SHORT.ffd` preserve solver and far-field evidence. The STEP and AEDT project remain local in private `hfss_work/full_board_gate_v31_short/` and are not pushed.
