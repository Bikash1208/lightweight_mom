# Phase 20 R4 item 2 — grid and cost estimate (STOP POINT: the full grid is not built yet)

**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET.** Cell `ris_cp2_octagon_madp_rt5880_15` (R3 finalist A, MADP, L_s = 0 nominal). Machine: 12 cores, 36 GiB. Measurements:
`scripts/research/phase20_r4_pilot.py`, `phase20_r4_symmetry.py` (`pilot.json`, `symmetry.json`, `SYMMETRY.md`). **Caution: the machine's wall clock jumps (R4 item 0), all times below are CPU/monotonic seconds.**

## Findings that shape the plan
1. **No symmetry reduction** (`SYMMETRY.md`): all four states, full phi in [0, 360). Exact savings only: the four states share one impedance matrix per (f, theta, phi); at theta = 0 phi is a pure frame rotation (one solve per f).
2. **CP-2 is light on memory** (small metal fraction): unknowns 2352 / 5354 / 9456 at mesh 60 / 90 / 120, peak 0.9 / 3.4 / 9.0 GB — mesh 120 (the R3 finest mesh) fits easily, unlike the R2 cell (15.6 GB). (A mesh 160 point would be ~17k unknowns, ~25-30 GB: possible singly, not planned.)
3. **Cost per (f, theta, phi) point, all four states, measured on this machine** (single process, all cores): mesh 60: **1.4 s**; mesh 90: **14.2 s**; mesh 120: **70.2 s** (Floquet cutoff n/2 + 1).
4. **Pilot of the mesh-correction map** (10 GHz, theta 30 at phi = 0, 45, ..., 315 and theta 45 at phi = 0, 90, 180, 270; r_RR, all states; correction = mesh 120 minus mesh 60):

| state | correction of the r_RR phase, theta 30 (mean, range over phi) | theta 45 | correction of \|r_RR\| (dB), theta 30 |
|---|---|---|---|
| 1 | +7.5 deg (6.9 ... 8.2) | +7.4 (6.1 ... 8.7) | -0.17 |
| 2 | +2.5 (2.3 ... 2.7) | +2.6 (2.6) | -0.03 |
| 3 | +6.8 (6.3 ... 7.1) | +6.6 (5.8 ... 7.5) | -0.25 |
| 4 | +3.1 (2.9 ... 3.3) | +3.0 (2.8 ... 3.2) | -0.05 |

   The correction depends on the **state** (2.5-7.5 deg) but hardly on **phi** (spread <= 1.3 deg) or **theta** (30 vs 45 agree to <= 0.4 deg mean): a smooth map over (f, theta) per state, sampled sparsely, is justified at 10 GHz.
   The mesh 90 -> 120 step is 1.4-3.1 deg (not yet converged at 90; R3's observed order was unusable, first-order fallback applies). Whether the map stays smooth in **frequency** (the discretisation shifts the effective resonance, so near resonance the correction follows the phase slope) is NOT yet known: the correction set below samples f every 0.25 GHz to measure it.

## Proposed grid (the brief's, with the reduction that is allowed)
f = 8.0-12.0 GHz step 0.1 (41); theta = 0, 5, ..., 45 deg (θ > 30 labelled `EXTRAPOLATED_BEYOND_R3_TARGETS`); phi step 15 deg over [0, 360) (24; theta = 0: one point); states 1-4; TE and TM.
Points: 41 x (1 + 9 x 24) = **8897 (f, theta, phi) points** x 4 states.

| step | points | time per point | wall time (estimate) | peak memory |
|---|---|---|---|---|
| A. grid at mesh 60 | 8897 | 1.4 s | **3.5 h sequential** (less with parallel single-thread workers; 11 workers should give ~1-1.5 h, to be measured) | < 1 GB per worker |
| B. correction set: mesh 90 and 120 at f every 0.25 GHz (17) x theta {0, 15, 30, 45} x phi {0, 90} (theta 0: one) = 17 x 7 = 119 points | 119 | 14.2 + 70.2 s | **2.8 h** | 9 GB (one at a time) |
| C. held-out validation (item 4): 50 random off-grid points at mesh 60 + 20 of them at mesh 120 | 50 + 20 | 1.4 s / 84 s | 0.5 h | 9 GB |
| D. tolerance corners (item 5): f {9.5, 10, 10.5} x theta {0, 30} x phi {0, 45, 90} x 4 corners at mesh 60 | ~60 | 1.4 s | < 5 min | < 1 GB |
| **Total** | | | **about 7 h sequential, 4-5 h with parallel workers (A)** | 9 GB |

**The full grid is within the 24 h limit; no reduced grid is needed.** Alternatives considered and rejected: mesh 90 for the whole grid = 8897 x 14.2 s = 35 h (exceeds the limit; mesh 90 on f step 0.25 would be 14.5 h and still not converged: 1.4-3.1 deg from mesh 120); mesh 120 for the whole grid = 173 h.

## Risks to flag before building
* The correction map may not be smooth in f near the resonance; if the f-sampling shows otherwise the correction becomes a per-frequency lookup (still B's 17 x 7 points) and the uncertainty in the library says so.
* Model uncertainty (staircased octagon, L_s = 0 assumed, MADP datasheet at 1 GHz / 1 MHz) is NOT reduced by any of this; the library is the MoM table of the registered cell, ~10 deg phase uncertainty remains until the user's HFSS run (R3's package) calibrates it.
* The library is stored in the incidence-plane (TE, TM) frame at azimuth phi; R5 must use the same convention (`docs/theory/CONVENTIONS.md`).
