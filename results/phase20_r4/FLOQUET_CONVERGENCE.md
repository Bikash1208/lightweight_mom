# Phase 20 R4 — Floquet-cutoff finding (affects how every number in the RIS work must be read)

**Labels: MODEL_ONLY.** Found while building the CP-2 library (R4 item 3) and cross-checked by an independent grid audit (read-only agents + my own checks).

## What happened
The solver accepts any Floquet cutoff F with 2F + 1 >= (cells per period); R2/R3 and the first R4 grid used the MINIMUM (F = n/2 + 1). At isolated (f, theta, phi) points that cutoff leaves
**narrow numerical lines**: a spurious resonance a few tens of MHz wide, with skirts, whose position moves with F. Example (mesh 60, theta 10, phi 45, state 3 level): F = 60 shows a dip to -8 dB at
10.610 GHz and is 22 degrees off 50 MHz away; F = 90 is smooth but still 11 degrees off F = 120 at 10.56 GHz; **F = 120 and F = 150 agree to < 1 degree** over 10.56-10.62 GHz. Other examples: mesh 60 /
F 31 at 8.2 and 11.35 GHz (up to 10 dB off, every angle), mesh 72 / F 37 at the same frequencies. F = 2n (n = cells per period) is the converged choice; 0.67 n and 0.83 n are not safe.

## What it does NOT change
With F = 2n the mesh convergence is **regular** (close to second order), which the minimal-F data hid ("observed order unusable" in R3 was partly this noise). 10 GHz, finalist A, phi 0, r_RR phase
(states 1-4) and opposite-handed level, `converged_study.json`:

| theta | mesh / F | r_RR phase (deg) | opposite-handed (dB) |
|---|---|---|---|
| 0 | 60 / 31 (R3 convention) | -113.1 -17.6 82.9 171.2 | -23.3 -16.4 -18.8 -20.9 |
| 0 | 60 / 120 | -112.2 -15.4 84.1 173.2 | -23.3 -15.8 -18.7 -20.7 |
| 0 | 90 / 180 | -106.8 -13.4 87.5 175.1 | -16.5 -15.0 -13.5 -19.1 |
| 0 | 120 / 240 | -104.7 -13.5 90.6 175.8 | -18.5 -17.7 -14.4 -20.8 |
| 30 | 60 / 120 | -110.0 -14.6 79.9 169.7 | -14.5 -16.9 -17.3 -13.0 |
| 30 | 120 / 240 | -101.5 -12.8 86.4 172.4 | -15.9 -15.4 -18.0 -12.4 |

Mesh 60 -> 120 correction at converged F: +7.5, +1.9, +6.5, +2.6 deg (states 1-4) at theta 0 — the same 2-7.5 deg the minimal-F pilot gave. The R3 headline numbers (ladder deviation 12.1 deg,
step errors 15.7-18.0 deg, amplitudes -0.9 to -1.1 dB, opposite-handed -13.6 / -14.0 dB) therefore move by a few degrees at most at the points checked (10 GHz); they were not recomputed. The R3 verdicts
stand within their ~10 deg uncertainty; the **minimal-F HFSS package** `hfss_exports/phase20_r3_cp_cell/mom_reference.csv` (0.5 GHz steps, F minimal) may contain an affected frequency and should be
regenerated at F = 2n before it is used for a comparison (listed in STATUS).

## Cost (all four states per point, this machine, full BLAS)
mesh 60: F = 31: 1.4-6 s, F = 120: 12 s; mesh 90: F = 46: 14 s, F = 180: 121 s; mesh 120: F = 61: 70 s, F = 240: 616 s. This is why the approved grid (8897 points, 0.1 GHz x 5 deg x 15 deg, 41 x 217 points)
was not built at F = 2n (about 18-30 h plus about 24 h of corrected points, past the brief's 24 h stop rule): the library uses the reduced grid f every 0.25 GHz, theta every 7.5 deg, phi every 30 deg (1241 points),
correction points at f every 0.5 GHz, theta {0, 30, 45}, phi {0, 90} (meshes 60 / 90 / 120, 45 points).
