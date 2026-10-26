# Phase 20 R4 item 1 — symmetry reduction: tested, NOT usable

**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET.** Cell `ris_cp2_octagon_madp_rt5880_15` (R3 finalist A, MADP, L_s = 0), mesh 60 (the R3 search mesh), 9 / 10 / 11 GHz,
theta = 0 / 30 / 45 deg, phi = 0, 15, ..., 345 deg, all four one-hot states. Scripts `scripts/research/phase20_r4_symmetry.py`, `phase20_r4_mirror.py`; data `symmetry.json`, `mirror.json`.

## Claim tested
State k at incidence azimuth phi equals state 0 at phi - 45 k deg in the incidence-plane TE/TM frame (the whole structure and the incident direction turned by 45 k deg; the
TE/TM frame turns with the plane of incidence, so the 2x2 matrices should be identical).

## Why it can only be approximate (and was expected to need a check)
* The three OFF stubs are not turned onto themselves: turning state 0 by 45 k deg moves the stub set {0, 45, 90, 135} to {45 k, ..., 45 k + 135}, so a stub appears where the rotated cell has none
  (an OFF diode is not an absent stub: 0.025 pF plus the pads and line are metal).
* A 45 deg turn maps the square lattice onto a rotated lattice (only k = 2, 90 deg, is a lattice symmetry).

## Result (216 points x 3 states): the relation does NOT hold to mesh-staircase tolerance
| quantity (state k at phi vs state 0 at phi - 45 k) | median | 95th percentile | max |
|---|---|---|---|
| r_RR phase error | 12.8 deg | 20.9 deg | 63.2 deg |
| r_RR level error | 0.29 dB | 2.2 dB | 6.6 dB |
| max entry of the 2x2 difference \|S^k - S^0\| (S ~ 1) | 0.18 | 0.37 | 1.36 |

By state (r_RR phase error, median / max): k = 1: 4.7-5.8 / 5.5-53 deg (theta 0 -> 45); k = 2: 15.7 / 16 deg at theta 0 and 15.1 / 23 deg at 45 deg; k = 3: 14.2 / 15 deg at theta 0, 13.3 / 63 deg at 45 deg.
Even k = 2 (a 90 deg turn, a true lattice symmetry) is 15-17 deg off, almost independent of phi: that is the OFF-stub asymmetry, the same effect as the R3 ladder steps of 94 / 108 / 80 deg instead of 90.

## Mirror symmetry (state k about its own stub axis, phi -> 90 k - phi, cross terms change sign)
Not usable either: the OFF stubs have no mirror partners ({0, 45, 90, 135} is symmetric about no member); the 2x2 difference is \|dS\| median 0.18, 95th percentile 0.24, max 0.26 for all four states.
Large phase errors on individual entries (up to 160 deg) occur where an entry is small, but the entry-wise difference of 0.2 of the incident amplitude is the meaningful measure.

## Decision
No reduction: **all four states are solved over the full phi range [0, 360)** (no rotation, no mirror shortcut). The only savings that are exact: the four states share one impedance matrix per
(f, theta, phi) (`compute_s_parameters_oblique_over_loads`, tested identical to separate solves), and at theta = 0 phi only relabels the TE/TM frame (one solve per (f, state), the rest follows by
the exact frame rotation `cartesian_reflection` / `jones_to_te_tm`).
