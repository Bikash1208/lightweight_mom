# Dual-beam (+-30 deg) 28 GHz nylon lens, design results (Phase 19 Part C) -- model-only

Regenerate: `scripts/research/phase19_design.py` (~15 min), then `scripts/research/phase19_local_periodicity.py` (~10 min).
HFSS comparison: `NOT_AVAILABLE_FROM_CURRENT_DATA`.

Setup: 150 mm octagon (across flats), 5 mm pitch, F = 180 mm (1.2 D), cos^q horn q = 14.38 (10 dB edge taper), x-pol,
25 mm total / 0.6 mm plate, through-hole pockets, nylon 2.45 / 0.0185, RCWA M = 8. 756 cells: 520 pillar, 236 hole.

## Synthesis (phase-only, exact pointing)
- Weighted Gerchberg-Saxton with a tight window and the ideal single-beam profile as in-window magnitude, mirror-symmetric (D2) by
  construction, window radius chosen so the beam peak lands on 30.00 deg. Beam peak = 0.394 of the ideal single-beam peak
  (binary baseline 0.403 at 30.35 deg); efficiency into the two beam windows (1.2 lambda/D)
  0.817 vs 0.798 for the binary 0/pi phase.
- **Finding:** with the beam directions held exact, continuous phase gives essentially the binary level (about -1 dB below an ideal equal
  split). Unconstrained GS reaches the ideal 0.5 but only by moving the beams to ~32 deg, which is not the requested design.

## Cell assignment (direct solve)
Every cell's octagon radius is solved by secant iteration of the full 2x2 Jones response at that cell's own (theta, phi) (library only as the
starting guess; family = pillar for delay <= 230 deg, hole above, overlap region). Phase residual: max 0.76 deg, rms 0.086 deg,
0 cells > 1 deg. Global phase offset psi0 = 228 deg.

## Predicted pattern (Huygens exit-aperture model, **not full-wave**)
Per-beam peak gain 25.6 dBi at 30.00 deg; level difference -0.00 dB (mirror symmetry);
SLL -9.7 dB; max cross-pol -19.3 dBi; aperture-averaged transmission 0.817;
75 % of horn power reaches the aperture and 60 % is radiated
(the difference is power in evanescent harmonics of the strongly phase-stepped aperture).

## Local-periodicity check (1x4 strips, exact supercell RCWA vs isolated-cell prediction)
- 432 of 1452 internal edges (29.8 %) are pillar|hole; median neighbour phase step across them
  122 deg (p90 177 deg).
- Order-0 error / cell amplitude: same-family pairs 0.01 (21 deg step), 0.07-0.09 (60 deg), 0.21 (109 deg), 0.47-0.54 (170 deg); pillar|hole pairs 0.14 (62 deg),
  0.26 (112 deg), 0.37-0.44 (150-170 deg), and two dense-hole|thin-pillar pairs 0.63 (26 deg) and 0.70 (83 deg). The error follows the neighbour phase step; the
  family boundary matters mainly when a dense hole cell sits next to a nearly empty pillar. Supercell solver truncation: order-0 changes 0.1 deg / 0.2 % from Mx = 16 to 32.
- Consequence: the HFSS result is expected to be worse than the numbers above. This design has not been optimised against neighbour coupling.
