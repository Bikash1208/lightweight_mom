# Step 1 — coherent feed reflection, A-2q, closed gate at the ring (MODEL_ONLY)

`scripts/research/sidequest_gate_reflection.py` -> `gate_reflection.json` (25 s). lambda_g = 2 pi / (1.065 k0) = 26.7 mm; ring front radius
r0 = 13 mm; sweep r0 +- lambda_g/4 (+-6.7 mm = 2 pi of round-trip phase); reflection phase of the gate is absorbed in the radius.

Model: cylindrical-harmonic multiple scattering of the surface wave, 121 harmonics. Matched point probe (unitarity with a matched
port forces s_p = 0: the returning n = 0 wave is absorbed, `LAUNCH` of it into the port; n != 0 harmonics pass the origin freely and
re-emerge towards the ring, i.e. "the reflected wave that passes the probe into the open quadrant" is included). Blocked half of the
ring (psi in [90, 270] deg) reflects with Gamma, T = 0; open half passes, IL = 0. Exterior field = half-aperture re-expansion
(diffraction at the sector edges included; it moves ~0.08 of the source power across the sector edge, and the Kirchhoff-style
re-expansion is about 6 % off in power vs the interface flux, so power numbers carry about +-0.06 of the source power).
Closure check: flux through the ring = 1 - |B_0|^2 to 1e-3 (the harmonic solve is self-consistent).

S11 below is the **reflection added by the gate to an ideally matched feed** (bare-feed Gamma_0 = 0). With the HFSS bare baseline
|Gamma_0| = 0.178 (-15 dB) of unknown phase the worst-case coherent sum is |Gamma_0| + |S11_gate|.

| Gate block state | best radius shift | S11 gate-only | worst case with -15 dB bare feed | open-quadrant power vs no reflection | gain vs open-circuit-free baseline |
|---|---|---|---|---|---|
| closed, \|Gamma\| = 1 | +0.8 ... +3.3 mm (r = 13.8-16.3 mm) | -13.4 dB (-12.9 ... -13.4) | -8.2 dB | 0.75 vs 0.42 (x1.8) | **+2.2 dB** (20.7 vs 18.4 dBi), SLL -13.5 to -14.8 dB |
| closed, same, worst radius (-3.3 mm) | | -5.6 dB | -3 dB | 0.08 (x0.2) | -7.5 dB |
| lossy, \|Gamma\| = 0.5 | +1.7 ... +3.3 mm | -17.9 dB | -10.3 dB | 0.60 (x1.4) | +1.5 dB (power only) |
| lossy, \|Gamma\| = 0.3 | +2.5 mm | -21.7 dB | -12.6 dB | 0.54 (x1.3) | +1.0 dB (pattern computed) |
| lossy, \|Gamma\| = 0.1 | any | -30 dB | -15 dB | 0.47 (x1.1) | +0.4 dB |

(Gamma = 1, 0 mm: -11.9 dB, 0.73, +2.3 dB; the in-phase window is wide: for dr from -0.8 to +6.7 mm S11 <= -10.5 dB and p_open >= 0.67.)

## Findings
1. The earlier bound (0.56 amplitude) was far too pessimistic: a half-ring reflector returns only its n = 0 harmonic to the probe
   (|B_0| ~ 0.5 x 0.62 -> -12 dB), the rest passes through the centre into the open half.
2. **A closed (|Gamma| = 1) gate at r_front ~ 14-16 mm gives S11 = -13 dB (gate only) AND +2.2 dB gain** in the open quadrant (reflection
   adds in phase). Radii near r_front = 9.7-11.2 mm (dr -3.3 to -2.5) must be avoided (destructive, -6 dB S11, -7.5 dB gain).
   Gate position tolerance for the good window: about +-3 mm = +-0.12 lambda_g.
3. With the -15 dB bare feed of unknown phase, a closed gate cannot guarantee S11 <= -10 dB (worst case -8.2 dB); a lossy block with
   |Gamma| <= 0.5 does (worst case -10.3 dB) and still keeps +1.5 dB. Design target for the cell: **block-state |Gamma| in 0.3-0.5 with
   T <= -20 dB**, or |Gamma| ~ 1 if the bare feed can be re-matched to <= -20 dB.
4. Not modelled: gate IL in the pass state, finite gate thickness (3 rows = 19 mm: the ring front radius is used), TM/TE splitting,
   angle dependence of the real gate reflection.
