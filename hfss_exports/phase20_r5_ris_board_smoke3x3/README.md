# ris_smoke3x3 -- 3 x 3 CP-2 RIS board for HFSS (Phase 20 R5-6)

**State: INFRASTRUCTURE / MODEL_ONLY / NOT_HFSS_VALIDATED.** Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION, AEDT_CALLS_VERIFIED_3X3_SMOKE_2026-10-24. **SMOKE TEST ONLY -- NO PREDICTIONS, NOTHING TO COMPARE WITH.** This 3 x 3 board (36 diodes) exists to check that `hfss_build_board.py` builds the model in AEDT and that the two calls that have never
been run there work: `AssignLumpedRLC` (36 boundaries; ON sheets R only, OFF sheets C only) and the incident-wave assignment. The coding is the plain 2-bit steering coding for 30 deg in the phi = 0 plane
(any valid coding would do); it carries no design meaning at 3 x 3 cells. What to look at after the run: (1) the build finishes without the message "Lumped RLC assignment failed" or "Incident wave not assigned";
(2) HFSS > Boundaries lists `Load_1_RLC` ... `Load_36_RLC`, nine with a resistance of 5.2 ohm and 27 with a 0.025 pF capacitance, each with its current line along the stub;
(3) the excitation list shows the incident wave with the intended direction (the specular direction is theta_i = 0 here, i.e. normal incidence). Report what failed (exact message) so the script can be fixed.

**Result (user, 2026-10-24): verified (3x3 smoke build, 2026-10-24)** -- the build completed with no errors; `Load_1`..`Load_36_RLC` are present with R = 5.2 ohm (ON) and C = 0.025 pF (OFF); the incident wave was assigned. Still to confirm by the user: the current-line direction of the boundaries. The circular-polarisation incidence (two plane waves) remains a manual step.

Board: 45 x 45 mm, 9 cells (15 mm pitch), 36 diodes (MADP-000907-14020, PROVISIONAL_DATASHEET: ON R_s = 5.2 ohm, OFF C_T = 0.025 pF, L_s = 0 assumed, flip-chip pads not modelled),
RT/duroid 5880 eps_r 2.2, tan d 0.0009, h = 1.575 mm, PEC ground, via omitted. State counts (states 1-4): 3 / 3 / 3 / 0. Design frequency 10.5 GHz.

## Files
`hfss_build_board.py` (AEDT script: Tools > Run Script), `layout.csv` (one row per cell: position, state 1-4), `diodes.csv` (one row per `Load_<k>`: position, size, current-line angle, ON/OFF, element), `board.step` (STEP with `Substrate`, `Metal`, `Ground`, `Load_<k>`; gitignored, regenerate with the exporter).

## Conventions
* x along columns, y along rows, board centred on the origin; cell (row j, column i) has state s = layout.csv `state_1to4`; its stub `s-1` (angle 45 (s-1) deg) has the diode ON, the other three OFF.
* `Load_<k>`: k = 4 (j nx + i) + stub + 1. Each is a sheet bridging the 0.75 mm diode gap; the current line runs ALONG the stub (angle in `diodes.csv`). **One single-type Lumped RLC boundary per sheet**: ON = R only, OFF = C only (L_s = 0).
* Metal is zero-thickness (sheets at z = h, Perfect E), like the MoM; the `board.step` metal is the 35 um slab (state which you used).
* Incident wave: propagation-relative **RHCP** at the tangential wave vector of the specular direction (theta_i, phi_i) = (0, 0) deg. In HFSS: two plane waves of equal amplitude, theta- and phi-polarised, 90 deg apart in phase
  with the sign that gives right-handed rotation about the propagation direction (check with the bare grounded slab: the reflected wave must be left-handed there, `docs/theory/CONVENTIONS.md`). The script assigns one linear wave as a placeholder.
* Radiation box padded by 15 mm; far-field sphere theta 0-90 deg, phi 0-360 deg. Report the scattered far field: co-handed (RHCP) and opposite-handed (LHCP) gain, normalised to the power incident on the board area.

## Honest limits
Model uncertainty ~10 deg on the cell phases (staircased octagon in the MoM), the local-response assumption of the forward model (checked at array level by the R5-4 server runs), no via, no bias network, lumped diode without pad geometry, finite 150 mm board edge effects in HFSS that the aperture model does not have.
The two new AEDT calls (`AssignLumpedRLC`, `AssignPlaneWave`) are verified in AEDT (3x3 smoke build, 2026-10-24); the current-line direction is for the user to confirm, and the circular-polarisation incidence (two waves) stays a manual step.
