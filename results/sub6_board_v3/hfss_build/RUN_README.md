# RUN_README - BOARD_V3_ELS (run by the user, one solve at a time)

Project `BOARD_V3_ELS.aedt` (this folder), geometry `BOARD_V3.step` (SHA256 `c68e7d1ddb000fef87b7e60e44c67dc0be4fc03a46a8d82c0a2e5940e39ba6e5`, exporter commit 2559529), diode BAR90-02ELS (0.30 x 0.60 mm sheet, Zs = 2 Z), **circular board R = 180 mm (D = 360 mm)**, 696 sites (636 static + 60 switched). Three designs, identical except the diode sheet:

| Order | Design | Sheet (Ohm/sq) | Meaning |
|---|---|---|---|
| 1 | `OFF_ELS10` | 21.6 - j297.5 | OFF, -10 V reverse bias |
| 2 | `ON_ELS` | 2.8 + j14.6 | ON, 10 mA forward |
| 3 (optional) | `OFF_ELS0` | 41.8 - j284.4 | OFF at 0 V, sensitivity only |

## Run order, solver and memory
- Run ONE solve at a time: OFF_ELS10, then ON_ELS, then (optional) OFF_ELS0.
- Setup `Setup1`: 5.8 GHz adaptive, MaxDeltaS 0.02, at least 6 and at most 20 passes, first-order basis; then the interpolating sweep 5.0-6.5 GHz. Use at least 8 cores; save the project after each solve.
- **Matrix solver: ITERATIVE** (`DrivenSolverType` = Iterative Solver, set through the API and read back after closing and reopening the saved project, in all three designs). To use the direct solver instead: Setup1 > Advanced > Solver Type > Direct Solver.
- Memory: the small board (R = 100 mm) needed about 110 GB with the direct solver for 2.56 million tetrahedra. This board is about 3.5 times larger; see `readback.md` for the initial-mesh tetrahedra count (mesh-only, nothing solved) and the estimate derived from it. Start only with at least 130 GB free.

## Meshing: automatic only (S3 mesh rule)
- There are NO custom mesh operations in any design (read back after reopening the saved file: none). HFSS automatic meshing is used: mesh method Auto, the default initial mesh, lambda refinement ON (DoLambdaRefine = True), and the adaptive passes refine the mesh.
- **Please report the tetrahedra count at every adaptive pass** (and the initial-mesh count), together with the Max Mag. Delta S per pass, so we can see whether any mesh seeding is actually needed. A mesh operation is added only for a specific failure (no convergence, a failed initial mesh, or a clearly unresolved small feature such as a diode sheet or the port), and then the operation, objects, size, reason and before/after passes, delta S and tetrahedra are reported.
- The first build of this project had four operations (substrate 3 mm, copper faces 1 mm, diode sheets 0.15 mm, pin and port 0.3 mm); they were deleted and the old file is kept as `BOARD_V3_ELS_meshops_old.aedt`.

## Checks before reading ANY gain
1. Convergence: the final Max Mag. Delta S is at or below 0.02 and at least 6 passes were run. If 20 passes are reached without it, report the passes, the final delta S and the tetrahedra, and do not use the far field.
2. Accepted power at 5.8 GHz is at least 50 % (accepted = 1 - |S11|^2 for the single port; the Antenna Parameters report `AcceptedPower` must agree). If it is below 50 %, stop and report.

## Reports and expressions to create after the solve (they need a solution to exist)
- S11: `dB(S(Port1,Port1))` vs frequency, solution `Setup1 : Sweep`.
- Powers at 5.8 GHz: Antenna Parameters, solution `Setup1 : LastAdaptive`: on `Sphere1` read `IncidentPower`, `AcceptedPower`, `RadiatedPower` (total radiated power, theta 0-180 deg). On sphere `Lower` (theta 90-180 deg, phi 0-360 deg) read `RadiatedPower`: **that is the power lost below the board**; the power into the upper hemisphere is `RadiatedPower(Sphere1)` minus `RadiatedPower(Lower)`.
- RHCP, IEEE definition, relative: `dB(abs((rETheta+i*rEPhi)/sqrt(2)))` vs Theta on `Sphere1`, Freq = 5.8 GHz, Phi = 0 deg and Phi = 180 deg.
- RHCP, **absolute gain** (use this one for dBi): `10*log10(2*pi*mag((rETheta+i*rEPhi)/sqrt(2))^2/(376.73*AcceptedPower))`, solution `Setup1 : LastAdaptive`, on `Sphere1`, Freq = 5.8 GHz, vs Theta for Phi = 0 deg and 180 deg.
- Far-field export: right-click `Sphere1` > Export Far Field as `.ffd`, 5.8 GHz, one file per design.

## What to report, per design
- RHCP level at theta = 24 deg and 26 deg for phi = 0 deg and 180 deg (the sphere has a 2 deg grid, so ask for 24 and 26 explicitly).
- OFF_ELS10 only: OFF balance of beam A against beam B in dB (A at phi 180, B at phi 0 in this rot180 layout, so compare the two phi cuts).
- ON_ELS: suppression of beam B in dB relative to OFF_ELS10.
- Total radiated power and the power lost below the board.
- The S11 over 5.0-6.5 GHz and the `.ffd` exports.
- For each solve: passes, final delta S, tetrahedra per pass, peak RAM, solve time.

## Notes
- Layout (jobs1_iter1_study.py, rot180, k = 1, N = 60, placement `inner`, cap 0.5 Np/m): the study's own read-out for this board is A/B OFF -9.9 / -17.3 dB (**A-B balance +7.4 dB, not within 3 dB**), B suppression in ON 8.8 dB, converted power 0.203 (OFF) / 0.118 (ON). So an A/B balance well above 3 dB and an ON suppression below 10 dB are what the model predicts; model only, not HFSS-validated.
- No FR4 (RT5880 3.175 mm / ground / RT5880 0.254 mm bias layer), vias only (no bias network); static A cells are FR4-calibrated; 12 x 12 mm SMA window through the bias layer; feed block from `results/sub6_feed/fix1/`.
- Nothing has been solved in this project.
