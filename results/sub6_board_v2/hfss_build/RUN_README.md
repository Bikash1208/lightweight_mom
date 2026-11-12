# RUN_README - BOARD_V2_ELS (run by the user, one solve at a time)

Project `BOARD_V2_ELS.aedt` (this folder), geometry `BOARD_V2.step` (commit 5d0c0be), diode BAR90-02ELS (0.30 x 0.60 mm sheet, Zs = 2 Z). Three designs, identical except the diode sheet:

| Order | Design | Sheet (Ohm/sq) | Meaning |
|---|---|---|---|
| 1 | `OFF_ELS10` | 21.6 - j297.5 | OFF, -10 V reverse bias |
| 2 | `ON_ELS` | 2.8 + j14.6 | ON, 10 mA forward |
| 3 (optional) | `OFF_ELS0` | 41.8 - j284.4 | OFF at 0 V, sensitivity only |

## Run order and memory
- Run ONE solve at a time: OFF_ELS10, then ON_ELS, then (optional) OFF_ELS0. Two concurrent solves ran out of memory earlier.
- Expect about 110 GB of RAM per solve (a first adaptive pass with 2.56 million tetrahedra asked for 110 GB; later passes may need more). Start only with at least 130 GB free.
- Setup `Setup1`: 5.8 GHz adaptive, MaxDeltaS 0.02, at least 6 and at most 20 passes, first-order basis, direct solver; then the interpolating sweep 5.0-6.5 GHz. Use at least 8 cores; save the project after each solve.

## Checks before reading ANY gain
1. Convergence: the final Max Mag. Delta S is at or below 0.02 and at least 6 passes were run. If 20 passes are reached without it, report the passes, the final delta S and the tetrahedra, and do not use the far field.
2. Accepted power at 5.8 GHz is at least 50 % (accepted = 1 - |S11|^2 for the single port; the Antenna Parameters report `AcceptedPower` must agree). If it is below 50 %, stop and report.

## Reports and expressions to create after the solve (they need a solution to exist)
- S11: `dB(S(Port1,Port1))` vs frequency, solution `Setup1 : Sweep`.
- Powers at 5.8 GHz: Antenna Parameters on `Sphere1`, solution `Setup1 : LastAdaptive`: `IncidentPower`, `AcceptedPower`, `RadiatedPower`.
- RHCP, IEEE definition, as an output variable `RHCP` on `Sphere1`: `(rETheta+i*rEPhi)/sqrt(2)`; gain cut: `dB(abs((rETheta+i*rEPhi)/sqrt(2)))` vs Theta, solution `Setup1 : LastAdaptive`, Freq = 5.8 GHz, Phi = 0 deg and Phi = 180 deg (two curves).
- Far-field export: right-click `Sphere1` > Export Far Field (or `ExportFieldsToFile`) as `.ffd`, 5.8 GHz, one file per design.

## What to report, per design
- RHCP level at theta = 24 deg and 26 deg for phi = 0 deg and 180 deg (the sphere has a 2 deg grid, so ask for 24 and 26 explicitly).
- OFF_ELS10 only: OFF balance of beam A against beam B in dB (A at phi 180, B at phi 0 in this rot180 layout, so compare the two phi cuts).
- ON_ELS: suppression of beam B in dB relative to OFF_ELS10.
- Total radiated power and the power lost below the board (RadiatedPower; the part of the sphere with theta > 90 deg).
- The S11 over 5.0-6.5 GHz and the `.ffd` exports.
- For each solve: passes, final delta S, tetrahedra, peak RAM, solve time.

## Notes
- The board is the 200 mm circular board with 60 diodes, no FR4 (RT5880 3.175 mm / ground / RT5880 0.254 mm bias layer), vias only (no bias network). The static A cells are FR4-calibrated (see `manifest.json` in the exporter output).
- Earlier OFF result: the old-sheet OFF solve (22.23 - j335.4) was stopped before it finished; the old-sheet ON solve failed with out-of-memory. Neither produced results.
