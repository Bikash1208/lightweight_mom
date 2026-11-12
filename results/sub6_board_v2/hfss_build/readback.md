# BOARD_V2_ELS read-back (UNSOLVED project, AEDT 2025.1)

STEP: `BOARD_V2.step` = commit 5d0c0be, 9271376 bytes, SHA256 `516c344bf0298129f171c38e88a9e4c12ea47d320c7aadd33b5f553782f5bca9`. 326 objects: METAL_PEC x 262 (pec), SUBSTRATE and BIAS_RT5880 (RT5880: eps_r 2.2, tan d 0.0009), AIR (vacuum), PORT, D1-D60 diode sheets.

| Design | Purpose | Diode sheet Rs / Xs (Ohm/sq), all 60 sheets | validate_simple |
|---|---|---|---|
| OFF_ELS10 | -10 V reverse bias; run first | 21.6 / -297.5 (60 of 60 read back) | 1 (1 = PASSED) |
| ON_ELS | 10 mA forward; run second | 2.8 / 14.6 (60 of 60 read back) | 1 (1 = PASSED) |
| OFF_ELS0 | 0 V; optional sensitivity | 41.8 / -284.4 (60 of 60 read back) | 1 (1 = PASSED) |

Identical in all three designs (read back per design in `readback.json`):
- Port: Lumped Port `Port1` on `PORT` (annulus r 0.635-2.05 mm at z = -3.175), integration line (0.635, 0, -3.175) -> (2.05, 0, -3.175), 50 Ohm, renormalised, no de-embedding.
- Boundaries: radiation boundary on AIR (['Rad__5BYA20']), one impedance boundary `BAR90_ELS` on the 60 sheets, lumped port.
- Setup `Setup1`: 5.8GHz, MaxDeltaS 0.02, MinimumPasses 6, MaximumPasses 20, BasisOrder 1, direct solver (default; the `UseIterativeSolver` property is not present in this setup, so it reads `None`).
- Sweep `Sweep`: Interpolating, 5.0-6.5 GHz, 0.05 GHz step, SaveFields False, SaveRadFields False (far fields come from the 5.8 GHz adaptive solution `Setup1 : LastAdaptive`).
- Far field: infinite sphere `Sphere1`, theta 0-180 deg and phi 0-360 deg, 2 deg grid, Theta-Phi, linear polarisation (theta to 180 deg so the power below the board is captured).
- Mesh operations: ['mesh_substrate', 'mesh_diodes', 'mesh_cu_faces', 'mesh_pin_port'] (substrate 3 mm, copper faces 1 mm, diode sheets 0.15 mm, pin and port 0.3 mm).
- Nothing has been solved: no analyze call was made; the project was saved unsolved and closed cleanly (no lock file).

Not defined in the project: the RHCP output variable and the reports (S11, powers, RHCP cuts). AEDT refused to create them before a solution exists (the first attempt failed with "Invalid commands" and broke Save, so they were removed). Their definitions are in `RUN_README.md`.

Not read back: the interpolating sweep tolerance (AEDT default).
