# BOARD_V2_ELS read-back (UNSOLVED project, AEDT 2025.1)

STEP: `BOARD_V2.step`, 9271376 bytes, SHA256 `516c344bf0298129f171c38e88a9e4c12ea47d320c7aadd33b5f553782f5bca9`. 200 mm circular, 200 sites (140 static + 60 switched). Objects: METAL_PEC x 262 (pec), SUBSTRATE and BIAS_RT5880 (RT5880: eps_r 2.2, tan d 0.0009), AIR (vacuum), PORT, D1-D60 diode sheets.

| Design | Purpose | Diode sheet Rs / Xs (Ohm/sq), all 60 sheets | validate_simple |
|---|---|---|---|
| OFF_ELS10 | -10 V reverse bias; run first | 21.6 / -297.5 (60 of 60 read back) | 1 (1 = PASSED) |
| ON_ELS | 10 mA forward; run second | 2.8 / 14.6 (60 of 60 read back) | 1 (1 = PASSED) |
| OFF_ELS0 | 0 V; optional sensitivity | 41.8 / -284.4 (60 of 60 read back) | 1 (1 = PASSED) |

Identical in all three designs (read back per design in `readback.json`, after closing and reopening the saved file):
- Port: Lumped Port `Port1` on `PORT` (annulus r 0.635-2.05 mm at z = -3.175), integration line (0.635, 0, -3.175) -> (2.05, 0, -3.175), 50 Ohm, renormalised, no de-embedding.
- Boundaries: radiation boundary on AIR (['Rad__5BYA20']), one impedance boundary `BAR90_ELS` on the 60 sheets, lumped port.
- Setup `Setup1`: 5.8GHz, MaxDeltaS 0.02, MinimumPasses 6, MaximumPasses 20, BasisOrder 1, **DrivenSolverType = Direct Solver**; lambda refinement: DoLambdaRefine = True, SetLambdaTarget = False, Target = 0.3333, DoMaterialLambda = True, UseMaxTetIncrease = False.
- Sweep `Sweep`: Interpolating, 5.0-6.5 GHz, 0.05 GHz step, SaveFields False, SaveRadFields False (far fields come from the 5.8 GHz adaptive solution `Setup1 : LastAdaptive`).
- Far-field spheres: `Sphere1` theta 0deg-180deg, phi 0-360deg, step 2deg.
- **Mesh operations: none** (S3 mesh rule: automatic meshing only; mesh method Auto was never changed, the default initial mesh is used, and the API gives no read-back of the mesh-method switch). The four operations of the first build (`mesh_substrate` 3 mm, `mesh_cu_faces` 1 mm, `mesh_diodes` 0.15 mm, `mesh_pin_port` 0.3 mm) were deleted in all three designs; the old project is kept as `BOARD_V2_ELS_meshops_old.aedt`.
- Nothing has been solved: no analyze call was made on this project (the `existing_analysis_sweeps` list only names `Setup1 : Sweep`).

Not defined in the project: the RHCP output variable and the reports (S11, powers, RHCP cuts); AEDT refuses to create them before a solution exists. Their definitions are in `RUN_README.md`.
