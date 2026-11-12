# BOARD_V3_ELS read-back (UNSOLVED project, AEDT 2025.1)

STEP: `BOARD_V3.step`, 19946802 bytes, SHA256 `c68e7d1ddb000fef87b7e60e44c67dc0be4fc03a46a8d82c0a2e5940e39ba6e5`. R = 180 mm (D = 360 mm) circular, 696 sites (636 static + 60 switched). Objects: METAL_PEC x 758 (pec), SUBSTRATE and BIAS_RT5880 (RT5880: eps_r 2.2, tan d 0.0009), AIR (vacuum), PORT, D1-D60 diode sheets.

| Design | Purpose | Diode sheet Rs / Xs (Ohm/sq), all 60 sheets | validate_simple |
|---|---|---|---|
| OFF_ELS10 | -10 V reverse bias; run first | 21.6 / -297.5 (60 of 60 read back) | 1 (1 = PASSED) |
| ON_ELS | 10 mA forward; run second | 2.8 / 14.6 (60 of 60 read back) | 1 (1 = PASSED) |
| OFF_ELS0 | 0 V; optional sensitivity | 41.8 / -284.4 (60 of 60 read back) | 1 (1 = PASSED) |

Identical in all three designs (read back per design in `readback.json`, after closing and reopening the saved file):
- Port: Lumped Port `Port1` on `PORT` (annulus r 0.635-2.05 mm at z = -3.175), integration line (0.635, 0, -3.175) -> (2.05, 0, -3.175), 50 Ohm, renormalised, no de-embedding.
- Boundaries: radiation boundary on AIR (['Rad__VUL9BK']), one impedance boundary `BAR90_ELS` on the 60 sheets, lumped port.
- Setup `Setup1`: 5.8GHz, MaxDeltaS 0.02, MinimumPasses 6, MaximumPasses 20, BasisOrder 1, **DrivenSolverType = Iterative Solver**; lambda refinement: DoLambdaRefine = True, SetLambdaTarget = False, Target = 0.3333, DoMaterialLambda = True, UseMaxTetIncrease = False.
- Sweep `Sweep`: Interpolating, 5.0-6.5 GHz, 0.05 GHz step, SaveFields False, SaveRadFields False (far fields come from the 5.8 GHz adaptive solution `Setup1 : LastAdaptive`).
- Far-field spheres: `Sphere1` theta 0deg-180deg, phi 0-360deg, step 2deg; `Lower` theta 90deg-180deg, phi 0-360deg, step 2deg.
- **Mesh operations: none** (S3 mesh rule: automatic meshing only; mesh method Auto was never changed, the default initial mesh is used, and the API gives no read-back of the mesh-method switch). The four operations of the first build (`mesh_substrate` 3 mm, `mesh_cu_faces` 1 mm, `mesh_diodes` 0.15 mm, `mesh_pin_port` 0.3 mm) were deleted in all three designs; the old project is kept as `BOARD_V3_ELS_meshops_old.aedt`.
- Nothing has been solved: no analyze call was made on this project (the `existing_analysis_sweeps` list only names `Setup1 : Sweep`).

Not defined in the project: the RHCP output variable and the reports (S11, powers, RHCP cuts); AEDT refuses to create them before a solution exists. Their definitions are in `RUN_README.md`.

## Initial mesh only (automatic meshing, no mesh operations, nothing solved)
```
{"generate_mesh_ok": true, "mesh_seconds": 690.4, "tetrahedra_found": [], "meshstats_error": null, "mesh_ops_present": []}
```

## Memory estimate (mesh-only, `Analyze` never run)
- Initial mesh (automatic meshing, no mesh operations), design OFF_ELS10: **739,263 tetrahedra** after a 690 s mesh-only run (AIR 235,297; SUBSTRATE 205,139; BIAS_RT5880 60,621; ground 91,600; the other copper bodies about 130 each).
- Calibration: on BOARD_V2_ELS (first build) 2,562,312 tetrahedra at adaptive pass 1 made the direct solver ask for about 110 GB, i.e. about 43 KB per tetrahedron.
- So the direct solver would need about **32 GB** at 739k tetrahedra, growing with the adaptive passes (about 110 GB once the mesh reaches ~2.6 million tetrahedra); the pass-1 count can differ from the initial-mesh count because port and lambda refinement follow.
- `DrivenSolverType` = Iterative Solver is set in all three designs (read back after closing and reopening the file); it normally needs much less memory than the direct solver (not measured here).
