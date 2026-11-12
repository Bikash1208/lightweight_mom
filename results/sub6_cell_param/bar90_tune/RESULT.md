# BAR90 island tuning — boundary diagnosis and partial PEC150 result

## STATUS

STOPPED after one successful PEC150 eigenmode row. Fix A (bottom radiation) failed again. Fallback B (PEC floor 150 mm below bottom Cu) passed boundary readback and design validation and solved `BASE_OFF_along_80`. AEDT then rejected the required eigenmode field-grid exports needed to determine the above/below-ground field-energy split. That readback failure invokes the stated hard stop. No other tuning row started.

## GATES

| Gate | Result |
|---|---|
| Fix A, bottom radiation readback | FAIL: `AssignRadiation` rejected |
| Fix B, PEC150 boundary and design validation | PASS: five required boundaries present; design validation PASSED |
| First solve | PASS: 12 modes, saved profile, 163.8 s |
| Below-ground field-energy split | FAIL: `ExportOnGrid` rejected both tested eigenmode selectors |
| G1, all per-row readbacks | INCOMPLETE because field-energy readback failed |
| G2, clean SW branches in both directions/states | NOT EVALUABLE from one phase row |
| G3, T1–T4 rankings | NOT RUN |

## DATA

### Fix A cause and face readback

The exact exception chain in `build_abort.json` is:

> `AttributeError: Error in Executing Method InvokeAedtObjMethod.`  
> `ansys.aedt.core.internal.errors.GrpcApiError: Failed to execute gRPC AEDT command: AssignRadiation`  
> `ansys.aedt.core.internal.errors.AEDTRuntimeError: Failed to create boundary Radiation RAD_BOTTOM`

AEDT did not expose a more specific geometric cause in that traceback. A boundary-edge incompatibility is possible because the bottom perimeter shares edges with the lattice side faces; this is an inference, not a demonstrated cause.

The rebuilt Fix A face readback (`builds/BASE_A/bottom_face_inspection.json`) is:

| Item | Readback |
|---|---|
| Target face | ID **310**; parent `AIR`, material `vacuum` |
| Centre and normal | (0, 0, **−33.499**) mm within numerical precision; normal (0, 0, −1) |
| Its own boundary before assignment | None |
| Its edges | 319, 320, 321, 322 |
| Shared with lattice side faces | x faces 312 and 314 share edges 321 and 319; y faces 311 and 313 share 320 and 322 |
| Assignment order | `BAR90_02LRH` → `Lattice_X` → `Lattice_Y` → `PEC_LID` → attempted `RAD_BOTTOM` |

The cited [sub6_island_cell_packages.py:141-149](../../../scripts/research/sub6_island_cell_packages.py) is a **manual setup sheet**, with no AEDT API call or programmed face-selection method. Line 145 states the order and faces: lattice pairs on the four complete AIR side faces, top Perfect E, bottom radiation. The rebuild followed that stated order, selected the AIR bottom face by its centre, and called `assign_radiation_boundary_to_faces([310], name="RAD_BOTTOM")`; AEDT returned the same error. Thus no claim of matching an undocumented API call is made.

### Fallback B

The fresh fallback project has the same copper, RT5880 layers, diode, via placement and period as Fix A. AIR bottom moved to **−153.499 mm**, exactly 150 mm below the bottom Cu at −3.499 mm. `PEC_BOTTOM_150` was assigned to face 310; boundary readback contains `BAR90_02LRH`, `Lattice_X`, `Lattice_Y`, `PEC_LID`, and `PEC_BOTTOM_150`. The base build report records AEDT **2025.2**, 12 eigenmodes, 4.0 GHz minimum, 15 maximum passes and `Design validation check PASSED` (`builds/BASE_PEC150.json`). The approved island via local x = −1.34 mm and y centre-line remain in effect; analytic margins for all seven Round 1 configurations are in `geometry_gate.csv`.

The only solved row, `BASE_OFF_along_80`, is marked `BOTTOM=PEC150`. It took **163.797 s**; the profile says adaptive passes converged, final 9,496 tetrahedra and 5.68 GB adaptive-meshing memory. At 80°, the light-line frequency is 5.742574 GHz. Its modes below that line are:

| Mode | f GHz | Q |
|---:|---:|---:|
| 1 | 5.152092534 | 1317.54 |
| 2 | 5.524203070 | 520.04 |
| 3 | 5.740186597 | 1,117,228.65 |

There is no two-phase slope yet, so none can be labelled RES/SW/HYBRID or used for a 5.8 GHz β interpolation. All 12 raw modes are retained in `rows/BASE_OFF_along_80.json` and `modes_partial.csv`.

## FLAGS

1. The required PEC150 below-ground filter was attempted on a private copy of the saved base solution. `ExportOnGrid` failed with both `Mode=1` and the exact mode-1 frequency/0° phase selectors. The exact second traceback is in `field_split_abort.json`. No reliable field-energy split or box-mode flag was produced; modes were not silently discarded.
2. An AEDT lock from the first field-export trial remained after its own process exited. Automatic approval review rejected removal of that exact lock file, so a private copy of the solved project was used for the second read-only export attempt. No lock or other user's file was deleted.
3. No `GetMessages` call was made. The field-export trials' own AEDT process trees were stopped by exact recorded PIDs. No other AEDT session or project was touched.
4. There are no Round 1 sensitivities, Round 2 combinations or winner thickness checks. No feed or board model was changed. No `.py` file was pushed.

## FILES

`RESULT.md`, `geometry_gate.csv`, `build_abort.json`, `builds/BASE_A/bottom_face_inspection.json`, `builds/BASE_PEC150.json`, `builds/BASE_PEC150/bottom_face_inspection.json`, `rows/BASE_OFF_along_80.json`, `profiles/BASE_OFF_along_80.prof`, `modes_partial.csv`, `field_split_abort.json`.
