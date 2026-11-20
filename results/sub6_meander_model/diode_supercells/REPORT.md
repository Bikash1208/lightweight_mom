# Diode-cell rotated supercells — stopped at ON-A

## STATUS

**ABORTED before a converged row.** No diode performance conclusion is available. The ordered queue stopped at ON-A; ON-B, SHORT, SHORT unrotated, OFF-D rotated and OFF-D unrotated were not launched. No board solve was run.

## GATES

| Row / attempt | AEDT 2026 R1 outcome | Mode result |
|---|---|---|
| ON-A, series lumped RLC | `AssignLumpedRLC` failed during boundary assignment before solve: `AEDTRuntimeError: Failed to create boundary Lumped RLC RLC_0`. | None |
| ON-A, impedance sheet retry | Build, validation and boundary readback passed; own solve exceeded 3600 s and its exact recorded AEDT PID was stopped by the watchdog. | No converged mode list; no alpha |

The impedance sheet is 0.3 mm along current × 0.5 mm wide. Its readback is `Resistance = 3.3333333333333335`, `Reactance = -232.66666666666666` ohm/square, namely ON-A element `2.0 - j139.6` ohm multiplied by 0.5/0.3. All eight sheets have impedance boundaries. The 0.5-mm strip, two 4.11-mm slots and 0.3-mm bridge gap were built in all eight cells with rotations 0° to 157.5°; each copper area read back as about 69.24 mm². The substrate, 40.035-mm air/PML and lattice readbacks passed. No settings or geometry were changed between the RLC attempt and the impedance retry.

The retained partial HFSS profile reached adaptive pass 6 at 27,080 tets, max delta frequency 1.10471%, with 44 cores granted. This is **not** a converged solution. The watchdog recorded `TIMEOUT_3600S` and terminated only the four process IDs belonging to the ON-A attempt. The controller subsequently failed to query the terminated AEDT process (`AttributeError: 'NoneType' object has no attribute 'objectID'`). Both AEDT PIDs from this job were verified gone. The earlier RLC attempt's AEDT process did not close after the API error and was terminated by its exact recorded PID.

## DATA

The anchor reference remains 5.6817 GHz, Q 67.4 at supercell phase 30° (earlier AEDT 2025 R2). The diode cell cannot be compared with it because ON-A produced no converged modes. No PASS, SHIFTED or FAIL physics label is assigned. K is undefined.

## FLAGS

The preferred lumped RLC API failed before the solver, so the permitted impedance-sheet alternative was used. The raw controller traceback and partial solver profile are preserved. No separate AEDT message-log file was emitted by this failed solve; the prior instruction forbidding `GetMessages` was respected. The old anchor reference was solved in AEDT 2025 R2, whereas this attempt ran in AEDT 2026 R1.

## FILES

`ON_A_RLC_attempt1_controller.log`, `ON_A_RLC_attempt1_status.json`, `ON_A_IMP_attempt2_controller.log`, `ON_A_IMP_attempt2_status.json`, `ON_A_IMP_attempt2_build.json`, `ON_A_IMP_attempt2_status_external_heartbeat.json`, `ON_A_IMP_attempt2_partial.profile`, `queue_status.json`, and `queue.done` contain the record. The private AEDT project remains local and is not pushed.
