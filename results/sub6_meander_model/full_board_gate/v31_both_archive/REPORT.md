# v3.1 BOTH — model-only AEDT archive

**STATUS:** COMPLETE. Archived a copy of the converged v3.1 BOTH model in HFSS 2026 R1. No solve or geometry edit was run.

**GATES:** The source convergence record says normal completion after six adaptive passes. Copied model validation passed; 518 PEC metal objects and total metal volume were unchanged. CAD import history was purged from 522 imported bodies, leaving zero STEP references in the saved project and archive. The archive contains the AEDT project and its preview only, with no simulation results or external STEP. ZIP CRC passed.

**DATA:** Local saved project `C:\lwmom_meander_model\hfss_work\full_board_gate_v31_both_archive\BOARD_R150_V31_BOTH.aedt` (6,107,726 bytes). Published `BOARD_R150_V31_BOTH.aedtz` (677,593 bytes). Port1 is a 50-ohm lumped annular sheet at z = −3.175 mm: outer radius 2.05 mm, nominal inner radius 0.635 mm (faceted area readback implies 0.6363 mm). Read-back integration line: (0.635, 0, −3.175) to (2.05, 0, −3.175) mm; renormalization 50 ohm, no de-embedding. Setup1 is 5.68 GHz, direct solver, basis order 1, max 12/min 6 passes, max ΔS 0.02. Mesh and all boundaries were retained.

**FILES:** `BOARD_R150_V31_BOTH.aedtz`, `audit.json`, `archive_status.json`; the `.aedt` remains local. `audit.json` includes raw port and setup readbacks.
