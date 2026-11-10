# Fixed BASELINE / G1 mesh trend and HFSS cell packages

MODEL_ONLY / NOT_HFSS_VALIDATED. Diagnostic mesh sweep only; no geometry tuning.

| Cell | Mesh | gap cells / mm | alpha OFF/ON Np/m | suppression dB | beta OFF/ON /k0 | extrapolated OFF/ON |
|---|---:|---:|---:|---:|---:|---|
| BASELINE | 96 | 8 / 0.966667 | 2.012370 / 0.225305 | 9.5094 | 1.0940043 / 1.0167730 | False / True |
| BASELINE | 128 | 11 / 0.996875 | 2.154554 / 1.807468 | 0.7629 | 1.0965649 / 1.2433828 | False / True |
| BASELINE | 160 | 14 / 1.015000 | 2.049919 / 0.136760 | 11.7578 | 1.0942488 / 1.0107136 | False / True |
| BASELINE | 192 | 17 / 1.027083 (input) | ALLOCATION_FAILED | — | — | — |
| G1 | 96 | 8 / 0.966667 | 2.824425 / 0.033603 | 19.2455 | 1.1119022 / 1.0035652 | True / True |
| G1 | 128 | 11 / 0.996875 | 3.075198 / 0.092803 | 15.2031 | 1.1246478 / 1.4297169 | True / True |
| G1 | 160 | 14 / 1.015000 | 2.699653 / 0.627012 | 6.3403 | 1.1070242 / 1.3455619 | True / True |
| G1 | 192 | 17 / 1.027083 (input) | ALLOCATION_FAILED | — | — | — |

The physical diode gap remains 1.00 mm. Nearest-cell discretization gives 8/11/14/17 cells at meshes 96/128/160/192. All Floquet cutoffs are 2n; the original F128 baseline is excluded from this standardized trend.

Bloch nodes remain kt/k0=1.017 and 1.097 in each direction. Beta outside that interval is EXTRAPOLATED, not validated by a direct solve at that beta. ON-state instability is intentionally reported and is not a stop condition for this diagnostic job.

The table alpha is first-order holographic leakage. It is NOT alpha_x from a uniform-cell eigenmode Q. No physical convergence or HFSS validation is claimed.

Four checked HFSS packages: hfss_exports/sub6_island_baseline_g1/README.md. Each contains one shared OFF/ON STEP, one-page manual setup, state table, eight-mode phase-run template, predictions and geometry gates. All metal is one METAL_PEC object; one 1.00 x 0.35 mm Impedance sheet; no AEDT scripts.

Calibration method found: scripts/research/sub6_final_cal.py; actual HFSS record results/sub6_final_board/hfss_data.json; archived setup commit 21fe444:hfss_exports/sub6_island_choke/HFSS_SETUP.md. The recorded calibration was ISL_vias2, not no-via. It used periodic eigenmode f-vs-phase and Q to obtain beta and alpha_x=omega/(2Qvg). Historical ON slope used MoM anchoring; new packages request a directly measured local HFSS slope from >=3 phase points.

User-confirmed variant (b) retains both 2.0 x 0.3 mm bottom stubs, original via centres, 0.8-mm FR4, 0.6-mm anti-pads, distant radiation bottom 30 mm below the lowest copper. +90-degree frame: across=ph_x, along=ph_y. PEC lid at +150 mm. No near-stub PEC floor.

Expected directional no-via beta/loss estimates and their extrapolation flags are in each model_reference.csv/predictions.json. Actual 3D via/FR4/stub values are NOT_AVAILABLE_FROM_CURRENT_DATA. No-via holographic leakage values are not substituted for Q loss. No HFSS solve was performed by S2.

Geometry differences are documented: BASELINE uses physical board 8.696 x 2.900 slot and 4.348 x 1.300 island; G1 CAD endpoints are snapped to 0.001 mm. MoM grid realizations and actual gap lengths are separately reported. G1 via barrel fits with only 0.003 mm x-edge clearance; no fabrication clearance qualification is claimed.

A tool-helper refresh interrupted the first mesh-160 attempt; completed nodes were retained and only missing nodes restarted in a hidden detached process. No parameter or worker-count change was made. Worker pools use all 64 logical workers, single-threaded BLAS; no worker cap.

Optional mesh192 was attempted with eight nodes and no worker cap. The process failed before any mesh192 checkpoint completed: NumPy could not allocate a 22.9 GiB complex128 array of shape 39220 x 39220. No mesh192 number is inferred. See mesh192.err and mesh192_status.json. This is an execution-resource failure, not a physical ON-state result or an imposed worker/RAM budget.

Desktop delivery: C:/Users/Bikash/Desktop/sub6_board_hfss/iter1_cell_checks_S2_2026-11-10. Only STEP/JSON/CSV/MD copied; no Python files. Existing user-modified delivery files are preserved.

Validation: four geometry packages passed all pre/re-import gates; eight CAD-to-model quantization mappings passed. Full fast suite: 2953 passed, 129 skipped, six accepted ENV_ONLY hash failures, two 60-second app polling timeouts. Both timeout tests passed isolated rerun (2 passed in 41.42 s). Targeted geometry/FFT checks: 5 passed. No src/ or existing tests/fixtures changes.
