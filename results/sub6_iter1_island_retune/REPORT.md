# BAR50-02V island retune (S2, 2026-11-10)

MODEL_ONLY / NOT_HFSS_VALIDATED. One diode; physical gap 1.00 mm.

All numbers below are computed with the confirmed series-L circuit.

| ID | mesh | alpha OFF/ON Np/m | suppression dB | beta OFF/ON /k0 | delta beta/k0 | OFF-corner suppression min/max dB |
|---|---:|---:|---:|---:|---:|---:|
| G1 | 96 | 2.824425 / 0.033603 | 19.2455 | 1.1119022 / 1.0035652 | -0.1083370 | 19.2219 / 19.2904 |
| G1 | 128 | 3.075198 / 0.092803 | 15.2031 | 1.1246478 / 1.4297169 | 0.3050691 | 15.1868 / 15.2316 |
| G2 | 96 | 2.290183 / 0.091986 | 13.9615 | 1.0982812 / 1.0073802 | -0.0909010 | 13.9260 / 14.0302 |
| G2 | 128 | 2.357577 / 0.502142 | 6.7164 | 1.0998624 / 1.3596743 | 0.2598119 | 6.6737 / 6.7980 |
| G3 | 96 | 2.006898 / 0.088151 | 13.5730 | 1.0937846 / 1.0069865 | -0.0867981 | 13.5337 / 13.6495 |
| G3 | 128 | 2.133791 / 1.562447 | 1.3535 | 1.0959997 / 1.2648079 | 0.1688081 | 1.3081 / 1.4403 |

Status: **HARD_STOP_MESH_NOT_CONVERGED**.

Dimensions below are actual mesh-96 dimensions; mesh-128 realizations are in mesh_report.json.

| ID | slot L x W mm | island L x W mm | bottom/end clearance mm | diode x mm |
|---|---:|---:|---:|---:|
| G1 | 9.183333 x 2.658333 | 3.866667 x 1.087500 | 0.604167 / 2.658333 | 0.000000 |
| G2 | 8.941667 x 2.658333 | 4.350000 x 1.087500 | 0.604167 / 2.295833 | 0.000000 |
| G3 | 8.700000 x 2.900000 | 3.866667 x 1.329167 | 0.604167 / 2.416667 | 0.000000 |

G1 maximum relative alpha change: 176.176%; HARD_STOP_ALPHA_CHANGE_GT_5_PERCENT.
G2 maximum relative alpha change: 445.889%; HARD_STOP_ALPHA_CHANGE_GT_5_PERCENT.
G3 maximum relative alpha change: 1672.462%; HARD_STOP_ALPHA_CHANGE_GT_5_PERCENT.

Physical intent uses gap exactly 1.00 mm and island width = slot width - bottom clearance - 1.00 mm; these widths differ from the grid-realized widths above. See physical_intent_mm in mesh_report.json. No STEP is exported from these mesh masks.

Sensitivity: C={0.105,0.12} pF crossed with Rp={3000,5000} ohm; ON is unchanged.
Baseline nominal C=0.11 pF, Rp=4000 ohm. Both states include L=0.6 nH; ON R=3 ohm.
All cases use Floquet order 2n; accepted gap realization 8 cells (96), 11 cells (128).
Sheet extraction uses the existing two-node self-consistent mean-sheet/first-order leakage model; raw tensors and extrapolation flags are preserved.

Resonance/stop-band clearance is NOT_ESTABLISHED: frequency screening was not launched after the mandated >5% mesh-change hard stop.

No candidate is qualified for a board or HFSS export. The table contains raw two-node model outputs, not validated propagation constants.

Post-stop sheet audit (existing tensors only): nominal fixed-point residuals are <=1e-5, but all ON states extrapolate outside the 1.017..1.097 kt/k0 node interval. Mesh-128 ON beta/k0=1.265..1.430 is far outside that interval; direct Bloch validation at those values has not been run. See sheet_iteration_audit.json.

Validation: production baseline alpha/beta reproduced within 2.12e-12 relative; existing FFT assembly audit <=1.03e-15 relative sampled matrix error; exact load-update/direct-solve tensor audit <=1.08e-12. No src/ solver changes.

Fast suite in previous S2 env: 2954 passed, 129 skipped, six accepted ENV_ONLY S11 hash failures, one 300-second app timeout. App test rerun passed in 68.55 seconds; four new geometry/circuit checks passed. Original regression tests and fixtures unchanged. The initially reported four ENV_ONLY failures are documented in ENV_DIAGNOSIS.md.

HFSS_SETUP.md and PACKAGE_9CELL_PLAN.md corrected: OFF sheet impedance 5.424-j79.320 ohm/sq includes series 0.6 nH; former -86.97 omitted it. ON remains 1.05+j7.65. No prior HFSS run is relabelled.

Original C:/lwmom_final modified CSV: only sheet_z_thick column moved after sheet_centre_y; all 240 rows' named values unchanged. Original artifacts preserved.
