# Results

Generated outputs, one folder per phase/script. Tracked research outputs (small
JSON/CSV/MD/plots/STEP) stay in git as the project's own record; a few genuinely heavy, fully
regenerable folders are gitignored (noted below).

| Folder | Produced by | Regenerate with |
|---|---|---|
| `phase8a_feasibility/` | Phase 8A feasibility gate | `scripts/research/run_phase8a_feasibility.py` |
| `phase8b_formulation/` | Phase 8B holographic formulation | `scripts/research/run_phase8b_formulation.py` |
| `phase8c_scattering/` | Phase 8C scattering architecture validation | `scripts/research/run_phase8c_validation.py` |
| `phase8e_control_bridge/` | Phase 8E control-variable bridge | `scripts/research/run_phase8e_control_bridge.py` |
| `phase9_validation/` | Phase 9 control-variable/local-periodicity validation | `scripts/research/run_phase9_validation.py` |
| `phase10_environment_model/` | Phase 10 environment-aware model | `scripts/research/run_phase10_environment_model.py` |
| `phase11_heterogeneous_mom/` | Phase 11 heterogeneous-MoM coupling study | `scripts/research/run_phase11_heterogeneous_mom.py` |
| `phase12a_ground_plane/` | Phase 12A PEC ground plane | `scripts/research/run_phase12a_ground_plane.py` |
| `phase12b_mesh_period_study/` | Phase 12B mesh/period study | `scripts/research/run_phase12b.py` |
| `phase12c_grounded_locality/`, `phase12c_hfss_export/` | Phase 12C grounded-supercell study + HFSS exports | `scripts/research/run_phase12c.py` |
| `phase12d_local_periodicity/` | Phase 12D step-limit + HFSS ingestion | `scripts/research/run_phase12d.py` |
| `phase13_part0/` | Phase 13 Part 0 carry-overs | `scripts/research/run_phase13_part0.py` |
| `validation_cases/` | Phase 6 STEP/JSON/CSV validation-case exports | `scripts/research/export_validation_cases.py` |
| `unit_cell_library/` | Phase 7 production unit-cell library | `examples/build_unit_cell_library.py` (long-running, resumable). **Gitignored**: `_solve/`, `database/`, `geometry/` (~554 MB) — only `manifest.json`/reports/plots are committed. |
| `nylon_cells_28ghz/` | Phase 19 Part B: 28 GHz nylon octagon cell library (RCWA, model-only) | `examples/build_nylon_cell_library_28ghz.py` (~3 min) |
| `transmitarray_dual30/` | Phase 19 Parts C-D: dual-beam ±30° nylon lens design, local-periodicity and error-bar records (model-only) | `scripts/research/phase19_design.py` (~15 min), then `phase19_local_periodicity.py` (~10 min); see its README |
| `phase20_r1_demo/` | Phase 20 R1: Oblique-incidence and lumped-load MoM demo outputs and convergence checks | `scripts/research/phase20_r1_demo.py` |
| `phase20_r2/` | Phase 20 R2: PIN diode circuit models (Skyworks SMP1345-079LF) and bias-dependent reflection sweeps | `scripts/research/phase20_r2_demo.py` |
| `phase20_r3/` | Phase 20 R3: 2-bit Pancharatnam-Berry octagon unit cell (CP-2) optimization, 4-state phase ladder, and supercell locality data | `scripts/research/phase20_r3_cp_package.py`, `scripts/research/phase20_r3_opt.py`, `scripts/research/phase20_r3_locality.py` |
| `phase20_r3b/` | Phase 20 R3b: Re-optimization study outcome; baseline CP-2 cell retained, optimisation stopped (`OUTCOME.md`) | `scripts/research/phase20_r3b_cp_package.py` |
| `phase20_r4/` | Phase 20 R4: Discrete state library (v1 and full-grid v2 at Floquet cutoff F = 2n = 120), cutoff audit, and HFSS compare plots. **Gitignored**: `chunks/` | `scripts/research/phase20_r4_build_library.py`, `scripts/research/phase20_r4_assemble.py`, `scripts/research/phase20_r4_hfss_compare.py` |
| `phase20_r5/` | Phase 20 R5: RIS design mode records, reference-plane calibration (`PLANE_CORRECTION.md`), R5-3 beam tables (`R5_3_REPORT.md`), and server MoM tolerance analysis (`TOLERANCE_GAP.md`) | `scripts/research/phase20_r5_rebuild_delta.py`, `scripts/research/phase20_r5_tolerance_mom.py`, `scripts/research/phase20_r5_tolerance_compare.py` |
| `phase23/` | Phase 23: square vs hexagon vs hexagon + sidelobe optimiser comparison on the held-out grid; the rectangle bit-identical gate | `scripts/research/phase23_compare_v2.py`; gate `scripts/research/phase23_rect_snapshot.py` |
| `phase25/` | Phase 25: Butler-fed sector scanner forward model (`forward_validation.json`), supercell locality checks, and aperture-taper sensitivity (model-only) | `scripts/research/phase25_forward_validation.py`, `scripts/research/phase25_locality_package.py`, `scripts/research/phase25_block_locality.py` |
| `phase26/` | Phase 26: Leakage-controlled tensor-CP hexagon design (`fourhex_A_cells.json`, `fourhex_B_cells.json`), launcher calibration, and substrate screening tables | `scripts/research/phase26_fourhex_design.py`, `scripts/research/phase26_launcher_package.py`, `scripts/research/phase26_substrate_study.py` |
| `infra/` | Server sync provenance records, mirror checks, and cross-machine equivalence artifacts | `scripts/research/equivalence_compare.py`, `scripts/research/phase20_r4_mirror_check.py` |
| `universal_database/` | Phase 8D-REV universal database architecture demo | `scripts/research/run_phase8d_rev_universal_database.py` |
| `hfss_compare_ledger.json` (not a folder) | Every HFSS-vs-MoM comparison run through the app's HFSS compare page | append-only; see `src/hfss_compare/service.py` |
| `submission_work/` | **Untracked scratch** (Oct 5–6 2026), not part of the repo | its reusable scripts moved to `scripts/adhoc/` in Phase 17; see that folder's own README |

Tests that read any of the folders above treat them as read-only fixtures/benchmarks — do not
hand-edit a results file to make a test pass; fix the generator or the test.
