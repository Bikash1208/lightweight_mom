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
| `universal_database/` | Phase 8D-REV universal database architecture demo | `scripts/research/run_phase8d_rev_universal_database.py` |
| `hfss_compare_ledger.json` (not a folder) | Every HFSS-vs-MoM comparison run through the app's HFSS compare page | append-only; see `src/hfss_compare/service.py` |
| `submission_work/` | **Untracked scratch** (Oct 5–6 2026), not part of the repo | its reusable scripts moved to `scripts/adhoc/` in Phase 17; see that folder's own README |

Tests that read any of the folders above treat them as read-only fixtures/benchmarks — do not
hand-edit a results file to make a test pass; fix the generator or the test.
