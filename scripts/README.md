# Scripts

Not for day-to-day use — see `examples/` for that. This folder holds the one-off drivers that
built this project's own validation record, kept for reproducibility and reference, not as a
supported CLI surface.

## `scripts/research/` — the per-phase CLI each phase's write-up (`docs/history/`) was produced from

| Script | Phase |
|---|---|
| `run_phase8a_feasibility.py` | 8A — scientific feasibility gate against the Phase 7 library |
| `run_phase8b_formulation.py` | 8B — holographic formulation / control-variable question |
| `run_phase8c_validation.py` | 8C — unit-cell scattering architecture validation |
| `run_phase8d_rev_universal_database.py` | 8D-REV — universal unit-cell database architecture |
| `run_phase8e_control_bridge.py` | 8E — control-variable / holographic-scattering bridge |
| `run_phase9_validation.py` | 9 — control-variable + local-periodicity validation |
| `run_phase10_environment_model.py` | 10 — environment-aware electromagnetic response model |
| `run_phase11_heterogeneous_mom.py` | 11 — heterogeneous MoM / nonlocal coupling |
| `run_phase12a_ground_plane.py` | 12A — optional PEC ground plane |
| `run_phase12b.py` | 12B — mesh alignment, duplicate-geometry audit, grounded period study |
| `run_phase12c.py` | 12C — grounded 8 mm local-periodicity/coupling/spot-check |
| `run_phase12d.py` | 12D — local-periodicity step limit + real HFSS ingestion |
| `run_phase13_part0.py` | 13 Part 0 — carry-overs from the Phase 12D addendum |
| `export_validation_cases.py` | 6 — exports the 10 validation-case STEP/JSON/CSV folders |
| `export_hfss_comparison.py` | 6A — the HFSS comparison export for an already-exported case |

Each writes to its own `results/<phase_name>/` folder (see `results/README.md`). Run from the
project root, e.g. `.venv/bin/python scripts/research/run_phase12d.py`.

## `scripts/adhoc/` — rescued scripts from an untracked working folder

Reusable scripts kept from `results/submission_work/` (an untracked Oct 5–6 scratch folder) that
generated committed outputs still referenced elsewhere. See `scripts/adhoc/README.md` for exactly
what each one does and its status.
