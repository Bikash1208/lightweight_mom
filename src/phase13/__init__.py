"""Phase 13 -- production repo + local web app.

Part 0 (carry-overs from Phase 12D / PHASE12D_ADDENDUM.md, this package):
  - `convergence.py`   -- Richardson extrapolation of the MoM solve in mesh cell size
                          (PHASE12D_ADDENDUM.md task 1). Reuses `sparameters.compute_s_parameters`
                          unmodified at a family of aligned mesh resolutions; does not touch
                          `solver.py`/`sparameters.py`.
  - `resonance_restatement.py` -- quantifies the 0.2 mm-mesh resonance-region error in the
                          Phase 12B/12C 8 mm results the addendum flagged (task 2). Read-only:
                          writes new numbers here, never edits earlier phases' outputs.
  - `hfss_resonance_ingest.py` -- ingests G4 (`g4.s36p`, zero-thickness sheet run) and G6
                          (`g6.s36p`) alongside G1 (task 3), comparing HFSS to both the fixed
                          0.2 mm MoM and the Richardson-converged MoM.
  - `order_power_prediction.py` -- the missing comparison PHASE12D's Task B left out: lookup-table
                          *predicted* per-Floquet-order power (from the isolated-cell curve, via
                          local periodicity) vs the full-supercell MoM ground truth, for
                          N = 6/8/12, so the aperture-level finding can be attributed to the design
                          point (186 deg deg achieved phase span) rather than to the lookup-table
                          forward model.
  - `gate.py`, `run.py`, `report.py`, `specification.py` -- orchestration, exactly the pattern used
    by every previous phase.

Part 2 (logging/diagnostics/errors/labels/manifest) and Part 3 (the web app, field import,
dispersion library, HMSA/reflectarray/pattern/export/verify/compare services) live in their own
sibling packages (`src/logging_setup.py`, `src/diagnostics.py`, `src/errors.py`, `src/labels.py`,
`src/manifest.py`, `src/fields/`, `src/dispersion/`, `src/app/`) -- see
`docs/phases/PHASE13_WEBAPP.md` for the full map. Everything under this package (`src/phase13/`)
is analysis-only, following the project's "new phases reuse earlier phases unmodified" rule; it
does not edit any earlier-phase or core module.
"""
