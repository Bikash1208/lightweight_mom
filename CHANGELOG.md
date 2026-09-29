# Changelog

Reconstructed per phase (this project did not keep a changelog before Phase 13; entries below
summarize each phase's own write-up in `docs/phases/` and `docs/PROJECT_HISTORY.md`). No version
tag existed before Phase 13; `0.13.0` (`pyproject.toml`) is the first.

## 0.13.0 - Phase 13 - production repo + local web app
- Part 0 carry-overs (`docs/phases/PHASE13_WEBAPP.md`): Richardson mesh-convergence extrapolation
  at resonance; resonance-region error re-stated for Phase 12B/12C/12D's own meshes; G4/G6 real
  HFSS ingestion with a mesh/thickness/residual error budget; the missing lookup-table-predicted-
  vs-measured per-order-power comparison for Phase 12D's Task B.
- Part 1: git repository initialized; `pyproject.toml` packaging (`pip install -e .[dev,app,cad]`);
  docs moved into `docs/{theory,phases,briefs}/` plus `docs/PROJECT_HISTORY.md`; `configs/` for
  config-file presets; ruff lint/format; pytest markers (`slow`/`cad`/`checkpoint`/`hfss_data`) and
  a fast (<5 min) suite; GitHub Actions CI (written, not run); `requirements-lock.txt`; `Makefile`.
- Part 2: structured logging (`LWMOM_LOG_LEVEL`, per-run ID), solver diagnostics capture, a typed
  error-code taxonomy, a fixed validation-label vocabulary, and a run manifest (inputs, resolved
  config, code version, timings, labels, warnings, output hashes) attached to every job.
- Part 3: local FastAPI web app (`make app`) covering unit cell / library / reflectarray / pattern
  prediction / verify / export / HFSS compare / HMSA designer pages; a real HFSS field-calculator
  `.fld` reference-field importer; an HFSS-eigenmode-table dispersion-library import
  (`DISPERSION_FROM_HFSS`) with a documented, tested MoM-eigenmode scaffold
  (`DISPERSION_FROM_MOM`, not implemented) left for later; both hologram formulations from Phase
  8B wired into an HMSA designer; job queue with progress/cancel/config-hash caching.

## 0.0.0 (unreleased, pre-Phase-13) - Phases 1-12D
- **Phases 1-5**: core periodic MoM solver -- geometry/config/units, dyadic layered Green's
  function, rooftop MoM (`Z I = V`), S-parameter extraction at configurable reference planes,
  circular-polarization basis and the unit-cell database schema.
- **Phase 6/6A**: STEP export of representative unit cells; independent full-wave-solver
  validation-export packages (no EM code changed, no automatic HFSS agreement claimed).
- **Phase 7**: production 324-geometry x 201-frequency unit-cell library with a query/phase-
  inversion API (resumable, checkpointed).
- **Phase 8A-8E**: holographic-synthesis feasibility gate; object-wave/multi-beam formulation;
  unit-cell scattering architecture validation; universal unit-cell database architecture; the
  physical control-variable/holographic-scattering bridge (reference-wave -> object-wave ->
  modulation -> control-variable pipeline).
- **Phase 9**: control-variable and local-periodicity validation against real heterogeneous
  supercell solves.
- **Phase 10**: environment-aware neighborhood/surrogate/inverse/beam validation.
- **Phase 11**: physics-based heterogeneous-MoM nonlocal-coupling study (no ML/fitted surrogates).
- **Phase 12A**: optional PEC ground plane (`backing="pec"`), additive to the Green's function and
  config schema.
- **Phase 12B**: mesh alignment / duplicate-geometry audit / grounded period study.
- **Phase 12C**: grounded-supercell local-periodicity/coupling/fine-mesh spot-check; HFSS-ready
  STEP exports.
- **Phase 12D** (+ addendum): design-relevant local-periodicity step-limit re-measurement;
  aperture-level Floquet-order-power test (decided Phase 13's forward-model recommendation); real
  HFSS Touchstone (`.sNp`) ingestion (G1); mesh non-convergence at resonance discovered and
  handed to Phase 13 Part 0 as a carry-over.
