# Changelog

Reconstructed per phase (this project did not keep a changelog before Phase 13; entries below
summarize each phase's own write-up in `docs/history/`/`docs/history/`). No version tag existed
before Phase 13; `0.13.0` (`pyproject.toml`) is the first. See `docs/STATUS.md` for current
validation numbers and open next steps; this file is the per-phase "what changed" record.

## Unreleased
- `hfss_exports/four_port_1bit_300mm_gap0/` and `_gap28/`: 4-port 1-bit coded boards with 300.8 mm
  sub-apertures for a minimum-gap/isolation study in HFSS (model-only predictions in each README).
- `scripts/adhoc/coded4.py`: optional `port_side_mm` argument; fine-grid (0.5 x 1 deg) per-port
  metrics and x/y-cut figure (`port_cuts.png`). No solver changes.
- `examples/make_paper_board_26ghz.py`: board-layout figure title uses the actual frequency
  (was hard-coded "26 GHz"; unchanged output at 26 GHz).

## 0.17.0 - Phase 17 - handoff-ready repo (no physics/solver changes)
- Onboarding rewrite: `README.md`, `docs/START_HERE.md` (new), `docs/STATUS.md` (new, single
  source of truth for current validation state), `CLAUDE.md` (rewritten to ≤80 lines),
  `docs/ARCHITECTURE.md` (added the three-layer + research-phase-package module map and the app's
  3-step data-flow diagram), `docs/USER_GUIDE.md` (brought up to date: stepper navigation, port
  defaults, level matching, pointing correction, Coded = model-only), `docs/validation/README.md`
  (new, indexes every validation run with numbers and provenance).
- Docs reorganized: completed phase write-ups and briefs moved to `docs/history/` (git history
  preserved via `git mv`); `docs/briefs/` keeps only active work + `_parked/`; added
  `docs/README.md`, `docs/briefs/README.md`, `docs/references/README.md`.
- `examples/` split: user-facing scripts stay; one-off research drivers moved to
  `scripts/research/`. `results/submission_work/`'s reusable scripts moved to `scripts/adhoc/`;
  scratch logs/probes removed (listed in the Phase 17 report, user-approved).
- `hfss_exports/README.md` added, indexing every package (status: awaiting HFSS / validated);
  the four untracked packages from the Phase 16 live-review round
  (`three_port_25deg`, `four_port_1bit_gap0/28`, `fss_square_loop_10mm`) committed (README/CSV/
  JSON/PNG only — `*.step` files stay gitignored alongside the already-tracked phase12/two_beam
  STEP files).
- App footer now shows its own version/commit (`HMSA Designer · <tag or short sha>`) so a
  screenshot or bug report is traceable; `make app` prints a restart-after-backend-changes hint.

## 0.16.0 - Phase 16 - follow-ups and a live-review round
- Default-on level matching (`hologram.pointing_correction.balance_beam_levels`) for any port or
  design with 2+ beams: iterates internal weights so achieved levels approach the requested ones,
  up to 3 passes, best kept; reported before/after in the beam cards.
- Optional pointing correction (off by default) wired into the single-aperture job for all three
  formulations, with the same before/after report as the multi-port path; later fixed to measure
  on a fine 0.1°-resolution local grid (the shared plotting grid was too coarse to resolve a
  small-aperture squint, which stalled correction for some beams).
- Co-pol side-lobe level (SLL) exclusion-region fix: first a grid-seam bug (a beam at phi=0/180
  read its own far shoulder as its sidelobe), then a direction-space rewrite with a nearest-beam
  (Voronoi) split so one requested beam's own main lobe is never counted as another beam's
  sidelobe.
- Multi-port layout fixes: active-aperture-side ports with the guard gap outside; the 2x2+1
  preset's large port centred correctly on the 2x2 block (was offset by gap/2) and defaulting to
  the four small ports' combined beams; the isolation-gap check tolerant of its own default-value
  rounding and reworded as information, not a fault; feed markers written in the STEP export's own
  frame and checked against each port.
- The per-port/per-beam results table's "Achieved" direction and "Level" now use the same fine
  local measurement as pointing correction and level matching, so the table and the refinement
  summary always agree.
- Restored deleted-but-tracked `hfss_exports/` files; noted in `hfss_exports/phase12_grounded_8mm/
  README.md` that the raw `.s36p` exports must be re-exported from HFSS to re-enable their tests.

## 0.15.0 - Phase 15 - polarization, vortex, multi-feed, multi-port, coded cells
- Root-cause fix for a vortex HFSS found: a monopole feed's radial surface current carries an
  orbital-angular-momentum term in circular polarization that the project's scalar pattern model
  couldn't see. Added a per-beam polarization-compensation term (natural/LHCP/RHCP/HP/VP/slant,
  or a deliberate OAM vortex request) and a genuine VECTOR far-field model via this project's own
  layered Green's function (`src/pattern/vector_aperture_pattern.py`), replacing the scalar model
  for HMSA pattern prediction.
- Multi-feed ("multipole") reference-field construction (superposition/segmentation).
- Multi-port sub-apertures (Mishra et al., IEEE TVT 2026 layout): per-port synthesis, layout
  presets, a first-order edge-distance isolation estimate.
- Coded (1-bit/N-bit) design method: quantizes the phase target to a small number of reflective
  states instead of a continuous value.
- Corresponding Step 2/3 app UI for all of the above, verified live in a real browser.

## 0.14.0 - Phase 14 - surface-wave eigenmode dispersion
- A direct MoM eigenmode search over the periodic array's own Z(kx,ky) was tried and found not to
  resolve a loading-dependent mode at any affordable mesh size (kept as a documented research
  item, `src/dispersion/mom_scaffold.py`).
- Working replacement: `src/dispersion/sheet_model.py` extracts a patch-sheet admittance from the
  already-validated normal-incidence S11 and solves the closed-form transverse-resonance
  condition for beta. Added a "Design method" choice (surface-wave/eigenmode vs. reflection-phase)
  to Step 2 — on the walk-through case the eigenmode method reaches all 5476 cells vs. the
  reflection-phase method's 4006 unreachable.

## 0.13.0 - Phase 13 - production repo + local web app
- Part 0 carry-overs (`docs/history/PHASE13_WEBAPP.md`): Richardson mesh-convergence extrapolation
  at resonance; resonance-region error re-stated for Phase 12B/12C/12D's own meshes; G4/G6 real
  HFSS ingestion with a mesh/thickness/residual error budget; the missing lookup-table-predicted-
  vs-measured per-order-power comparison for Phase 12D's Task B.
- Part 1: git repository initialized; `pyproject.toml` packaging (`pip install -e .[dev,app,cad]`);
  docs moved into `docs/{theory,phases,briefs}/` plus `docs/history/PROJECT_HISTORY.md`; `configs/` for
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
