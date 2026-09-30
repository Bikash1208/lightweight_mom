# CLAUDE.md

Guidance for Claude when working in this repository.

## What this is

`lightweight_mom` — a lightweight Python Method-of-Moments solver for **periodic metasurface unit cells**
(not a general full-wave solver), plus a phased research pipeline toward holographic metasurface synthesis:

```
geometry → periodic MoM (Z I = V) → surface current → fields → S11/S21 → unit-cell database → holographic synthesis
```

Status: Phases 1–11, 12A, 12B, 12C, 12D, 13, 14 and 15 done. Phase 13 curated the repo to production-push grade (git,
`pyproject.toml`, `docs/` layout, logging/diagnostics/error-taxonomy/validation-labels/run-manifest, ruff, CI) and
built a local web app (FastAPI backend, static HTML/vanilla-JS frontend) that produces reflectarray/HMSA designs NOW,
labelled honestly, rather than waiting on the open questions below. A live UI review (`docs/briefs/
PHASE13B_UI_REVIEW_FIXES.md`) then a redesign brief (`docs/briefs/PHASE13C_REDESIGN.md`) rebuilt the app around one
3-step HMSA wizard (home page) plus a "More tools" page for the rest, added server-rendered PNG plots for the
wizard's own flow (still tables elsewhere — a documented, narrower scope gap now, see
`docs/phases/PHASE13_WEBAPP.md` "Part 3B/3C"), and fixed several correctness bugs the review found, including one
real physics bug (Formulation B's predicted aperture field was missing the reference-wave phasor). A second live
walk-through (`docs/briefs/PHASE13D_POLISH_AND_PHYSICS.md`) then made `src/dispersion/mom_scaffold.py`'s
surface-wave eigenmode search real (validated against the analytic bare-slab TM0), fixed a two-cause pattern
side-lobe bug (a missing obliquity factor, and `detect_local_maxima` not wrapping the periodic phi axis), made STEP
export opt-in with estimated progress, and added a real DXF export — see `docs/phases/PHASE13_WEBAPP.md` "Part 3D".
Phase 14 (`docs/phases/PHASE14_EIGENMODE.md`) then built a real surface-wave dispersion source: a direct MoM
eigenmode search over the periodic array's own Z(kx,ky) was tried and found not to resolve a loading-dependent mode
at any affordable mesh size (kept as a documented research item, `src/dispersion/mom_scaffold.py`); the working
replacement (`src/dispersion/sheet_model.py`) extracts a patch-sheet admittance from the already-validated
normal-incidence S11 and solves the closed-form transverse-resonance condition for beta. The app now offers a
"Design method" choice (surface-wave/eigenmode vs. reflection-phase) in Step 2; on the walk-through case the
eigenmode method reaches all 5476 cells vs. reflection-phase's 4006 unreachable.
Phase 15 (`docs/phases/PHASE15_POLARIZATION_VORTEX.md`) fixed a root-cause physics gap the user found in HFSS: every
monopole-fed design radiated a vortex (doughnut) beam because the feed's radial surface current, decomposed into
circular polarization, carries an OAM term the app's scalar pattern model couldn't see. Added a per-beam
polarization-compensation term in hologram synthesis (Beam.desired_polarization/oam_charge — natural/LHCP/RHCP/HP/
VP/slant, plus deliberate OAM vortex requests; `src/hologram/polarization_compensation.py`), a genuine VECTOR
far-field model via this project's own layered Green's function (`src/pattern/vector_aperture_pattern.py`, replacing
the scalar Fourier-transform model for HMSA — measures null depth + phase winding, not just assumes vortex vs
pencil), multi-feed ("multipole") reference-field construction (superposition/segmentation,
`src/hologram/multi_feed.py`), and the corresponding Step 2/3 app UI, verified live in a real browser.
Holographic synthesis, geometry inversion and full-array STEP export remain recorded as **LOCKED** scientific
findings (see `README.md` "Before Phase 13's app is trusted for a real design → still open") even though the app now
produces designs/exports anyway per the user's explicit Phase 13 decision — every such output carries these labels.

## Commands

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev,app]"        # add ",cad" for STEP export/plotting
.venv/bin/pip install -r requirements-lock.txt                        # OR: reproduce the exact pinned environment
.venv/bin/python -m pytest -m "not slow" -q                           # fast suite (<5 min), run from project root
.venv/bin/python -m pytest -q                                         # full suite (tens of minutes)
.venv/bin/python -m pytest tests/test_mom_matrix.py -q                # single file
make app                                                              # starts the local web app (FastAPI, localhost)
.venv/bin/python examples/print_config.py                             # load configs/config.yaml, print in mm/GHz
.venv/bin/python examples/run_sweep.py                                # sweep → results/sweep_example (resumable)
.venv/bin/python examples/run_phase11_heterogeneous_mom.py            # each phase has examples/run_phaseN_*.py
.venv/bin/python examples/run_phase12a_ground_plane.py                 # Phase 12A: optional PEC ground plane (backing=pec)
.venv/bin/python examples/run_phase12b.py                               # Phase 12B: mesh alignment, duplicate audit, period study
.venv/bin/python examples/run_phase12c.py                               # Phase 12C: grounded supercells, local-periodicity/coupling/spot-check, HFSS exports
.venv/bin/python examples/run_phase12d.py                               # Phase 12D: step-limit re-measurement, aperture Floquet-power test, real HFSS ingestion
.venv/bin/python examples/run_phase13_part0.py                          # Phase 13 Part 0: convergence, resonance restatement, HFSS G4/G6, model-vs-design
```

- Optional extras (`pyproject.toml`): `cad` (cadquery, matplotlib — STEP export, plots, Phase 6/6A/7 modules,
  **never imported from the solver core**), `app` (fastapi, uvicorn — the Part 3 web app), `dev` (pytest, ruff,
  pre-commit, httpx). `requirements.txt`/`requirements-cad.txt` still exist for backward compatibility but
  `pyproject.toml` extras are the source of truth going forward.
- Pytest markers (`pytest.ini`): `slow` (real MoM solves; excluded from the fast suite via `-m "not slow"`), `cad`,
  `checkpoint`, `hfss_data` (the last three self-skip when their optional dependency/data file is absent — see
  `docs/DEVELOPING.md` for what's excluded from git and how to regenerate or obtain it).
- Tests using the full Phase 7 checkpoint (`results/unit_cell_library/_solve/database.json`, ~470 MB, ~20 s load,
  excluded from git) go through the session fixture in `tests/conftest.py` and skip if it's absent.
- `examples/build_unit_cell_library.py` is long-running (324 geometries × 201 freqs); runs are resumable via
  append-only `checkpoint.jsonl`.

## Read before touching code

- **`docs/theory/CONVENTIONS.md` is the single source of truth.** If code and it disagree, that is a bug: fix one and
  add a test. Key points: SI internally; `exp(+jωt)`; `Im(kz) ≤ 0`; `G = −Z`; `S_ij[received, incident]`, index 0 = x,
  1 = y, same global axes at both ports; phase stored in degrees on (−180, 180]; IEEE right-hand-rule circular naming.
- Derivations: `docs/theory/GREEN_FUNCTION.md`, `docs/theory/MOM_FORMULATION.md`, `docs/theory/S_PARAMETERS.md`,
  `docs/theory/DATABASE.md`. Per-phase docs (`docs/phases/`): `PHASE8A_…` through `PHASE11_HETEROGENEOUS_MOM.md`,
  `PHASE12A_GROUND_PLANE.md`, `PHASE12B_MESH_PERIOD_STUDY.md`, `PHASE12C_GROUNDED_LOCALITY.md`,
  `PHASE12D_LOCAL_PERIODICITY.md`, `UNIT_CELL_LIBRARY.md`, `UNIVERSAL_UNIT_CELL_DATABASE.md`,
  `VALIDATION_EXPORT.md`, and now `PHASE13_WEBAPP.md`. Original task briefs: `docs/briefs/` (`*_PROMPT.md`,
  `PHASE12D_ADDENDUM.md`).
- `docs/PROJECT_HISTORY.md` holds the full phase-by-phase narrative, the module layout table, and the assumptions /
  resolved / still-open lists (what used to be the whole of `README.md` before Phase 13 shortened it to a front
  door). `docs/ARCHITECTURE.md`, `docs/USER_GUIDE.md`, `docs/DEVELOPING.md` cover the Part 3 app specifically.

## Layout

- `src/` core modules (flat): `config`, `units`, `constants`, `materials`, `geometry`, `periodic_bc`, `excitation`,
  `green_function`, `basis_functions`, `mom_matrix`, `solver`, `sparameters`, `sweep`, `runner`, `database`,
  `polarization`, `phase`, library/feasibility modules (Phase 7/8A).
- `src/phase8b/ … src/phase12d/` — one subpackage per earlier phase, each with `run.py`, `gate.py`, `report.py`,
  `specification.py` etc. `tests/test_<phase>_<module>.py` mirrors them.
- `src/phase13/` — Part 0 carry-overs only (`convergence.py`, `resonance_restatement.py`,
  `hfss_resonance_ingest.py`, `order_power_prediction.py`, `part0_run.py`); analysis-only, reuses everything
  earlier unmodified. Parts 2/3's own new capabilities live in their own top-level packages, not under `phase13/`:
  `src/logging_setup.py`, `src/diagnostics.py`, `src/errors.py`, `src/labels.py`, `src/manifest.py`, `src/fields/`
  (HFSS field-calculator `.fld` import), `src/dispersion/` (surface-wave dispersion library), `src/hologram/`,
  `src/reflectarray/`, `src/pattern/`, `src/verify/`, `src/export/`, `src/hfss_compare/` (app-facing
  synthesis/prediction/verification/export services), `src/library_sweep.py` (generalizes `library_builder.py`'s
  sweep orchestration to an arbitrary, including grounded, config), `src/app/` (the FastAPI backend and static
  frontend — see `docs/ARCHITECTURE.md`).
- `examples/` — CLIs (`run_phaseN_*.py`) that write to `results/<phase_name>/`.
- `results/` — generated outputs (CSV/JSON/plots/reports). Treat existing library data as read-only benchmarks.
  Some large/regenerable subfolders and the real HFSS `.s36p`/`.s4p` exports are `.gitignore`d — see
  `docs/DEVELOPING.md`.
- `configs/` — config-file presets (`config.yaml`, `config_phase12c_grounded_8mm.yaml`); `docs/` — all non-root
  documentation (`theory/`, `phases/`, `briefs/`, plus `PROJECT_HISTORY.md`/`ARCHITECTURE.md`/`USER_GUIDE.md`/
  `DEVELOPING.md`).
- `tests/reference_layered.py` — an independent field-matching solver used only by tests.

## Code conventions

- Import as `from src.xxx import ...`; `pytest.ini` sets `pythonpath = .`, so run from the project root.
- Config is strict: no electrical value has a code default; unknown keys are errors; `mesh` and `reference_planes`
  sections are required. Mesh rule: `2·floquet_order + 1 ≥ max(cells_x, cells_y)`.
- New phases **reuse earlier phases unmodified** (wrap or duck-type, e.g. `SupercellUnitCell`); don't edit frozen
  earlier-phase modules to suit a new phase. Phase 12A, 12B, 12C and 12D are the deliberate exceptions: each is a core
  *capability* addition (an optional PEC ground plane; mesh alignment/duplicate-geometry detection; grounded
  heterogeneous supercells + supercell STEP export + one-port HFSS ingestion; a general N-port Touchstone reader), not
  an analysis-only phase, so each extends core modules directly (`green_function.py`/`geometry.py`/`config.py`/
  `excitation.py`/`sparameters.py`/`database.py`/`sweep.py`/`runner.py`/`cad_export.py`/`discretization.py`/
  `phase9/hfss_reference.py`) — always additively (a new optional field/parameter/function defaulting to, or beside,
  the pre-existing behaviour), verified bit-identical for every pre-existing config/sweep/result
  (PHASE12A_GROUND_PLANE.md, PHASE12B_MESH_PERIOD_STUDY.md, PHASE12C_GROUNDED_LOCALITY.md,
  PHASE12D_LOCAL_PERIODICITY.md). Phase 12C's own grounded-supercell wrapper
  (`phase12c.grounded_supercell.GroundedSupercellUnitCell`) does NOT touch Phase 9's frozen `SupercellUnitCell` —
  it wraps it, the same pattern every other phase already uses for reuse. Phase 12D's own aperture-level
  Floquet-order-power capability (`phase12d.floquet_power`) is new orchestration/analysis built entirely on
  already-existing, unmodified core functions (`basis_functions.RooftopBasis.current_coefficients`,
  `green_function.LayeredGreenFunction.sheet_fields`/`stack_fields`, `sparameters.decompose_air_fields`) — it does
  not itself edit any core module (only `phase9/hfss_reference.py`'s new `read_touchstone` does, additively).
- Outputs are deterministic (no timestamps in files); checkpoints are append-only and resumable.
- Tests check mathematical/numerical properties (reciprocity, energy conservation, symmetry, convergence),
  not just "it runs". Add tests with every change.

## Scientific-honesty rules (project-wide)

- Never fabricate data: no invented HFSS/measured results — report `NOT_AVAILABLE_FROM_CURRENT_DATA` instead.
- Gates return explicit named states (e.g. `CONTROL_VARIABLE_UNVALIDATED`, `DATABASE_COVERAGE_MISSING`,
  `NO_VALIDATED_FORMULATION`), never a bare pass/fail; don't upgrade a status without the evidence.
- Never choose S11 vs S21 (control variable) by phase coverage; never rescale/interpolate a library to a period,
  topology or substrate it wasn't solved for.
- Phase 11 (`src/phase11/`) is physics-only: no ML / fitted surrogates (a test greps the source for this).
- Report measured numbers, not assumed ones; don't fit decay laws or models before checking the data.

## Physical model assumptions

Rectangular lattice, one cell; air | zero-thickness PEC pattern at z = 0 | one isotropic substrate | air, or
(Phase 12A, opt-in via `substrate.backing: pec`) air | pattern | substrate | PEC ground at z = −h — no superstrate
either way; a grounded stack has no bottom port (S12/S21/S22 undefined). Normal incidence; only the (0,0) Floquet
order is a port mode; ports at the explicit reference planes in `config.yaml`. Only one library is actually built
(still ungrounded): 4.8 mm period, rectangular PEC patch, Rogers RT/Duroid 5880, 10–15 GHz. Phase 12C's own working
design point (grounded, P = 8 mm, same substrate) is recorded as a real, loadable config
(`config_phase12c_grounded_8mm.yaml`), not hard-coded — see `src/phase12c/design_point.py`.
