# Phase 13 handoff — read this first if resuming in a new session

This session hit a hard block partway through Phase 13: a safety classifier started refusing
every code-execution action (pytest runs, `python -c`, even `ps aux`), explicitly stating it was
unrelated to the actual command content and would keep firing for the rest of that conversation.
The user chose "start a fresh session to continue." This file is the exact state to resume from.

**Read `PHASE13_WEBAPP_PROMPT.md` (in this same `docs/briefs/` directory) for the original brief.**
Everything below assumes it.

## What's done and fully test-verified (safe to build on without re-checking)

### Part 0 — carry-overs from PHASE12D_ADDENDUM.md
All done. `src/phase13/`: `convergence.py` (Richardson mesh-convergence extrapolation, reproduces
the addendum's G4 numbers: phi_inf=-174.66deg, order=0.96), `resonance_restatement.py` (scales that
fit to Phase 12C/12D's own meshes), `hfss_resonance_ingest.py` (G4/G6 real HFSS ingestion + error
budget: mesh~5.9deg, thickness~3.7deg, residual~3-4deg), `order_power_prediction.py` (the missing
lookup-table-predicted-vs-measured comparison: N=6/8/12 diffs of 0.002/0.018/0.008 — confirms the
~50% intended-order-power ceiling is a DESIGN limit, not a lookup-table forward-model failure),
`part0_run.py` + `scripts/research/run_phase13_part0.py`. **Already run for real**; output at
`results/phase13_part0/part0_carry_overs.json`. All tests pass (`tests/test_phase13_*.py`).

### Part 1 — repo curation
All done. git initialized (local commits only — `git config user.name "biksmac"` /
`user.email "bikashkumarshanti@gmail.com"`, set in this session since none existed; the user should
confirm/change this before any real push). `.gitignore` excludes `.venv/`, `__pycache__/`,
`.pytest_cache/`, `results/unit_cell_library/{_solve,database,geometry}/`, `*.s36p`/`*.s4p`,
`hfss_exports/_src_snapshot.tgz`. `pyproject.toml` (distribution name `lwmom`, import path stays
`src.*` — deliberately NOT renamed, justified in the file's own comment: 200+ files import
`src.*`, a mechanical rename was judged not worth the risk/time here; revisit if the user disagrees).
Docs moved: `docs/theory/` (CONVENTIONS/GREEN_FUNCTION/MOM_FORMULATION/S_PARAMETERS/DATABASE),
`docs/history/` (every PHASE*.md write-up + UNIT_CELL_LIBRARY/UNIVERSAL_UNIT_CELL_DATABASE/
VALIDATION_EXPORT), `docs/briefs/` (every `*_PROMPT.md` + `PHASE12D_ADDENDUM.md` + this file).
`docs/history/PROJECT_HISTORY.md` is the OLD full README.md content (all internal links fixed to the new
paths); root `README.md` is now short, per the brief. `configs/` holds `config.yaml` and
`config_phase12c_grounded_8mm.yaml` (moved from root; the 2 real code references —
`src/library_schema.py`, `src/phase12c/design_point.py` — and their tests were updated and
re-verified). `ruff` configured (pragmatic ignores: E501/E702, documented why) and run once
(`ruff check --fix` + `ruff format`) across the whole tree as its own pass. `pytest.ini` has
`slow`/`cad`/`checkpoint`/`hfss_data` markers. `.github/workflows/ci.yml` written (not run).
`Makefile`, `requirements-lock.txt` (pip freeze of the working `.venv`), `CHANGELOG.md`,
`CONTRIBUTING.md` all written.

**Known gap, be honest about it**: the `slow` marker is only applied to Phase 13's OWN new
expensive tests (~8 of them). It has NOT been retrofitted onto the many pre-existing real-MoM-solve
tests scattered across Phases 9-12D's test files, so `pytest -m "not slow"` right now excludes
almost nothing and still takes close to the full ~32 minutes the complete suite takes (a background
run of it in this session was still going after 25+ minutes when the block hit, never confirmed
finished). **The brief's "<5 min fast suite" target is NOT yet met.** To finish this: run
`pytest -m "not slow" --durations=40 -q`, take the worst ~20-40 offenders (almost certainly
supercell/library-scale solves in `tests/test_phase9_*`, `test_phase10_*`, `test_phase11_*`,
`test_phase12c_*`, `test_phase12d_*`), and add `@pytest.mark.slow` (or a module-level
`pytestmark = pytest.mark.slow`) to them. Report the ACTUAL achieved fast-suite time honestly, do
not claim <5 min without measuring it.

### Part 2 — logging/diagnostics/errors/labels/manifest
All done and tested. `src/logging_setup.py` (LWMOM_LOG_LEVEL, rotating JSON file log + console,
contextvar-based run-ID propagation), `src/errors.py` (typed `LwmomError` + 10 stable error codes),
`src/labels.py` (the fixed validation-label vocabulary, plain and parameterized), `src/diagnostics.py`
(`instrumented_compute_s_parameters` — wraps `sparameters.compute_s_parameters` UNMODIFIED, adds
timing/condition-estimate/grating-margin/surface-wave-margin/energy-conservation/reciprocity
diagnostics from the outside, translates any failure into `errors.solve_failed`), `src/manifest.py`
(`RunManifest` with git hash/dirty flag, dependency versions, output-file hashing, and
`build_diagnostic_bundle` — the "download diagnostic bundle" zip). Tests:
`tests/test_logging_setup.py`, `test_errors.py`, `test_labels.py`, `test_diagnostics.py`,
`test_manifest.py` — all pass.

## Part 3 — the web app: what's written and test-verified so far

Physics/math pieces, all done and passing their own tests **before the block hit**:
- `src/dispersion/analytic_tm0.py` — TM0 grounded-slab dispersion, closed-form, reproduces the
  brief's own cited beta values (213/269/326 rad/m at 10/12.5/15 GHz for RT/Duroid 5880 1.575mm)
  almost exactly. `tests/test_dispersion_analytic_tm0.py`.
- `src/dispersion/hfss_table.py` — HFSS eigenmode dispersion CSV import (`DISPERSION_FROM_HFSS`),
  linear interpolation + inverse lookup. `tests/test_dispersion_hfss_table.py`.
- `src/dispersion/mom_scaffold.py` — the MoM eigenmode capability, deliberately NOT implemented
  (brief explicitly allows shipping this as a documented, tested scaffold that raises a specific
  error code). `tests/test_dispersion_mom_scaffold.py`.
- `src/fields/fld_parser.py` — the real HFSS field-calculator `.fld` format parser. Round-trips
  the user's actual example file (`inputs/reference_fields/efl_phase_EzZ_example.fld`, 151x151,
  detects x-major ordering, finds the single feed-point NaN). `tests/test_fld_parser.py` (16 tests).
- `src/fields/reference_field_grid.py` — builds a `phase8e.reference_field_input.ReferenceFieldGrid`
  from a parsed `.fld` file; grid-regularity check; 2-D phase unwrap + local-wavenumber estimate
  (empirically reproduces the brief's own cited ~-225 rad/m plateau over r=10-140mm on the REAL
  example file); TM0 comparison; cells-per-fringe/aliasing warning (reproduces the brief's own
  "~3.5 cells/fringe at 8mm period" worked example). `tests/test_reference_field_grid.py` (10 tests).
- `src/hologram/hmsa_synthesis.py` — both Phase 8B hologram formulations (A: Fong impedance-
  modulation, needs a dispersion table; B: phase-difference, needs only a normal-incidence-phase
  library) wired end to end: reference field -> object beam -> control map -> geometry, with
  `IMPEDANCE_UNREACHABLE`/`EXPERIMENTAL`/`LOCAL_PERIODICITY_APPROXIMATION` labels.
  `tests/test_hologram_hmsa_synthesis.py` (4 tests).
- `src/reflectarray/designer.py` — standard reflectarray required-phase-map formula (feed distance
  + steering term), `cos^q` illumination taper, cell selection via
  `phase8d_rev.response_matching.nearest_response_match` (reused unmodified), per-cell phase-
  gradient check against Phase 12D's own 5.8deg/cell limit (`STEEP_GRADIENT` label).
  `tests/test_reflectarray_designer.py` (6 tests).
- `src/pattern/aperture_pattern.py` — far-field/beam-pointing/SLL/directivity from an aperture
  field, built on `phase8e.beam_reconstruction.far_field_pattern` (reused unmodified). Verified
  against array theory (uniform aperture -> broadside; linear phase gradient -> correctly steered
  beam). `tests/test_pattern_aperture_pattern.py` (4 tests).
- `src/verify/supercell_check.py` — the "Verify" page: cost-estimate/refusal gate fit from Phase
  12D's own 3 measured runtimes (N=6/8/12: 7.0/18.5/124.2s), wraps
  `phase12d.aperture_test.run_aperture_test` + `phase13.order_power_prediction.predict_order_powers`
  unmodified. `tests/test_verify_supercell_check.py` — fast tests pass; the one `@pytest.mark.slow`
  end-to-end test (real N=6 solve, ~7s) was written but never actually run before the block hit —
  run it first thing in the new session.
- `src/export/array_export.py` — full-array STEP export (reuses
  `phase12c.grounded_supercell.build_grounded_supercell` + `cad_export.export_supercell_step`
  unmodified; square-grid-only, documented limitation), layout CSV/JSON, HFSS-ready README.
  Manually verified at a realistic 20x20 (400-cell) array: 6.98 MB STEP file, 5.26s export time —
  no streaming/tiling needed at this scale. `tests/test_export_array_export.py` (3 tests, all
  passed before the block).

**Written but never actually run/verified** (the block hit right after writing these — read them
carefully before trusting them):
- `src/hfss_compare/service.py` (upload validation + comparison + on-disk ledger, reusing
  `phase9.hfss_reference.read_touchstone`/`phase12d.hfss_ingest.compare_single_cell` unmodified)
  and its test file `tests/test_hfss_compare_service.py` (4 tests). This is thin glue over already-
  validated functions, so the risk is low, but it has literally never been executed — run it first.

## Part 3 — NOT STARTED at all

This is the majority of the remaining work:

1. **The FastAPI backend itself**: `src/app/main.py` (app factory + `main()` entry point that
   `pyproject.toml`'s `lwmom-app` script and `Makefile`'s `app` target already reference but which
   does not exist yet), a job queue (`src/app/jobs.py` — in-process, e.g. `ThreadPoolExecutor` +
   a dict keyed by job_id/config-hash, progress/cancel/cache-hit; "keep it simple" per the brief,
   this is a single local machine, no need for Celery/Redis), pydantic request/response schemas,
   and one FastAPI router per page wrapping the services above (most of the actual physics/service
   layer already exists above — the routers are mostly thin HTTP glue).
2. **The frontend**: plan was a single static HTML shell + vanilla JS (fetch calls to the FastAPI
   JSON endpoints) + server-side matplotlib-rendered PNG plots (served from a `/plot` endpoint per
   job) — deliberately NO client-side charting library and NO Node/build step, so the app stays
   fully offline with zero new frontend dependencies (this was the intended justification for
   `docs/ARCHITECTURE.md`'s "why no frontend build step" section — not yet written). Not started.
3. **The 8 pages**: unit cell, library, reflectarray designer, pattern prediction, verify, export,
   HFSS compare, HMSA designer. The service-layer logic for reflectarray/pattern/verify/export/
   hfss_compare/hologram already exists (above); "unit cell" and "library" pages need thin wrappers
   around `sparameters.compute_s_parameters`/`diagnostics.instrumented_compute_s_parameters` and
   the existing Phase 7 library query API (`unit_cell_library.py`) respectively — not yet written.
4. **App-level tests**: API test per endpoint, job lifecycle (queue/progress/cancel/cache-hit),
   error-code propagation to the UI, manifest completeness, and the two required end-to-end smoke
   tests (design + export a small reflectarray; import the example `.fld` + a synthetic dispersion
   table + build + export a hologram). None written.
5. **Docs**: `docs/USER_GUIDE.md`, `docs/ARCHITECTURE.md`, `docs/DEVELOPING.md`,
   `docs/history/PHASE13_WEBAPP.md` (what was built, the Part 0 numbers above, known limitations,
   validation-label definitions) — none written yet. `README.md`'s "Run the app" section already
   references `make app` / `lwmom-app` as if they exist; keep them working once `src/app/main.py`
   is real, or fix the README if the final command differs.
6. Final steps once the app exists: fresh-clone check (clone to a temp dir, `pip install -e .`,
   run the fast suite), retrofit the `slow` marker (Part 1's known gap, above), full regression
   suite run with final counts/runtime, and the final report to the user (test counts+runtimes,
   app start command, fresh-clone result, repo size, what's excluded from git and how to
   regenerate it, open questions: license — none exists, brief says ask the user; remote — none
   added, brief says don't add one; package name — kept `src.*`, distribution name `lwmom`,
   justified above).

## Practical notes for the new session

- `git log --oneline` to see what's committed so far; there may be uncommitted work-in-progress
  from this session's very end (the block hit mid-commit-attempt) — check `git status` first and
  commit or fold in anything pending before starting new work.
- A background test run (`pytest -m "not slow" -q`, checking Parts 0-2 + everything in Part 3 up to
  `array_export`) was started in this session and never confirmed finished when the block hit — its
  actual pass/fail result is unknown. Re-run it before trusting that nothing regressed.
- Every module written above follows the established project convention: reuse earlier
  phases/modules unmodified, typed errors from `src/errors.py`, labels from `src/labels.py`,
  docstrings explaining WHY not just WHAT. Keep matching that style for the remaining app code.
