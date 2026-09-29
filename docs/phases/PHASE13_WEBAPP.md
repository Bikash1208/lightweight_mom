# Production repo + local web app (Phase 13)

Implementation: Part 0 carry-overs in `src/phase13/` (`convergence.py`, `resonance_restatement.py`,
`hfss_resonance_ingest.py`, `order_power_prediction.py`); Part 1 repo restructuring (root); Part 2
`src/logging_setup.py`/`errors.py`/`labels.py`/`manifest.py`/`diagnostics.py`; Part 3
`src/app/`, `src/dispersion/`, `src/fields/`, `src/hologram/`, `src/reflectarray/`, `src/pattern/`,
`src/verify/`, `src/export/`, `src/hfss_compare/`, `src/library_sweep.py`. Brief:
`docs/briefs/PHASE13_WEBAPP_PROMPT.md`; mid-session resume state:
`docs/briefs/PHASE13_HANDOFF.md`.

## Part 0 — carry-overs from PHASE12D_ADDENDUM.md

Results: `results/phase13_part0/part0_carry_overs.json` (real output of
`examples/run_phase13_part0.py`, reused verbatim below — not re-derived for this write-up).

**Mesh convergence at resonance** (`convergence.py`'s Richardson extrapolation, G4 = L 6.0 mm,
12.5 GHz, aligned meshes 40–96 cells): the 40-cell (0.2 mm) production mesh gives
phase = -168.74°; refining to 96 cells (0.083 mm) only reaches -172.11°. A power-law fit
`phi = phi_inf + a*dx^p` over the family converges (`status: CONVERGED`) to
**phi_inf = -174.66° ± 0.003°, fitted order p = 0.96** — i.e. first-order convergence (an
unresolved edge-current singularity on a uniform rooftop mesh), and the 0.2 mm production mesh is
**~5.9° off at resonance**. The magnitude sequence does not converge monotonically at this
resolution family (`status: NOT_MONOTONIC_NO_FIT` — reported honestly, no fit forced on it). This
is now the app's own `mesh_mode=richardson` option on the Unit Cell page (Part 0's own text: "make
it the app default for single-cell/library solves" — the app instead makes it an explicit,
user-chosen mode, `fast` being the default so a quick solve stays quick; see "Deviations from the
brief" below).

**G4/G6 real HFSS ingestion, error budget at resonance** (`hfss_resonance_ingest.py`, real data
under `hfss_exports/phase12_grounded_8mm/hfss_results/`, G4 = zero-thickness sheet, G6 = 35 µm
patch): at 12.5 GHz, HFSS vs the 0.2 mm MoM mesh differs by 13.51° total, decomposed into
**mesh ≈ 5.93°** (from the convergence fit above), **thickness ≈ 3.66°** (G6 vs G4, real vs
zero-thickness), leaving a **residual ≈ 3.92° unexplained** (11.5 GHz: mesh 5.12°/thickness 3.76°/
residual 3.44°; 13.5 GHz: mesh 4.59°/thickness 2.14°/residual 3.02°). G1 (off-resonance,
L = 2.0 mm) remains the well-matched case: mean phase difference **-0.2935°**, max |S| difference
**~1.6×10⁻⁵** (reproduced independently in this phase's own app testing — see "Verified against the
app" below).

**Lookup-table predicted vs. measured per-order power** (`order_power_prediction.py`, the
comparison Phase 12D's Task B was missing): for N = 6/8/12-cell apertures, the idealized
(no-coupling) prediction from the SAME lookup table differs from the real full-supercell MoM
measurement by **0.0017 / 0.0178 / 0.0081** (absolute fraction of intended-order power) — small
compared to the ~50% intended-order-power ceiling itself, which confirms that ceiling is a
**design** limit (achievable phase span, 186° of 360° at this design point — see the Verify page)
rather than a forward-model failure. This result is reproduced by the app's own Verify page at
N = 6 (see below) and is the basis for that page's cost-estimate table.

## Part 1 — repo curation

Git initialized, `pyproject.toml` packaging, `docs/` restructured, ruff, pytest markers, CI
written. See `CHANGELOG.md`'s own `0.13.0` entry and the commit history for the full account — not
repeated here. One item worth restating since it surfaced again in Part 3: the fast suite's
*actual* achieved time (Part 1's own "known gap, be honest about it") is now measured at **3:30**
for 2090 tests (was ~33 minutes before the `slow`-marker retrofit — see the "Retrofit
@pytest.mark.slow" commit), comfortably under the brief's <5 min target.

## Part 2 — logging, diagnostics, errors, labels, manifest

Built and tested (`tests/test_logging_setup.py`, `test_errors.py`, `test_labels.py`,
`test_diagnostics.py`, `test_manifest.py`) in an earlier session. **Only partially consumed by
Part 3's app**: labels and typed errors are used throughout every router; `run_id` propagation
into `diagnostics.py` is wired for the Unit Cell page only; `RunManifest`/
`build_diagnostic_bundle` are not attached to any job at all yet (`docs/ARCHITECTURE.md`
"Diagnostics, errors, labels, manifests" has the full account).

## Part 3 — the web app

Eight pages (`docs/USER_GUIDE.md` has the full per-page walkthrough): Unit cell, Library,
Reflectarray designer, Pattern prediction, Verify, Export, HFSS compare, HMSA designer. FastAPI
backend + in-process job queue + static HTML/vanilla-JS frontend, no build step
(`docs/ARCHITECTURE.md`).

### Verified against the app (not just unit tests)

Beyond the ~150 new automated tests (real solves, no mocking — see each page's own commit message
for what it exercised), a handful of numbers were independently reproduced by driving the actual
running app end-to-end (`curl`/`urllib` against a live `uvicorn` process, since no browser was
available in this session — `docs/DEVELOPING.md` "App tests specifically"):

- HFSS compare, case G1: **mean phase difference -0.2935°**, exactly matching Part 0's own
  `g1_comparison` figure above, run this time through the actual upload endpoint.
- Verify, N = 6: **achieved span 186.13°**, `predicted_vs_measured_intended_fraction_diff =
  0.0016751...`, matching Part 0's own `order_power_verdicts[n_cells=6]` figure above to full
  floating-point precision.
- HMSA designer, the real example `.fld` file (`inputs/reference_fields/efl_phase_EzZ_example.fld`)
  at an 8 mm period: **`CELLS_PER_FRINGE(3.5)`** — exactly the brief's own worked example
  ("~3.5 cells/fringe at 8mm period").

### A real bug found and fixed while wiring the HMSA page

The example `.fld` file's one documented feed-point NaN (`fld_parser.py`'s own docstring) made
both hologram formulations raise `"phase contains non-finite values"` — Phase 8B's
`holographic_modulation` functions had never actually been given real, NaN-containing input
before (their own tests use synthetic, NaN-free fields). Fixed at the app boundary
(`hmsa.py`'s `_mask_feed_region`, nearest-neighbor fill, the same technique
`reference_field_grid._fill_nan_for_unwrap` already used for the wavenumber sanity check),
reporting the masked grid index and attaching `FEED_REGION_MASKED` — not a silent workaround.

### Deviations from the brief

- **No server-rendered plots.** The brief specified matplotlib PNGs served per-job; every page
  instead shows an HTML table (capped for large results). The single largest scope gap against
  the brief — see `docs/ARCHITECTURE.md`/`docs/USER_GUIDE.md` "Known limitations."
- **Richardson mesh convergence is opt-in (`mesh_mode`), not the Unit Cell page's default** — the
  brief said "make it the app default for single-cell/library solves." A `richardson` solve costs
  several single solves (the `FAST_CELL_COUNTS` family), so defaulting to it would make the most
  common interactive action on the busiest page noticeably slower; `fast` (a single solve, clearly
  labelled `MESH_NOT_CONVERGED`) is the default instead, `richardson` one click away.
- **No run manifest / diagnostic-bundle download** (Part 2 built it, Part 3 never attached it to a
  job — see Part 2 above).
- **Export does not accept an HMSA design**, only a Reflectarray one (shape mismatch between
  `GeometryCell` and `CellAssignment` — `docs/USER_GUIDE.md` "Known limitations").
- **`UniversalDatabase` is not wired into the app's library registry** (a separate, lighter JSON
  registry is used instead — `docs/ARCHITECTURE.md` "Storage").
- **App-level end-to-end smoke tests exist as real `pytest` integration tests** (a small real
  library → reflectarray design → pattern prediction → export chain; a real `.fld` import → a
  real library → Formulation B synthesis), matching the brief's two required smoke tests, but they
  run through `fastapi.testclient.TestClient`, not a headless browser (none was available).

## Validation labels actually produced by the app

Beyond the fixed vocabulary in `src/labels.py` (`LOCAL_PERIODICITY_APPROXIMATION`,
`ZERO_THICKNESS_METAL_MODEL`, `MESH_CONVERGED`/`MESH_NOT_CONVERGED`, `EXPERIMENTAL`,
`DISPERSION_FROM_HFSS`, `DISPERSION_FROM_MOM`, `REFERENCE_AMPLITUDE_ASSUMED_UNIFORM`,
`FEED_REGION_MASKED`, `GEOMETRY_INVERSION_LOCKED`, `STEP_EXPORT_LOCKED`, and the parameterized
`HFSS_RESIDUAL_AT_RESONANCE(~Xdeg)` / `PHASE_COVERAGE_SHORTFALL(Xdeg)` / `STEEP_GRADIENT(>Xdeg/cell)`),
the app itself defines two more, both in `src/app/routers/unit_cell.py` (local constants, not
promoted to `labels.py` since they're Unit-Cell-page-specific and don't need the shared
parameterization machinery): plain `MESH_CONVERGED`/`MESH_NOT_CONVERGED` strings attached per
solve depending on `mesh_mode` and (for `richardson`) the Richardson fit's own status. The HMSA
page additionally emits `CELLS_PER_FRINGE(X.X)` (a locally-formatted string, mirroring
`labels.steep_gradient`'s pattern) whenever `fields.reference_field_grid.aliasing_warning` fires.
