# Production repo + local web app (Phase 13)

Implementation: Part 0 carry-overs in `src/phase13/` (`convergence.py`, `resonance_restatement.py`,
`hfss_resonance_ingest.py`, `order_power_prediction.py`); Part 1 repo restructuring (root); Part 2
`src/logging_setup.py`/`errors.py`/`labels.py`/`manifest.py`/`diagnostics.py`; Part 3
`src/app/`, `src/dispersion/`, `src/fields/`, `src/hologram/`, `src/reflectarray/`, `src/pattern/`,
`src/verify/`, `src/export/`, `src/hfss_compare/`, `src/library_sweep.py`. Brief:
`docs/history/briefs/PHASE13_WEBAPP_PROMPT.md`; mid-session resume state:
`docs/history/briefs/PHASE13_HANDOFF.md`.

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
(`docs/ARCHITECTURE.md`). **Superseded by Part 3B/3C below**: the flat page list is now a 3-step
HMSA wizard plus a "More tools" page for the rest; this section (and the numbers below) describes
Part 3 as originally delivered and is kept as the historical record.

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
  the brief at the time — **since closed for the 3-step wizard's own flow, not for the "More
  tools" pages**, by Part 3B/3C below (`src/app/plots.py`).
- **Richardson mesh convergence is opt-in (`mesh_mode`), not the Unit Cell page's default** — the
  brief said "make it the app default for single-cell/library solves." A `richardson` solve costs
  several single solves (the `FAST_CELL_COUNTS` family), so defaulting to it would make the most
  common interactive action on the busiest page noticeably slower; `fast` (a single solve, clearly
  labelled `MESH_NOT_CONVERGED`) is the default instead, `richardson` one click away.
- **No run manifest / diagnostic-bundle download** (Part 2 built it, Part 3 never attached it to a
  job — see Part 2 above). Still true after Part 3B/3C.
- **Export does not accept an HMSA design**, only a Reflectarray one — **fixed in Part 3B (B2)
  below**; both now chain into Export the same way.
- **`UniversalDatabase` is not wired into the app's library registry** (a separate, lighter JSON
  registry is used instead — `docs/ARCHITECTURE.md` "Storage").
- **App-level end-to-end smoke tests exist as real `pytest` integration tests** (a small real
  library → reflectarray design → pattern prediction → export chain; a real `.fld` import → a
  real library → Formulation B synthesis), matching the brief's two required smoke tests, but they
  run through `fastapi.testclient.TestClient`, not a headless browser (none was available).

## Part 3B/3C — live-review fixes and the 3-step redesign

A live browser review of the running Part 3 app (`docs/history/briefs/PHASE13B_UI_REVIEW_FIXES.md`) found
several real bugs; a follow-up brief (`docs/history/briefs/PHASE13C_REDESIGN.md`) then asked for the app to
be rebuilt around one 3-step HMSA flow, superseding Part 3's flat 8-page layout and absorbing the
review's must-fix items. Both are implemented; commit history has the full per-item account.

**13B fixes (A1–A10, B2–B4), all landed:**
- **A1** Richardson mode now displays the extrapolated S11_xx (not the finest fixed-mesh solve),
  with an error bar defined as the extrapolation-correction size (`|v_inf - finest|`), not
  `curve_fit`'s covariance (not credible from a 4-point fit).
- **A2** Hologram synthesis resamples onto the CHOSEN period's cell-centre grid (complex
  interpolation) instead of running one "cell" per raw imported-field sample.
- **A3** An out-of-library-frequency request now raises `LIBRARY_FREQUENCY_UNAVAILABLE` naming
  what IS available; the upload's own frequency is persisted and compared against the synthesis
  frequency, labelling (`REFERENCE_FIELD_FREQUENCY_MISMATCH`) rather than silently diverging.
- **A4** "Reachable" now requires the phase error to be within a threshold
  (`reflectarray.designer.PHASE_ERROR_UNREACHABLE_THRESHOLD_DEG`, a documented policy default, not
  a measured limit), with mean/RMS/max error stats and a `PHASE_COVERAGE_SHORTFALL` label.
- **A5** An HMSA cell with no achieved phase now refuses (`phase_unreachable`) instead of a
  placeholder `1+0j` aperture entry.
- **A6** The Library page's grounded "Resonant?" column was always false: the classifier looked for
  a magnitude PEAK, but a low-loss grounded cell's own resonance is a phase crossing through
  ±180° at essentially constant |S11| (measured: 0.9957–0.9959 across a 1 GHz band at the project's
  own G4 design point). `library_sweep.detect_resonance_phase_crossing` fixes this for `backing=
  "pec"` cells only, leaving the ungrounded classifier untouched.
- **A7** TM0 check: could not reproduce the reported `tm0_comparison: null` via either
  `TestClient` or a live `uvicorn` process with the same inputs the review describes — the backend
  path was already correct when given all three of eps_r/thickness/frequency. Re-verified as part
  of the new Step 1 wizard flow instead of chased further as a standalone bug.
- **A8** Export's `n_cells` was actually N (the side length) — an N×N design reported "N cells".
  Split into `n_cells_per_side`/`total_cells` everywhere.
- **A9** `peak_memory_kb` read ~1000× too large on macOS: `ru_maxrss`'s unit is
  POSIX-unspecified (Linux: kB, macOS/BSD: bytes). `diagnostics._ru_maxrss_kb` normalizes per
  `sys.platform`.
- **A10** Verify page: the predicted intended-order fraction is now returned explicitly
  (`VerificationResult.predicted_intended_fraction`) instead of a broken ternary leaving it blank;
  the specular/intended order is shown alongside the spurious orders; "N cells (aperture side
  length)" reworded to "N (cells per 2π gradient period)".
- **B2** HMSA Formulation B results carry a reflectarray-shaped `assignments` list (with both L_mm
  and W_mm from the library match, previously only L was kept), chaining into Export exactly like a
  reflectarray design.
- **B3** Reference-wave rescaling to β̄ (the guided wavenumber implied by the hologram's own mean
  surface reactance) is on by default for Formulation A, applying an additive phase correction that
  preserves the imported field's own measured detail; opting out reports the accumulated-error
  estimate instead of staying silent.
- **B4** Formulation A now runs on a phase-only field with an explicit, labelled uniform-amplitude
  assumption (every real reference-field export this project has is phase-only) instead of
  refusing outright; Formulation B now also carries `EXPERIMENTAL` (a reflectarray-style control
  variable standing in for a genuinely different, surface-wave one).

**A deeper bug found while building 13C's own end-to-end test, not on the 13B list:** Formulation
B's aperture field for Pattern prediction was built as `exp(j*achieved_modulation_phase)` alone.
The actual radiated field of a modulated hologram is `E_ref(x,y) * Modulation(x,y)` — designed so
the reference wave's own phase cancels back out, leaving the pure desired-beam phase. Omitting that
multiplication silently propagated the reference wave's own phase structure into the predicted
aperture instead. Every existing hologram test happened to use a uniform-zero reference phase, the
one case where the bug is invisible (multiplying by 1 changes nothing) — an idealized
(zero-quantization) case with the real example `.fld` file and a 20°/0° object beam predicted a
peak at 51°/130° with the old construction, 21°/0° with the fix. `hologram.hmsa_synthesis.
synthesize_formulation_b` now computes and returns the correct field on `HmsaDesign.aperture_field`.

**13C redesign:** the app's home page is now a 3-step wizard (Reference wave → Beams & board →
Pattern & export), replacing the flat page list; everything else (Unit cell, Library, Reflectarray,
Pattern, Export, Verify, HFSS compare) moved to a "More tools" page, one click away, unchanged
otherwise. New pieces the wizard needed that didn't exist in Part 3:
- **Server-rendered PNG plots** (`src/app/plots.py`, matplotlib — closing Part 3's largest
  documented scope gap, at least for the wizard's own flow): the Step 1 reference phase map
  (inlined as base64 — no job id exists yet for a synchronous upload), a Step 3 layout/hologram
  preview, a pattern u-v map, and one principal cut per requested beam. "More tools" pages still
  show tables, not plots — not revisited.
- **`POST /api/hmsa/ensure_library`**: the wizard never sends the user to a separate Library page.
  Builds (or reuses, both across requests via the job cache and across app restarts via the
  persisted library registry) a small square-patch (L=W only) library for the chosen board and
  frequency automatically. Used a fixed 6-point L family and a fixed 32-cell mesh in this Part —
  not an alignment search (`discretization.smallest_aligned_cells` could need an unboundedly large
  mesh for an arbitrary period/L combination) and not a resonance search. **Known consequence,
  measured while building the walkthrough test**: for a board/frequency combination whose
  resonance falls between two of those 6 points, phase coverage can be poor (the app's own
  `PHASE_COVERAGE_SHORTFALL` label fires) — e.g. an 8 mm period at 10.54 GHz showed mean/RMS/max
  phase error 33°/41°/82° with 686 of 1369 cells unreachable in one run. This is reported honestly
  (the label is the point), not hidden. **Part 3D (below) made this denser and mesh-aligned by
  construction** — the underlying limited *range* of the phase-lookup method itself (not the
  library's resolution within that range) is what Phase 14 addresses with a genuinely different
  control variable.
- **`POST /api/export/package`**: bundles a completed export job's STEP/layout/README, an optional
  pattern job's plots, and a `project_record.json` (the export/pattern RESULTS, not the full
  original request bodies — the app does not currently persist those anywhere a later request
  could read back, so "reload later" means re-inspecting the numeric outcome, not resuming an
  edit session) into one downloadable zip.
- **A label→plain-English wording layer** (`app.js`'s `LABEL_DESCRIPTIONS`/`describeLabel`): the
  validation-label vocabulary stays machine-readable in the raw job result and in each label
  pill's tooltip; the visible pill text is a plain sentence, per 13C's own code→text table.
- **Apple-like light/dark theme** (`style.css`, `prefers-color-scheme` + a `data-theme` override
  hook), shared by every page so "More tools" kept its existing look, restyled.

**New deferrals/simplifications this introduced (beyond the ones already listed below):**
- **Aperture size is shown read-only**, pre-filled from the uploaded field's own extent — 13C asked
  for it "editable," but no cropping/resizing logic exists for the reference field or the synthesis
  grid; claiming it was freely editable would have been a UI that didn't do what it said.
- **No fast/accurate toggle for the auto-built library** (Part 3D added one — a finer-mesh tier,
  not Richardson extrapolation, which the auto-library builder still does not support for a
  multi-size sweep).
- **DXF export was never built** (`ezdxf` is not a project dependency, despite an old docstring
  claiming otherwise, now corrected) — STEP + CSV/JSON layout only. (Part 3D added a real,
  hand-written DXF export, no new dependency.)
- **Screenshots of each step (light/dark) were not captured** — 13C's "Tests & done-when" asks for
  them saved to `docs/screenshots/`; no browser or screenshot tool was available in this session
  either (same constraint Part 3 itself recorded — see "Verified against the app" above). The
  `curl`/`urllib`-driven live-server methodology substituted for it, as it did in Part 3, and caught
  the aperture-field bug above. (Part 3D's own brief was itself written from a real human
  walk-through with real screenshots, now checked into `docs/screenshots/` — see Part 3D below for
  what the fixes changed and what's still unverified visually.)

## Part 3D — live walk-through findings, fixed (`docs/history/briefs/PHASE13D_POLISH_AND_PHYSICS.md`)

A real human walk-through (light mode, 800 px pane; the exact case used throughout this section:
monopole field, 151x151, 2 mm grid, β=224.96 rad/m, 10.54 GHz; RT/duroid 5880, 1.575 mm; 1 beam
θ=20°, φ=0°; 4 mm cells) found real screenshots (`docs/screenshots/`, checked into git) and a
punch list of physics bugs and wording problems. Fixed, in the order the brief asked for:

**1. The MoM surface-wave dispersion capability is real now, not a scaffold.**
`src/dispersion/mom_scaffold.py` (kept its filename; the module itself is no longer a scaffold)
finds the periodic grounded cell's real TM surface-wave pole by scanning the largest singular
value of the assembled Z(kt) across `sparameters.surface_wave_interval` — a sharp, well-separated
peak (the bare-slab Green's-function pole enters additively through the (0,0) Floquet harmonic, so
it drives `||Z||` to a genuine spike, not the smallest-singular-value dip the original scaffold's
plan expected and this session's own experimentation found unreliable). Validated against
`analytic_tm0.tm0_beta`: **224.9792 rad/m found vs. 224.9797 analytic (0.0005 rad/m off) for a
lightly-loaded cell**, `tests/test_dispersion_mom.py`. Derives an equivalent surface reactance from
the found β via the transverse-resonance circuit (`Y_s = -Y_TM,bare(β_loaded)`) and builds a
`DispersionTable` in the same shape `hfss_table.py` produces, so Formulation A can use either
source interchangeably. Still labelled `EXPERIMENTAL` (young, self-validated on one board) —
**Phase 14 replaces this module with a purpose-built eigenmode solver** (Brillouin-zone-aware,
isotropy-checked, mesh-convergence-checked) per its own brief; this module was Part 3D's
scoped-down, "prove the physics works" version, not the final one.

**2. Pattern side-lobe level was wrong for two independent reasons.** The bare aperture radiation
integral (`phase8e.beam_reconstruction.far_field_pattern`) has no element/obliquity factor by
design (its own docstring: "not new radiation-integral physics"); for a tangential aperture field
on an infinite PEC ground plane (this stack) the standard result carries a `cos(theta)` factor —
without it, power didn't fall off toward grazing incidence, showing an unphysical bright ring near
the horizon. Separately, `detect_local_maxima` didn't wrap the periodic phi axis at its array
boundary, so a beam sitting exactly on the phi=0 grid line had its phi=(360-step) neighbour read as
`-inf` instead of the true (very high) main-lobe value one grid step away — a false "sidelobe" one
grid step from the real peak. Together: the walk-through's own case went from a reported SLL of
**-1.5 dB (wrong) to -10.9 dB**, matching the true first sidelobe visible on the φ=0 cut.
`tests/test_pattern_aperture_pattern.py` has synthetic regression tests for each cause separately;
`tests/test_app_wizard_walkthrough.py` asserts `sll_db < -5.0` on the walk-through's own case.

**3. STEP export is opt-in, not automatic, and reports estimated progress.** It was always built
as part of the same synthesis chain the pattern prediction ran in — for the walk-through's 5476
cells that took 126 s and produced a 103 MB file, with the UI showing "Done." (from the pattern
step) before it actually finished. Now: `export.array_export.write_layout_dxf` (new, see item 4)
and CSV/JSON/README build automatically and fast; STEP is a separate "Build STEP (CAD)" action with
its own status line. While it runs, `app/routers/export.py`'s `_StepProgressEstimator` ticks a
percentage from Phase 12D's own measured runtimes (N=6/8/12 ≈ 7/18.5/124 s), log-log interpolated —
labelled "estimated" everywhere it's shown, never claimed as measured.

**4. Zip contents.** `export.array_export.write_layout_dxf` is a real, hand-written DXF R12 file
(`LINE` entities only, valid since R12; no new dependency — `ezdxf` still isn't one) of the patch
layer, always generated alongside CSV/JSON/README. Separately, `export/package`'s zip silently
never got the layout PNG: `_run_package_job` read `layout_png_url` from the *export* job's own
result, which never had that key (the layout preview is generated by the *synthesis* job) — fixed
by adding an optional `synthesize_job_id` to `ExportPackageRequest`. `tests/
test_app_wizard_walkthrough.py`'s zip-contents assertion now requires `layout.dxf` and
`layout.png` alongside what it already checked.

**5. Auto-library: denser and mesh-aligned by construction.** The old 6-point L family (fractions
0.15, 0.30, ..., 0.90 of the period) was not a multiple of the mesh pitch (`2 / cells`), so it
carried staircasing error on top of being sparse. The new family (13 points, step 0.0625) is a
multiple of `2 / AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[accuracy]` by construction — no
`smallest_aligned_cells` search needed. `ensure_library` also gained an `accuracy: "fast" |
"accurate"` field (a finer 48-cell mesh, not Richardson extrapolation, which the multi-size
sweep still doesn't support) and now reports `phase_span_deg`/`resolution_deg` so Step 2 can show
the achievable range before the user commits to a synthesis run. **This did not fix the walk-
through's underlying complaint**: the phase-lookup method's achievable range is a property of the
control variable itself (reflection phase at small electrical size), not the library's sampling
density — the walk-through case still shows ~36° of range and ~74% of cells unreachable with 4 mm
cells on this board, same order of magnitude as before. Phase 14's eigenmode/impedance method is
the actual fix; Part 3D's job was to stop hiding the range problem and stop adding staircasing
error on top of it, not to solve it.

**Wording** (`src/app/static/`, all pages): brand name "HMSA Designer" in every page's title/header
(the Python package stays `lwmom`); every " -- " in visible UI text became "; " or an em dash;
bracket-style units (`[GHz]`, `[deg]`) became bare suffixes or "(°)"; plot titles/axes in
`app/plots.py` now read "φ = 0° cut" / "θ (°)"; the magnitude and dispersion file inputs use the
same styled button (`.file-picker` in `style.css`) as the phase dropzone instead of the browser's
native grey button; the Step 1 TM0 line now names the actual substrate (`substrateDisplayName()`)
and re-runs silently whenever Step 2's board changes; Step 2 hides the impedance-method-only
fields (mean reactance/modulation depth) when no dispersion table is uploaded, suggests a cell size
from the board's own measured wavenumber (`cells_per_fringe`'s relation run in reverse), and shows
the achievable phase-lookup range as soon as board/cell size/frequency are set; Step 3 shows a
one-line design summary before running, is titled "Design antenna," explains *why* when more than
20% of cells can't reach their required value (a plain sentence naming the cause and the fix, not
just the `PHASE_COVERAGE_SHORTFALL` label), and gives the "Experimental." chip a specific reason
in the phase-lookup-method context; the stepper has real top margin now.

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

**Added in Part 3B/3C** (`src/labels.py`): `reference_wave_rescaled_to_beta_bar(beta_bar)` /
`reference_wave_not_rescaled(worst_case_error_deg)` (B3), `reference_field_frequency_mismatch
(upload_hz, synthesis_hz)` (A3), and a redefinition of `phase_coverage_shortfall` (A4 — now the
worst per-cell phase error among unreachable cells, not an unused achieved-vs-target-span
metric). The wizard's wording layer (`app.js`'s `LABEL_DESCRIPTIONS`) maps every label above,
including these, to the plain-English text 13C's own table asks for; unrecognised codes fall back
to showing the raw code rather than hiding it.
