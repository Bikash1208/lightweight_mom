# Phase 13C brief — redesign the app around one 3-step HMSA flow (supersedes the ordering of 13B)

Read `docs/history/briefs/PHASE13B_UI_REVIEW_FIXES.md` first — it lists bugs found in a live review. This brief keeps the
bugs that affect the main flow as **must-fix**, and changes what the app is: the user wants **one simple, polished
flow**, everything else tucked away. Physics/numerics unchanged. Commit per step.

## 1. Information architecture
Top level has exactly three steps (a stepper: 1 · 2 · 3), plus one quiet link to **More tools**:

**Step 1 — Reference wave**
- Upload the monopole reference-field file (HFSS `.fld`, phase; optional magnitude file). Drag-and-drop area.
- Ask for the frequency it was exported at (the file doesn't contain it).
- Show immediately: a PNG of the reference phase map, and one plain-language summary line, e.g.
  "Cylindrical surface wave from the centre · β = 225 rad/m · matches the bare substrate at 10.54 GHz ✓"
  (TM0 check — currently broken, 13B A7). Mismatch → friendly warning, not a code.

**Step 2 — Beams & board**
- *Beams*: number of beams (1–4), then a small table per beam: θ (°), φ (°), relative power (dB, default 0).
- *Board*: substrate (preset dropdown — Rogers RT/duroid 5880, RO4003C, FR-4… — plus Custom: εr, loss tangent,
  thickness), cell size (mm), aperture size (pre-filled from the reference file's extent, editable).
- *Frequency*: pre-filled from Step 1, locked unless the user unlocks it (mismatch must be explicit — 13B A3).
- **Advanced** (collapsed): method (impedance-modulated / phase-lookup), modulation depth, mean reactance, accuracy
  (fast / accurate = mesh-extrapolated), upload an HFSS dispersion table, reference-wave correction to the loaded
  wavenumber (on by default — 13B B3).
- The unit-cell data the synthesis needs is built automatically for the chosen board and frequency (square patches,
  aligned mesh, cached by config) with a progress bar and time estimate. The user never has to visit a "Library"
  page. Be honest about the method: the impedance method needs dispersion data (HFSS table today; MoM dispersion is
  not implemented) — if none is available, fall back to the phase-lookup method and say so in one plain sentence.

**Step 3 — Pattern & export**
- Run synthesis (resampled onto the cell grid — 13B A2), then show **PNG plots**: radiation pattern (2-D u–v map in dB)
  and principal cuts through each beam; the hologram/layout preview.
- Key results as a short list: achieved vs requested beam directions, directivity, side-lobe level, cells used,
  how many cells fall outside what the unit cell can realize (with real phase/impedance error stats — 13B A4).
- **Export**: one "Download package" (zip) + individual items: board CAD (STEP: substrate + patches + ground), layout
  (CSV and DXF), pattern plots (PNG), HFSS setup notes for the **full antenna** (not unit-cell Floquet advice —
  13B B10), and the project record (inputs + results, to reload later).

**More tools** (one page with a list, or a drawer): Unit cell, Unit-cell library, Reflectarray designer,
Full-wave check (Verify), HFSS comparison, Diagnostics & logs. These keep working; 13B fixes for them are lower
priority (do A1, A6, A8–A10 there; defer the rest in `docs/history/PHASE13_WEBAPP.md`).

## 2. Look and feel — Apple-like
- System font stack (`-apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", Inter, sans-serif`),
  light theme by default with automatic dark mode (`prefers-color-scheme`), off-white background, white cards with
  large radius (12–16 px) and soft shadows, generous spacing, one accent colour (system blue), one primary button
  per step, secondary actions as plain text buttons, segmented controls instead of dropdowns where there are ≤ 4
  choices, subtle transitions. No framework or build step needed — plain CSS is fine. Must work offline and at
  800 px width (the review ran in an 800 px pane).
- Plots (matplotlib → PNG) styled to match: clean, light grid, same accent colour, readable at card width, dark
  variants in dark mode.

## 3. Words — short, human, production-grade
- Field labels in plain English with units as suffixes: "Frequency" [GHz], "Cell size" [mm], "Beam 1 — θ" [°].
  No snake_case, no internal names, no phase numbers, no job IDs on screen (keep them in the project record).
- Replace validation-label codes on screen with one-line notes; keep the codes in the manifest/diagnostics:
  | code | on screen |
  |---|---|
  | LOCAL_PERIODICITY_APPROXIMATION | "Estimated from individual cells — confirm in HFSS before fabrication." |
  | CELLS_PER_FRINGE(3.5) | "Cells are large for this wave (3.5 per fringe). Try ≤ 4 mm for a cleaner pattern." |
  | FEED_REGION_MASKED | "Feed point excluded." |
  | PHASE_COVERAGE_SHORTFALL(x°) / unreachable | "N cells can't reach the needed value (worst error x°)." |
  | MESH_NOT_CONVERGED | "Fast accuracy — switch to Accurate in Advanced for final designs." |
  | REFERENCE_AMPLITUDE_ASSUMED_UNIFORM | "No magnitude file — uniform amplitude assumed." |
  | EXPERIMENTAL | "Experimental." |
- Errors: what happened + what to do, one sentence each ("This file has no frequency information. Enter the
  frequency it was exported at.").

## 4. Must-fix before calling this done (from 13B)
A1 (richardson display), A2 (resample to cells), A3 (frequency mismatch), A4 (real error stats + coverage), A5 (no
placeholder aperture), A7 (TM0 check), B2 (HMSA export), B3 (β̄ correction), B4 (method labelling; impedance method
runs phase-only with uniform amplitude).

## 5. Tests & done-when
- API-level walk-through test of the 3-step flow with the real example `inputs/reference_fields/efl_phase_EzZ_example.fld`
  (frequency 10.54 GHz, 1 beam at θ=20°, φ=0, RT/duroid 5880 1.575 mm, cell 4 mm) ending in a valid zip whose STEP
  has one patch per cell plus substrate and ground, and whose pattern peak lies within a stated tolerance of 20°.
- A two-beam case (θ=20°, φ=0 and θ=20°, φ=180) — both beams present in the pattern.
- No snake_case or internal code visible in any rendered page (a test that scans the served HTML/JS strings).
- Screenshots of each step (light and dark) saved to `docs/screenshots/` for review; fast suite green; commit.
