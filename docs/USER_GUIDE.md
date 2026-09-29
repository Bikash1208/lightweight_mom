# User guide — the local web app

`make app` (or `lwmom-app`) starts the backend and opens `http://127.0.0.1:8000` in your browser.
Everything below runs locally; nothing leaves the machine. No screenshots are included here — no
browser was available in the session that wrote this guide to capture any (see
`docs/phases/PHASE13_WEBAPP.md` "Known limitations").

Every page follows the same pattern: fill in a form, submit, watch a status line poll
(`queued` → `running` → `done`/`failed`), then see a table of results plus a row of coloured
label pills (`LOCAL_PERIODICITY_APPROXIMATION`, `MESH_NOT_CONVERGED`, ...) — read these; they are
the app's honesty mechanism, not decoration. A `failed` job's status line shows the error message;
the underlying error also has a stable `code` (e.g. `CONFIG_INVALID`, `SOLVE_FAILED`,
`IMPEDANCE_UNREACHABLE`) if you're scripting against the API directly instead of the UI.

## 1. Unit cell

Solve one rectangular-patch unit cell (grounded, `backing: pec`, one port; or ungrounded, two
ports) at one frequency or a comma-separated list. **Mesh mode**:
- `fast` — solves once at the mesh you specify; always labelled `MESH_NOT_CONVERGED` since it was
  never checked against a finer mesh.
- `richardson` — runs a small family of mesh sizes and Richardson-extrapolates to the zero-cell-
  size limit (Part 0's own addition); labelled `MESH_CONVERGED` or `MESH_NOT_CONVERGED` depending
  on whether the fit actually converged (not every geometry/frequency does — a flat or
  non-monotonic sequence is reported honestly, not forced).

"Export STEP" builds the same unit cell's physical geometry (substrate + zero-thickness-modelled
patch, drawn at a configurable thickness for CAD purposes only) and downloads a `.step` file.
Needs the optional `cad` extra (`pip install -e '.[cad]'`) — fails with `CAD_UNAVAILABLE`
otherwise.

## 2. Library

**Browse**: pick a registered library (the read-only Phase 7 benchmark — 324 cells, 4.8 mm,
ungrounded, flagged with its own duplicate-geometry caveat — or one you've built) and a frequency;
see its phase-vs-geometry table and coverage (largest phase gap, span).

**Run a new sweep**: solves every (L, W) combination in the ranges you give, at the frequency grid
you give, for *any* period/substrate/backing — including grounded, which the original Phase 7
library never covered (`docs/DEVELOPING.md` "Why a second library-sweep module"). A small sweep
(a handful of L/W values, one frequency) takes seconds; a fine sweep over many frequencies can
take minutes, same cost as running that many unit-cell solves. The result is registered under the
library ID you give it and immediately available to the Reflectarray/HMSA pages.

## 3. Reflectarray designer

Pick a library, a frequency, an array size (cells × cells) and period, a feed position and beam
direction, and a feed-taper exponent `q`. Computes the required phase at each cell (standard
feed-distance-plus-steering formula), picks the nearest-response cell from the library at each
point, and reports the achieved-vs-required phase error and the largest per-cell phase gradient
against Phase 12D's own measured 5.8°/cell off-resonance design limit (`STEEP_GRADIENT` if
exceeded). The result's job ID links directly to the Pattern and Export pages with the design
pre-selected.

## 4. Pattern prediction

Paste in a design job's ID (from Reflectarray or HMSA) and a frequency; computes the far field
from that design's synthesized aperture (main beam pointing, side-lobe level, a directivity
estimate) via aperture-integration Fourier optics. Always labelled
`LOCAL_PERIODICITY_APPROXIMATION` — each cell's contribution comes from its own isolated-cell
library entry, not a full-array solve (see Verify below for that check).

## 5. Verify

Runs a *real* full-supercell MoM solve of an N×N periodic-gradient aperture and compares its
actual per-Floquet-order power against the idealized (no inter-cell coupling) prediction from the
same lookup table — the one place the app checks its own forward model against ground truth
rather than trusting the library. Costed before running (Phase 12D's own measured runtimes: N=6
≈7s, N=8 ≈18.5s, N=12 ≈124s); sizes beyond N=20 are refused outright as untrusted extrapolation,
sizes between the largest measured point and that limit require an explicit "run anyway" click.

## 6. Export

Chains off a completed Reflectarray design (**not** an HMSA design — see "Known limitations"
below) to produce a full-array STEP file, layout CSV/JSON, and an HFSS-ready README (Floquet-port
mode count, PEC-vs-copper, zero-thickness-vs-drawn-thickness lessons from Phase 12). Grounded
designs only, and the array must be square (N×N) — the underlying supercell builder's own
limitation, inherited unmodified.

## 7. HFSS compare

Upload a Touchstone (`.sNp`) export and pick the matching single-cell case (G1–G6); overlays it
against this project's own MoM reference and reports mean/max phase difference and max magnitude
difference in dB, labelled `HFSS_VALIDATED_OFF_RESONANCE` if within tolerance (1.0° phase, 0.05 dB)
or with the measured residual otherwise. Every comparison is appended to a ledger you can browse
on the same page. Only case G1 has a real HFSS export bundled with this checkout (`hfss_exports/`,
excluded from git — `docs/DEVELOPING.md`); the others need your own export.

## 8. HMSA designer

**Step 1 — import**: upload an HFSS field-calculator `.fld` phase export (magnitude optional).
Reports grid size, detected point ordering, the one documented feed-point NaN, grid regularity,
and — if you give substrate parameters and a frequency — a sanity-check comparison of the
field's own measured surface-wave wavenumber against the analytic bare-slab TM0 value (never a
requirement, just a consistency check).

**Step 2 — synthesize**: give one or more object-beam directions (`theta,phi` in degrees, one per
line), a cell period, and pick a formulation:
- **B (phase-difference)** — needs only phase (works with the example file as-is); picks a cell
  from a library by nearest phase match, same machinery as the Reflectarray page. This is the
  formulation this project's own real data (a phase-only field export) actually uses.
- **A (impedance-modulated, Fong et al.)** — needs a *complex* reference field (upload a matching
  magnitude `.fld` too) plus fabrication-calibration constants (`X_avg`, modulation index) this
  project has not independently derived (labelled `EXPERIMENTAL`) and an HFSS eigenmode dispersion
  table (CSV: `geometry_param, beta_rad_per_m, x_s_ohm`) instead of a phase library.

Both report `CELLS_PER_FRINGE(...)` if the chosen period undersamples the reference field's own
measured fringe spacing (recommended ≥ 6 cells/fringe for a surface-wave hologram) and
`FEED_REGION_MASKED` for the interpolated feed-point NaN. Formulation B's result links to the
Pattern page the same way a Reflectarray design does.

## Known limitations (see also `docs/phases/PHASE13_WEBAPP.md`)

- **No plots.** Every page shows tables, not the server-rendered charts the original brief called
  for. Large results are capped (e.g. HMSA's 22,801-cell grid shows only its first 20 rows).
- **Export does not accept an HMSA design.** HMSA's per-cell result carries one
  `geometry_param` (from the library match), not the separate chosen length/width pair Export's
  job function expects from a Reflectarray design's `CellAssignment`.
- **No run manifest or diagnostic-bundle download**, despite `src/manifest.py` existing and being
  tested — no router currently attaches one to a job or exposes it for download.
- **Job results are held in memory for the process's lifetime** — no eviction; a long session with
  many large sweeps will grow memory usage. Restart the app to clear it.
- **Cancellation is cooperative and only checked between whole solves** in a loop; the Library
  sweep's own progress callback does not check it at all yet.
