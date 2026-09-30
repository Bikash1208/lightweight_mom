# User guide — the local web app ("HMSA Designer"; `lwmom` is the Python package)

`make app` (or `lwmom-app`) starts the backend and opens `http://127.0.0.1:8000` in your browser.
Everything below runs locally; nothing leaves the machine. A real walk-through's own screenshots
are checked into `docs/screenshots/`; fresh ones after Part 3D's fixes were not captured — no
browser or screenshot tool was available in the session that made those fixes (see
`docs/phases/PHASE13_WEBAPP.md` "Known limitations").

The home page is a 3-step HMSA (holographic metasurface antenna) designer wizard — the app's main
flow. Everything else (single unit-cell solves, a standalone reflectarray designer, full-wave
verification, HFSS comparison) lives one click away under **More tools**.

## The 3-step wizard

### Step 1 — Reference wave

Upload an HFSS field-calculator `.fld` phase export (phase; magnitude optional — drag-and-drop or
"choose a file"), and the frequency it was exported at (the file itself doesn't record this).
"Analyze" shows immediately: a PNG of the reference phase map, and a plain-language summary line
naming the actual substrate you'll pick in Step 2 and comparing the field's own measured
surface-wave wavenumber against that substrate's analytic bare-slab TM0 value (a sanity check,
never a requirement — a mismatch is shown as a warning, not blocked; this line re-checks itself
whenever Step 2's board changes). The one documented feed-point NaN a real export carries is
interpolated and reported (`FEED_REGION_MASKED`), never silently dropped.

### Step 2 — Beams & board

**Beams**: 1–4 object-beam directions (θ, φ in degrees) with a relative power in dB (converted to
a field-amplitude weight internally — power ratio = amplitude²).

**Board**: a substrate preset (Rogers RT/duroid 5880, RO4003C, FR-4, or Custom εr/loss
tangent/thickness) and a cell size (a suggested value appears once Step 1's field is analyzed,
targeting 8 cells per surface-wave fringe at this board's own measured wavenumber). Aperture size
is shown read-only, taken from the uploaded file's own extent — there is no cropping/resizing
logic, so this is informational, not editable (an honest scope reduction from what was asked for).
Frequency is pre-filled from Step 1 and locked; "Use a different frequency" unlocks it, and a
mismatch against Step 1's frequency is flagged. As soon as board, cell size and frequency are all
set, a line shows the phase-lookup method's achievable phase range for that combination (built via
the same auto-library call Step 3 would make, cached) — a narrow range (well under 360°) is common
for small cells and is flagged, not hidden.

**Advanced** (collapsed): an **Accuracy** control (Fast/Accurate — Accurate uses a finer mesh for
the board's own unit-cell data; slower, less staircasing error, not a substitute for HFSS
validation, not Richardson extrapolation). Upload an HFSS eigenmode dispersion table (CSV) to use
the impedance-modulated method (Formulation A, Fong et al.) — the recommended path when you have
real dispersion data; its own fields (mean reactance X̄, modulation depth) only appear once a table
is uploaded. Without one, the wizard falls back to the phase-lookup method (Formulation B, a
reflectarray-style control variable standing in for the surface-wave one a monopole-fed hologram
actually needs — labelled `EXPERIMENTAL`, with a reason shown in Step 3) and says so in one
sentence. Also here: whether to rescale the reference wave to the board's own guided wavenumber β̄
(on by default — using the raw bare-slab wave as-is accumulates a phase error across the aperture,
reported if you turn this off).

The unit-cell data synthesis needs is built automatically the first time you run Step 3 for a given
board/frequency (a small square-patch library, mesh-aligned and cached — later runs with the same
board reuse it instantly). You never visit a separate library page.

### Step 3 — Pattern & export

A one-line summary (beam count, board, cell size, method) appears before you run. "Design antenna"
chains: builds/reuses the board's unit-cell library → synthesizes the hologram → predicts the
far-field pattern → builds the layout export (CSV/JSON/DXF) and HFSS notes. Shows the layout
preview, the radiation pattern (a 2-D u–v map in dB, one principal cut per requested beam),
directivity, side-lobe level, achieved beam direction, and how many cells couldn't reach their
required value (with the mean/max phase error) — if that's more than 20% of the cells, a plain
sentence explains why (the phase-lookup method's limited range at small cell sizes) and what to do
about it (a dispersion table, or a larger cell).

STEP (CAD) is a separate, explicit "Build STEP (CAD)" action — it can take minutes and produce a
file well over 100 MB for a large array, so it no longer runs automatically as part of the chain
above; its own status line shows an estimated progress percentage (from Phase 12D's own measured
runtimes, log-log interpolated — labelled "estimated," never claimed as measured) while it runs.
**Download package** bundles whatever's been built (layout CSV/JSON/DXF, HFSS notes, STEP if you
built it, every plot, and a project record) into one zip; the individual files are also linked
separately.

Export requires a real per-cell L×W assignment, which only Formulation B (phase-lookup)
currently produces — Formulation A's achieved quantity is a surface reactance, not a geometry pair,
so it doesn't chain into Export (Pattern still works for both).

## More tools

Lower-level pages the wizard uses internally, useful directly for a reflectarray (not HMSA) design
or for inspecting a single cell/library/design without going through the full flow. Same
request → poll → result pattern throughout: a status line polls (`queued` → `running` →
`done`/`failed`), then a table of results plus label pills — read these; they're the app's honesty
mechanism, not decoration. A `failed` job's status line shows the error message; the underlying
error also has a stable `code` (e.g. `CONFIG_INVALID`, `SOLVE_FAILED`, `IMPEDANCE_UNREACHABLE`,
`LIBRARY_FREQUENCY_UNAVAILABLE`) if you're scripting against the API directly.

### Unit cell

Solve one rectangular-patch unit cell (grounded, `backing: pec`, one port; or ungrounded, two
ports) at one frequency or a comma-separated list. **Mesh mode**:
- `fast` — solves once at the mesh you specify; always labelled `MESH_NOT_CONVERGED` since it was
  never checked against a finer mesh.
- `richardson` — runs a small family of mesh sizes and Richardson-extrapolates to the zero-cell-
  size limit (Part 0's own addition). The table shows the EXTRAPOLATED S11, not the finest raw
  solve, with an error bar (the size of the extrapolation correction — a fit-covariance number
  from as few as 4 points was found not to be a credible uncertainty); labelled `MESH_CONVERGED`
  or `MESH_NOT_CONVERGED` depending on whether the fit actually converged (a flat or
  non-monotonic sequence is reported honestly, not forced).

"Export STEP" builds the same unit cell's physical geometry (substrate + zero-thickness-modelled
patch, drawn at a configurable thickness for CAD purposes only) and downloads a `.step` file.
Needs the optional `cad` extra (`pip install -e '.[cad]'`) — fails with `CAD_UNAVAILABLE`
otherwise.

### Library

**Browse**: pick a registered library (the read-only Phase 7 benchmark — 324 cells, 4.8 mm,
ungrounded, flagged with its own duplicate-geometry caveat — or one you've built, including any
the wizard auto-built) and a frequency; see its phase-vs-geometry table and coverage (largest
phase gap, span).

**Run a new sweep**: solves every (L, W) combination in the ranges you give, at the frequency grid
you give, for *any* period/substrate/backing — including grounded, which the original Phase 7
library never covered (`docs/DEVELOPING.md` "Why a second library-sweep module"). "Resonant?" is
correct for both grounded and ungrounded cells: a grounded stack's resonance is a phase crossing
through ±180° at nearly constant |S11| (a low-loss reflector barely absorbs anything even at
resonance), not the magnitude peak an ungrounded cell shows — the classifier picks the right test
based on the cell's own `backing`. A small sweep (a handful of L/W values, one frequency) takes
seconds; a fine sweep over many frequencies can take minutes, same cost as running that many
unit-cell solves. The result is registered under the library ID you give it and immediately
available to the Reflectarray page (and, via `ensure_library`, to the wizard).

### Reflectarray designer

Pick a library, a frequency, an array size (cells × cells) and period, a feed position and beam
direction, and a feed-taper exponent `q`. Computes the required phase at each cell (standard
feed-distance-plus-steering formula), picks the nearest-response cell from the library at each
point, and reports achieved-vs-required phase error (mean/RMS/max, and a `PHASE_COVERAGE_
SHORTFALL` label when a cell's error exceeds a documented threshold — not just "some candidate
existed, however poor a match") and the largest per-cell phase gradient against Phase 12D's own
measured 5.8°/cell off-resonance design limit (`STEEP_GRADIENT` if exceeded). The result's job ID
links directly to the Pattern and Export pages with the design pre-selected.

### Pattern prediction

Paste in a design job's ID (from Reflectarray or the wizard's HMSA synthesis) and a frequency;
computes the far field from that design's synthesized aperture (main beam pointing, side-lobe
level, a directivity estimate, a u-v map and per-beam cuts if you supply the beam directions) via
aperture-integration Fourier optics. Always labelled `LOCAL_PERIODICITY_APPROXIMATION` — each
cell's contribution comes from its own isolated-cell library entry, not a full-array solve (see
Verify below for that check).

### Verify

Runs a *real* full-supercell MoM solve of an N-cells-per-2π-gradient-period periodic aperture and
compares its actual per-Floquet-order power (specular/intended order shown alongside every
spurious one) against the idealized (no inter-cell coupling) prediction from the same lookup
table — the one place the app checks its own forward model against ground truth rather than
trusting the library. Costed before running (Phase 12D's own measured runtimes: N=6 ≈7s, N=8
≈18.5s, N=12 ≈124s); sizes beyond N=20 are refused outright as untrusted extrapolation, sizes
between the largest measured point and that limit require an explicit "run anyway" click.

### Export

Chains off a completed Reflectarray OR HMSA (Formulation B) design's cell assignments to produce
layout CSV/JSON/DXF and an HFSS-ready README (Floquet-port mode count, PEC-vs-copper,
zero-thickness-vs-drawn-thickness lessons from Phase 12) immediately, plus an optional full-array
STEP file (a separate, explicit action — can take minutes and produce a large file, with an
estimated-progress status line while it builds). Grounded designs only, and the array must be
square (N×N) — the underlying supercell builder's own limitation, inherited unmodified. The DXF is
a real, hand-written file (no new dependency — `ezdxf` still isn't one) of the patch layer only
(rectangles on one layer), not a full multi-layer PCB/Gerber export.

### HFSS compare

Upload a Touchstone (`.sNp`) export and pick the matching single-cell case (G1–G6); overlays it
against this project's own MoM reference and reports mean/max phase difference and max magnitude
difference in dB, labelled `HFSS_VALIDATED_OFF_RESONANCE` if within tolerance (1.0° phase, 0.05 dB)
or with the measured residual otherwise. Every comparison is appended to a ledger you can browse
on the same page. Only case G1 has a real HFSS export bundled with this checkout (`hfss_exports/`,
excluded from git — `docs/DEVELOPING.md`); the others need your own export.

## Known limitations (see also `docs/phases/PHASE13_WEBAPP.md`)

- **Plots exist for the wizard's own flow only.** The "More tools" pages still show tables, not
  the server-rendered charts the original brief called for; large results are capped (e.g. a
  library sweep's cell table).
- **The wizard's auto-built library uses a 13-point square-patch family, mesh-aligned by
  construction** — denser and free of staircasing error compared to Part 3C's original 6 points,
  but still not resonance-searched. More fundamentally, the phase-lookup method's *achievable
  range* at small cell sizes is a property of the control variable itself (reflection phase, not
  surface impedance) — density and alignment don't fix that; Step 2 now shows the achievable range
  up front so this is visible before you commit to a run, and the `PHASE_COVERAGE_SHORTFALL` label
  (plus a plain-sentence explanation in Step 3 when over 20% of cells are affected) still fires
  honestly when it happens.
- **No run manifest or diagnostic-bundle download**, despite `src/manifest.py` existing and being
  tested — no router currently attaches one to a job or exposes it for download.
- **Job results are held in memory for the process's lifetime** — no eviction; a long session with
  many large sweeps will grow memory usage. Restart the app to clear it.
- **Cancellation is cooperative and only checked between whole solves** in a loop; the Library
  sweep's own progress callback does not check it at all yet.
- **The export package's `project_record.json` captures results, not full inputs** — enough to
  inspect what a run produced, not to resume editing the original beam/board form.
