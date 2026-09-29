# User guide — the local web app

`make app` (or `lwmom-app`) starts the backend and opens `http://127.0.0.1:8000` in your browser.
Everything below runs locally; nothing leaves the machine. No screenshots are included here — no
browser was available in either session that wrote this guide to capture any (see
`docs/phases/PHASE13_WEBAPP.md` "Known limitations").

The home page is a 3-step HMSA (holographic metasurface antenna) designer wizard — the app's main
flow. Everything else (single unit-cell solves, a standalone reflectarray designer, full-wave
verification, HFSS comparison) lives one click away under **More tools**.

## The 3-step wizard

### Step 1 — Reference wave

Upload an HFSS field-calculator `.fld` phase export (magnitude optional — drag-and-drop or
"choose a file"), and the frequency it was exported at (the file itself doesn't record this).
"Analyze" shows immediately: a PNG of the reference phase map, and a plain-language summary line
comparing the field's own measured surface-wave wavenumber against the analytic bare-slab TM0
value for the substrate you'll pick in Step 2 (a sanity check, never a requirement — a mismatch is
shown as a warning, not blocked). The one documented feed-point NaN a real export carries is
interpolated and reported (`FEED_REGION_MASKED`), never silently dropped.

### Step 2 — Beams & board

**Beams**: 1–4 object-beam directions (θ, φ in degrees) with a relative power in dB (converted to
a field-amplitude weight internally — power ratio = amplitude²).

**Board**: a substrate preset (Rogers RT/duroid 5880, RO4003C, FR-4, or Custom εr/loss
tangent/thickness) and a cell size. Aperture size is shown read-only, taken from the uploaded
file's own extent — there is no cropping/resizing logic, so this is informational, not editable
(an honest scope reduction from what was asked for). Frequency is pre-filled from Step 1 and
locked; "Use a different frequency" unlocks it, and a mismatch against Step 1's frequency is
flagged.

**Advanced** (collapsed): upload an HFSS eigenmode dispersion table (CSV) to use the
impedance-modulated method (Formulation A, Fong et al.) — the recommended path when you have real
dispersion data. Without one, the wizard falls back to the phase-lookup method (Formulation B, a
reflectarray-style control variable standing in for the surface-wave one a monopole-fed hologram
actually needs — labelled `EXPERIMENTAL`) and says so in one sentence. Also here: mean reactance
X̄, modulation depth, and whether to rescale the reference wave to the board's own guided
wavenumber β̄ (on by default — using the raw bare-slab wave as-is accumulates a phase error across
the aperture, reported if you turn this off).

The unit-cell data synthesis needs is built automatically the first time you run Step 3 for a given
board/frequency (a small square-patch library, cached — later runs with the same board reuse it
instantly). You never visit a separate library page.

### Step 3 — Pattern & export

"Run synthesis" chains: builds/reuses the board's unit-cell library → synthesizes the hologram →
predicts the far-field pattern → exports the array. Shows the layout preview, the radiation
pattern (a 2-D u–v map in dB, one principal cut per requested beam), directivity, side-lobe level,
achieved beam direction, and how many cells couldn't reach their required value (with the
mean/max phase error). **Download package** bundles the STEP file, layout CSV/JSON, HFSS setup
notes, every plot, and a project record into one zip; the individual files are also linked
separately.

Export/STEP requires a real per-cell L×W assignment, which only Formulation B (phase-lookup)
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

Chains off a completed Reflectarray OR HMSA (Formulation B) design's cell assignments to produce a
full-array STEP file, layout CSV/JSON, and an HFSS-ready README (Floquet-port mode count,
PEC-vs-copper, zero-thickness-vs-drawn-thickness lessons from Phase 12). Grounded designs only, and
the array must be square (N×N) — the underlying supercell builder's own limitation, inherited
unmodified. No DXF export (an old docstring claimed otherwise; `ezdxf` was never a dependency and
no code reads or writes one).

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
- **The wizard's auto-built library uses a fixed 6-point square-patch family and a fixed mesh** —
  not alignment-searched, not resonance-searched. For a board/frequency combination whose
  resonance falls between two of those points, phase coverage can be poor (the app's own
  `PHASE_COVERAGE_SHORTFALL` label fires honestly when this happens — it is not hidden, but a
  denser or resonance-aware library would give better designs).
- **No run manifest or diagnostic-bundle download**, despite `src/manifest.py` existing and being
  tested — no router currently attaches one to a job or exposes it for download.
- **Job results are held in memory for the process's lifetime** — no eviction; a long session with
  many large sweeps will grow memory usage. Restart the app to clear it.
- **Cancellation is cooperative and only checked between whole solves** in a loop; the Library
  sweep's own progress callback does not check it at all yet.
- **The export package's `project_record.json` captures results, not full inputs** — enough to
  inspect what a run produced, not to resume editing the original beam/board form.
