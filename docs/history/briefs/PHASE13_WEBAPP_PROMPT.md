# Phase 13 brief — production-grade repo + local web app (reflectarray / metasurface design tool)

You are working in `lightweight_mom`. Read `CLAUDE.md`, `CONVENTIONS.md`, `README.md`, `PHASE12D_LOCAL_PERIODICITY.md`
and `PHASE12D_ADDENDUM.md` first. Work in the order below; each part ends with the full test suite green.

## Decision already made (by the user)
Build the app **now, with the capabilities that exist**. The project's gates (`GEOMETRY_INVERSION_LOCKED`,
`STEP_EXPORT_LOCKED`) stay recorded as scientific findings, but they no longer block tooling: **outputs carry honest
validation labels instead of being withheld.** Every design the app produces is explicitly a lookup-table /
local-periodicity design to be verified in full-wave (the standard reflectarray workflow), and the app says so.

Defaults (change only if something forces it, and say why): runs **locally** (Python backend + browser UI on
localhost, no external services, data never leaves the machine). **Two design modes, sharing the same library,
labels and export machinery:**
- **Holographic metasurface antenna (HMSA) — primary.** The reference wave is **imported from HFSS** as a field
  profile (format below). The user's example (`inputs/reference_fields/efl_phase_EzZ_example.fld`) is — confirmed by
  the user — a **short monopole at the centre of a bare grounded substrate (no patches)**: phase of Ez on a plane at
  z = 1.7 mm (RelativeCS2), NaN at the feed (0,0). Measured: a pure cylindrical wave, 2-D coherent fit
  **β = 224.96 rad/m** (coherence 0.9988 over r = 20–140 mm). For RT/Duroid 5880, 1.575 mm the analytic TM0 surface
  wave has β = 224.96 rad/m at **≈ 10.54 GHz** — so either the export frequency is ≈ 10.54 GHz (then HFSS and the
  analytic TM0 agree), or the substrate/frequency differ; the user supplies the frequency, and the app reports this
  consistency check.
- **Reflectarray** (space feed → pencil beams, S11-phase lookup) — secondary, same tooling.
Design point preset: grounded RT/Duroid 5880, h = 1.575 mm, P = 8 mm, square PEC patch, 10–15 GHz — every
parameter user-editable. Frequency is **not** stored in the field file; the user must enter it.

---

## Part 0 — Carry-overs from Phase 12D (small, do first)
1. **Mesh convergence at resonance** (`PHASE12D_ADDENDUM.md` task 1): MoM phase converges only first-order in cell
   size near resonance (G4: 0.2 mm mesh is ~6° off; extrapolated φ∞ ≈ −174.65° at 12.5 GHz). Add a solver option
   `convergence="richardson"` (aligned meshes, fitted order, φ∞ with error bar) and make it the app default for
   single-cell/library solves; keep the fixed-mesh path for speed with a visible "mesh not converged near resonance"
   label.
2. **HFSS data actually on disk**: `hfss_exports/phase12_grounded_8mm/hfss_results/` holds G1
   (`unit_HFSSDesign1.s36p`), G4 (`g4.s36p` = **zero-thickness sheet** run; the 35 µm run exists only as
   `G4_thick35um_comparison.csv`), G6 (`g6.s36p`, 35 µm patch). Phase 12D reported G4/G6 as pending — ingest them
   with the new Touchstone reader, compare against converged MoM, record the known budget at resonance
   (mesh ≈ 5.9°, thickness ≈ 3.7°, residual ≈ 3–4°), update the HFSS gate honestly.
3. **Task B of 12D measured the design, not the model.** 186° achieved span (why not the 279° Phase 12B measured at
   8 mm? find out — L range, mesh or frequency) makes ~50 % intended-order power inevitable for *any* method. Add the
   missing comparison: lookup-table **predicted** per-order power vs full-supercell MoM, for N = 6/8/12. That, not the
   absolute efficiency, says whether the lookup-table forward model is adequate. Report it; the app uses the result
   to set its accuracy labels.

## Part 1 — Curate the repository to production-push grade
Constraints: **no change to physics or numerics**; existing results reproduce bit-identically; all 2061+ tests pass.
- **Git**: the folder is not a repo. `git init`, a `.gitignore` that excludes `.venv/`, `__pycache__/`, `.pytest_cache/`,
  OS/IDE files, and large generated data (`results/unit_cell_library/_solve/` ~470 MB checkpoint, `*.s36p`/`*.s4p`
  ~24 MB each, other bulky outputs). Decide per results folder what is source-of-truth vs regenerable; keep small
  reference outputs tests need; document how to regenerate or fetch the rest (Git LFS only if clearly justified —
  explain). Logical commits with clear messages. **Do not push or add a remote** — leave that to the user.
- **Packaging**: `pyproject.toml` (name, version `0.13.0`, Python ≥ 3.11 but tested on the 3.14 venv, dependencies
  split into core / cad / app / dev extras), installable with `pip install -e .[dev,app]`. The package is imported as
  `src.*` everywhere; decide whether to rename to a real package name (e.g. `lwmom`) — if you do, do it mechanically
  in one commit with tests green and a compatibility note; if not, justify.
- **Layout**: move phase write-ups and briefs out of the root into `docs/` (`docs/theory/` for GREEN_FUNCTION,
  MOM_FORMULATION, S_PARAMETERS, CONVENTIONS, DATABASE; `docs/history/` for PHASE*_*.md; `docs/briefs/` for the
  *_PROMPT.md / ADDENDUM files). Fix every link. Keep root: README (short: what it is, install, run app, run tests,
  status, links), CLAUDE.md, CHANGELOG.md (reconstructed per phase), CONTRIBUTING.md, pyproject, config presets in
  `configs/`. Remove `src/postprocess.py` if still empty (or document it). No license file — ask the user (ISRO
  internal); note it in the final report.
- **Quality tooling**: ruff (lint + format) with a pragmatic rule set, applied in its own commit; type hints on the
  public API surface used by the app; `pre-commit` config; pytest markers `slow`, `cad`, `checkpoint`, `hfss_data` so
  a **fast suite runs in < 5 min** and the full suite stays available; CI workflow file (GitHub Actions: lint + fast
  suite on push; full suite manual/nightly) — written, not run.
- **Reproducibility**: pinned lock file (e.g. `requirements-lock.txt` or `uv.lock`), `make`/`just` targets (setup,
  test-fast, test-full, lint, app), and a fresh-clone check: clone to a temp dir, install, run the fast suite.

## Part 2 — Logging, diagnostics and error reporting (for future solver work)
- One logging setup (`logging` + optional JSON formatter), configurable level (`LWMOM_LOG_LEVEL`), rotating file logs
  under a user data dir, a **run ID** per job propagated through solver calls.
- **Solver diagnostics captured on every solve**, at DEBUG plus a structured record: mesh (cells, dx, aligned?),
  Floquet order, unknowns, assembly/solve time, peak memory, matrix condition estimate where cheap, convergence
  data (Richardson points, fitted order, error bar), grating-order / surface-wave margins, duplicate-geometry and
  mesh-alignment warnings, energy-conservation and reciprocity residuals.
- **Error taxonomy**: typed exceptions with stable codes (e.g. `CONFIG_INVALID`, `GEOMETRY_OUT_OF_CELL`,
  `MESH_UNALIGNED`, `FLOQUET_RANK_DEFICIENT`, `SOLVE_FAILED`, `NOT_CONVERGED`, `PHASE_UNREACHABLE`, `CAD_UNAVAILABLE`),
  user-facing message + remedy + debug detail. No silent fallbacks.
- **Validation labels** attached to every result object and shown in the UI, from a fixed vocabulary grounded in the
  project's evidence, e.g. `HFSS_VALIDATED_OFF_RESONANCE`, `HFSS_RESIDUAL_AT_RESONANCE(~4°)`,
  `MESH_CONVERGED`/`MESH_NOT_CONVERGED`, `LOCAL_PERIODICITY_APPROXIMATION`, `PHASE_COVERAGE_SHORTFALL(x°)`,
  `STEEP_GRADIENT(>…°/cell)`, `ZERO_THICKNESS_METAL_MODEL`, `EXPERIMENTAL`.
- **Run manifest** saved with every job: inputs, resolved config, code version/git hash, dependency versions,
  timings, labels, warnings, output file hashes. A "download diagnostic bundle" action (manifest + logs + inputs) for
  bug reports / future solver enhancement.

## Part 3 — The web app (local)
Stack: FastAPI (or equivalent) backend wrapping the existing public APIs — **no physics re-implemented in the app
layer**; a job queue with progress, cancel and result caching (keyed by config hash; reuse existing
database/checkpoint machinery) because solves take seconds to minutes; a browser UI (keep the build simple —
justify the choice; Plotly or similar for plots). One command starts it (`lwmom app` / `make app`) and opens
localhost. Works offline.

### Reference-field import (HFSS field-calculator export) — required
Format of the user's files (example: `inputs/reference_fields/efl_phase_EzZ_example.fld`, 22,803 lines):
```
Grid Output Min: [-150mm -150mm 1.7mm] Max: [150mm 150mm 1.7mm] Grid Size: [2mm 2mm 0mm] CS: "RelativeCS2" Unit: "mm"
X, Y, Z, Scalar data "CmplxPhase(ScalarZ(<Ex,Ey,Ez>))"
-1.5000000000000000e+02 -1.5000000000000000e+02 1.7000000000000000e+00  -3.1073796258277993e+00
...
```
- Parse line 1 (min/max/grid size with units, coordinate-system name, unit) and line 2 (the field-calculator
  expression); then whitespace-separated `X Y Z value` rows. Example: 151 × 151 points, 2 mm pitch, ±150 mm,
  z = 1.7 mm in "RelativeCS2", x-major ordering (x outer loop, y inner) — detect ordering, don't assume it.
- Interpret the expression: `CmplxPhase(...)` → phase in **radians** on (−π, π]; `CmplxMag(...)`/`Mag` → magnitude;
  `CmplxReal`/`CmplxImag` → parts; `ScalarX/Y/Z(<Ex,Ey,Ez>)` → component. Support combining a phase file with a
  matching magnitude file (or Re/Im pair) into a complex field, per component. Unknown expressions → clear error
  (`FIELD_EXPRESSION_UNSUPPORTED`), never a guess. Phase-only input is allowed: amplitude then defaults to uniform
  and the result carries a `REFERENCE_AMPLITUDE_ASSUMED_UNIFORM` label.
- NaN/missing values (the feed point) → masked, reported, excluded from synthesis (`FEED_REGION_MASKED`).
- Units → SI; record the CS name and plane height; the app asks the user where that CS origin sits relative to the
  metasurface origin (default: coincident, labelled) and the frequency.
- Validation & display: plot unwrapped phase/magnitude maps; grid-regularity check; estimate the local wavenumber
  β(r) = −dφ/dr (2-D gradient) and n_eff = β/k0, and compare with the analytic TM0 of the bare grounded slab
  (for RT/Duroid 5880, 1.575 mm: β ≈ 213 rad/m at 10 GHz, 269 at 12.5 GHz, 326 at 15 GHz — recompute, don't
  hard-code) as a sanity check, not a requirement.
- Resampling to unit-cell centres: complex (not phase) interpolation; warn when the phase advance per cell exceeds
  a documented limit (aliasing). **Note:** with β ≈ 225 rad/m, an 8 mm cell advances ~103° per cell and the
  interference fringes (period ≈ λ_sw ≈ 28 mm for broadside) get only ~3.5 cells — the app must compute
  cells-per-fringe and warn; surface-wave holograms normally use much smaller periods (≈ λ/6–λ/10). The period
  preset must not be silently assumed adequate for HMSA.
- Tests: parser round-trip on the example file (counts, grid, NaN handling, units, ordering), synthetic fields
  with known β, expression-interpretation table, malformed-file errors.

### HMSA synthesis — what exists and what doesn't (be explicit in the app)
- Hologram = interference of the imported reference wave with the object wave(s) (one or more beam directions,
  polarization). Implement both the classic modulated-surface-impedance form (Fong et al. 2010:
  X(x,y) = X̄[1 + M·Re(ψ_obj ψ_ref*)/|…|]) and the Phase 8B/8E phase-difference form, clearly named.
- **The per-cell control quantity for a surface-wave hologram is the unit cell's surface reactance / surface-wave
  index (dispersion), not normal-incidence S11 phase.** The solver has no dispersion/eigenmode capability yet
  (Phase 8B finding; the ground plane now exists). Provide, in order:
  1. **Import of an HFSS eigenmode dispersion table** (geometry parameter vs β or X_s at the design frequency; CSV)
     as a library type, so HMSA designs can be made now with user-supplied data (label `DISPERSION_FROM_HFSS`).
  2. **A MoM dispersion capability** (surface-wave eigenmodes of the periodic grounded cell: real kx > k0 in the
     existing Floquet framework, locate the modal zero/pole), validated against the analytic TM0 of the bare
     grounded slab and against any imported HFSS table before use (label `DISPERSION_FROM_MOM`, `EXPERIMENTAL`
     until validated). If it can't be finished cleanly in this phase, ship (1) and leave (2) as a documented,
     tested scaffold — don't block the app on it.
- **Reference wave on the patterned surface ≠ bare-slab wave.** The imported field is the bare-slab TM0 wave; once
  patches are present the guided wavenumber becomes β̄ (set by the mean surface impedance X̄ of the hologram). Using
  the bare-slab phase directly accumulates a phase error ≈ (β̄ − β_slab)·r across the aperture. The app must offer
  (and default to) re-scaling the reference phase to β̄ from the dispersion library (ψ_ref ≈ A(r)·e^{−jβ̄r}, amplitude
  from the imported field where available), show the accumulated error if the raw field is used, and label it.
  The imported bare-slab field is also a benchmark for the MoM dispersion capability (TM0 at the given frequency).
- Map X(x,y) → geometry per cell from the dispersion library; flag values outside the library's range
  (`IMPEDANCE_UNREACHABLE`), too-few cells per fringe, and the local-periodicity caveat.
- Pattern prediction for HMSA: aperture field from the modulated surface (documented approximation, e.g. the
  standard leaky-wave/aperture-integration model), labelled as approximate; export for HFSS verification.

Pages / features (MVP, in this order):
1. **Unit cell** — stack (grounded / ungrounded), substrate, period, patch L×W (and existing dipole/slot types),
   frequency sweep, polarization, mesh mode (auto-aligned + Richardson / fast). Outputs: S11 (and S21 when
   ungrounded) mag/phase, co/cross-pol, circular basis, grating/surface-wave warnings, validation labels, STEP export
   of the cell.
2. **Library** — L (and/or W) sweep → phase-vs-size curves per frequency, coverage, largest gap, resonance, loss;
   saved as a `UniversalDatabase` entry; browse existing entries (incl. the Phase 7 benchmark, flagged with its
   duplicate-geometry and ungrounded caveats).
3. **Reflectarray designer** — aperture (size / cell count), feed (position, cos^q pattern, polarization), one or more
   beam directions, design frequency. Computes required phase map (standard reflectarray formula, documented),
   selects cells from a chosen library (nearest phase, magnitude-aware), shows: layout, required vs achieved phase,
   quantization/coverage-shortfall map, per-cell phase-gradient map against the measured local-periodicity limit
   (12D: ≤ 5.8°/cell off-resonance; unresolved near resonance), and labels accordingly.
4. **Pattern prediction** — far field from the aperture (reuse Phase 8E aperture Fourier-optics code): patterns,
   beam pointing, SLL, directivity estimate, spillover/illumination efficiencies; bandwidth sweep over the library's
   frequencies. Labelled `LOCAL_PERIODICITY_APPROXIMATION`.
5. **Verify** — pick a sub-region / phase-gradient cut and run the full-supercell MoM check (12C/12D code): lookup
   prediction vs full solve, per-order power. Shows cost estimate before running; refuses sizes known to be
   infeasible (n ≥ 7 supercells) with the reason.
6. **Export** — full-array STEP (substrate + patches + ground; test on a realistic array size and report file size /
   time; stream or tile if large), layout CSV/JSON, DXF of the patch layer if `ezdxf` is available, and an HFSS-ready
   README (setup lessons from Phase 12: Floquet port ≥ 18 modes at 5 mm air, PEC vs copper, zero-thickness vs 35 µm,
   reference plane). Every export embeds the run manifest and the validation labels.
7. **HFSS compare** — upload Touchstone (.sNp) + choose the matching case → overlay vs MoM, error metrics, feed into
   the HFSS gate records.
8. **HMSA designer** — import reference field (above) → object beam(s) → hologram → impedance map → geometry via the
   dispersion library → layout, fringe/cell warnings, pattern estimate, export. Uses the same labels, manifests and
   exports as the reflectarray designer.

App-level tests: API tests for every endpoint (small, fast configs), job lifecycle (queue, progress, cancel, cache
hit), error-code propagation to the UI, manifest completeness, end-to-end smoke tests (headless browser if
available; otherwise API-level) that (a) design a small reflectarray and export it, and (b) import the example
`.fld` file, build a hologram with a synthetic dispersion table, and export it.

## Deliverables
- Parts 0–3 as above, each in its own commits. `docs/USER_GUIDE.md` (with screenshots if you can capture them),
  `docs/ARCHITECTURE.md` (backend/job/UI/solver layers, data flow, where logs and manifests go),
  `docs/DEVELOPING.md` (setup, tests, markers, adding a solver capability, reading diagnostic bundles).
- `docs/history/PHASE13_WEBAPP.md` with what was built, carry-over results (Part 0 numbers), known limitations and the
  validation-label definitions.
- Final report: test counts and runtimes (fast and full), app start command, fresh-clone result, repo size, anything
  excluded from git and how to regenerate it, open questions for the user (license, remote, package name if changed).

## Rules
- Physics/numerics untouched except the Part 0 convergence option (additive, default only in the app).
- Measured numbers only; no fabricated HFSS data; labels must be backed by recorded evidence.
- Never hide a limitation to make the UI look finished — surface it as a label or warning.
- Don't push, don't add remotes, don't add a license without the user.
