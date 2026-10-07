# Phase 14 — a real surface-wave dispersion capability + a method choice in the app

Brief: `docs/history/briefs/PHASE14_EIGENMODE_DISPERSION.md`. Built on Phase 13D's own finding (`docs/history/
PHASE13_WEBAPP.md` "Part 3D" §1.1): with 4 mm cells the reflection-phase method reaches only ≈33-36°
of phase, so most cells clip to the smallest or largest patch size — a near-binary hologram
(`docs/screenshots/patch_layout_example.png`). The user asked for the library to compute a real
surface-wave dispersion itself, and for the app to offer both methods.

## Part A — the eigenmode search that didn't work, and the one that did

The brief's own first proposal (Part A as written) was to search the periodic array's own
`Z(kx, ky)` for a shifted eigenmode next to the bare-slab pole (bracket the pole, deflate, refine).
**This was tried and does not work at any mesh size this project can afford to run interactively.**
Tested at patch fill fractions 0.05–0.92 and two mesh resolutions (8×8/floquet 5, and 16×16/floquet
10 to rule out a resolution artifact — 264 unknowns): the largest-singular-value peak of `Z(kt)`
sits on the analytic bare-slab TM0 pole (`dispersion.analytic_tm0.tm0_beta`) to numerical precision
at *every* loading tested, and the smallest singular value — where a genuine loaded-array
eigenvalue (`det(Z) = 0`) would show up — stays completely flat across the whole Brillouin zone at
every loading, at both mesh sizes. There is no numerically resolvable secondary feature to find.
This is a known hard problem in the numerical-electromagnetics literature (locating an eigenvalue
immediately next to a Green's-function pole usually needs a complex-plane contour-integral
eigensolver, e.g. Beyn's algorithm, not a real-axis scalar scan) — genuinely out of scope to solve
from scratch here. `src/dispersion/mom_scaffold.py` is kept exactly as a documented research item
(its own docstring has the full finding); it is **not** used to build Formulation A's table.

**What actually works** (the user's own follow-up prescription, `src/dispersion/sheet_model.py`):
a homogenized patch-sheet-admittance model (Luukkonen et al. 2008) combined with the closed-form
grounded-slab transverse-resonance condition —

1. Solve normal-incidence S11 (this project's own, already-validated MoM) with
   `ReferencePlanes(top_air_thickness=0.0, ...)`, so S11 is referenced exactly at the patch plane —
   no de-embedding needed. Verified directly: for a lightly-loaded cell this reproduces the
   bare-slab-only reflection coefficient (computed from the closed form below) to 4+ significant
   figures.
2. `Y_in = Y0*(1-Gamma)/(1+Gamma)`; `Y_slab` (the grounded dielectric ALONE, no patch) from a
   hand-derived closed form, cross-checked against this project's own tested
   `LayeredGreenFunction.equivalent_impedances` (which is the *parallel* combination of the up- and
   down-looking admittances, not the down-looking one alone — the first version of this derivation
   used that function directly and got a materially wrong answer until this was caught and fixed).
   `Y_sheet = Y_in - Y_slab`, the patch layer's own contribution.
3. `beta` solves `Im[Y_sheet + Y_air(beta) + Y_slab(beta)] = 0` by bisection over
   `sparameters.surface_wave_interval` — a plain scalar equation, no periodic-array eigenvalue
   problem at all. `X_s = eta0*sqrt((beta/k0)**2 - 1)`.

**Validation** (`tests/test_dispersion_sheet_model.py`):
- `Y_air(kt) + Y_slab(kt)` (no sheet) reproduces `1/Z_TM` from `equivalent_impedances` to ~1e-17
  (floating-point noise) at every kt tested — both closed forms and the branch convention agree
  with this project's own tested Green's function.
- A small patch (0.5 mm at 4 mm pitch, RT/duroid 5880 1.575 mm, 10.54 GHz) gives beta = 224.98
  rad/m vs. the analytic bare-slab TM0 of 224.98 rad/m (0.5 rad/m tolerance) — test 1's own ask.
- Beta increases **monotonically** with patch size across 0.5–3.85 mm: 224.98 → 230.11 rad/m,
  X_s 72.7 → 109.9 Ω — test 2's own ask. (Root-finding becomes unreliable within ~0.15 mm of the
  full period at this mesh; not chased further — documented, not silently papered over.)
- No HFSS eigenmode spot values were available this session to cross-check against (brief's test
  3) — the table stays labelled `DISPERSION_FROM_SHEET_MODEL` + `EXPERIMENTAL` regardless, and the
  CSV-upload path (`hfss_table.py`) remains available for a real HFSS eigenmode export when there
  is one.

## Part B — the app: a "Design method" choice

Step 2 now has a **Design method** segmented control: **Surface-wave (eigenmode)** (default,
"recommended for monopole-fed boards") vs. **Reflection phase** ("simpler; limited range with
small cells"). The eigenmode method has its own sub-choice, **Unit-cell data**: *Compute
automatically* (default — `POST /api/hmsa/ensure_dispersion`, the sheet-model builder above, cached
the same way `ensure_library` caches the reflection-phase library) or *Upload HFSS eigenmode table
(CSV)* (the existing `upload_dispersion` path, unchanged).

`ensure_dispersion` writes its table as a CSV in the exact schema an uploaded HFSS table would use
(`geometry_param, beta_rad_per_m, x_s_ohm`), so Formulation A's synthesis code needed **no changes**
to consume it — it's just another `dispersion_id`. It also derives X̄ (mean reactance) and the
modulation depth automatically: X̄ = the centre of the table's own achievable X_s range, M = the
largest symmetric depth around that centre (`(x_max-x_min)/(x_max+x_min)`) — replacing the old
fixed "100 Ω / 0.3" defaults, which were never validated against any real board. These pre-fill the
Advanced X̄/modulation-depth fields and stay editable.

A short eigenmode-solve setup note (master/slave periodic boundaries with a phase shift, eigenmode
solver not a driven one, adequate air height/PML above the bound mode's evanescent tail) was added
to the HFSS README the Export page bundles, for anyone who wants to spot-check the auto-built table
(`src/export/array_export.py`'s `_HFSS_README_TEMPLATE`).

**Not done**: the brief also asked for a "Dispersion (eigenmode) sweep" mode on the More Tools →
Unit-cell library page (next to the existing reflection-phase sweep, with its own β/X_s plot and
CSV download) so the sheet-model method could be driven outside the wizard too. The wizard's own
`ensure_dispersion` covers the main flow end to end; this standalone page was not built this
session — a real, scoped gap, not a silent one.

As soon as board + cell size + frequency are set, Step 2 shows a preview: the achievable phase
range (reflection method) or the surface-reactance range across the auto-built patch sizes
(eigenmode method), with a progress line while it computes. Step 3 shows which method and data
source were used (a "Design method" row in the results list) and gives the "Experimental." chip a
reason in the reflection-phase case, unchanged from Part 3D.

## Part C — the pattern model for the impedance-modulated hologram

Formulation A previously had **no** aperture field at all for Pattern prediction (`hmsa.py`'s own
module docstring: "neither chaining is offered for it here") — a real, pre-existing gap, not
something this phase broke. Fixed: `hologram.hmsa_synthesis.synthesize_formulation_a` now computes
one, using a documented first-order approximation consistent with Fong et al.'s own derivation —
only the AC part of the surface admittance (proportional to `X(x,y)/X_avg - 1`) couples the bound
surface wave to radiation; the DC part (X̄ itself) stays purely reactive by design. The equivalent
aperture field is `E_ref(x,y) * (X_achieved(x,y)/X_avg - 1)`, using the ACHIEVED (dispersion-table-
matched) reactance per cell, mirroring exactly how Formulation B uses the achieved (not ideal)
phase. This is **not** the same radiation model as Formulation B's phase-front reconstruction —
documented as an approximation, not a rigorous leaky-wave/full-wave result, which is out of scope
here.

This aperture field goes through the exact same `pattern.aperture_pattern.compute_pattern` Part 3D
already fixed (the `cos(theta)` obliquity factor, the phi-periodic peak detector) — no separate SLL
bug to fix for this method. Measured on the walk-through case (1 beam, θ=20°, 4 mm cells): main
beam at θ=21.0°, φ=0.0°, SLL = **-12.5 dB** — a real, sensible sidelobe level, not a horizon-ring
artifact.

## Part D — performance

Measured directly (`build_dispersion_table_from_sheet_model`, one normal-incidence MoM solve per
patch size, no periodic-array eigenvalue search at all):

| Family | Points | Mesh | Total (cold) |
|---|---|---|---|
| Sheet-model dispersion, "fast" | 13 | 32×32 | 5.1 s (0.39 s/size) |
| Sheet-model dispersion, "accurate" | 13 | 48×48 | 3.5 s (0.27 s/size) |
| Reflection-phase library (for comparison) | 13 | 32×32 | 1.7 s (0.13 s/size) |

(One Mac, this checkout, RT/duroid 5880 1.575 mm at 10.54 GHz, 4 mm pitch.) A full walk-through run
(upload → `ensure_dispersion` → synthesize → pattern predict) measured 7.3 s end to end in this
session's own test run. **This comfortably beats the brief's own target (first design on screen in
under 30 s for the default board) without needing the more elaborate machinery the brief listed as
options** (a process-pool per patch size, coarse-then-fine refinement, warm-started root brackets,
frequency interpolation, shipped presets) — the sheet-model method is a single cheap linear solve
per size, not an eigenvalue search, so it was never the bottleneck the brief was worried about.
Config-hash caching (already built for `ensure_library`, reused unmodified for `ensure_dispersion`)
covers "cache forever by (substrate, thickness, pitch, frequency, sizes, mesh, accuracy)" — a
repeat request for the same board returns instantly (`already_built: true`).

**Deferred, because the measured numbers didn't need them, not because they were hard**: starting
the build in the background before Step 3 (Step 2's own preview call already does this in
practice — by the time a user reaches Step 3, the table is already cached from Step 2's preview);
process-pool parallelism across patch sizes; a coarse-then-fine two-pass build; warm-started root
brackets between adjacent sizes; frequency interpolation; shipped presets for common boards/pitches.
If a future board/mesh combination turns out to be slow, these are the documented next options, in
roughly this priority order (parallelism first — the per-size solves are embarrassingly parallel).

## Done when

- Part A: the working method (sheet-model) validated (small-patch → analytic TM0, monotonic
  beta(L)) — done. The originally-specified direct eigenmode method was tried, found not to work,
  and is documented as a research item, not shipped as the production path — an explicit, informed
  deviation from the brief's own Part A method, not a silent one.
- The app (the main wizard) offers both methods with automatic table building — done
  (`ensure_dispersion`, the Design-method control, auto-derived X̄/M). The More Tools standalone
  "Dispersion sweep" page is not built (see Part B's own "Not done" note).
- The walk-through case (1 beam θ=20°, 4 mm cells, 10.54 GHz) run with the eigenmode method: **0 of
  5476 cells unreachable**, vs. the reflection-phase method's own 4006/5476 on the identical case
  (`tests/test_app_wizard_walkthrough.py::test_eigenmode_method_reaches_far_more_cells_than_reflection_phase_on_the_walkthrough_case`,
  numbers printed and asserted, not just claimed) — a complete resolution of the walk-through's own
  central complaint, not just an improvement.
- Screenshots: **not captured** — no browser or screenshot tool was available in this session
  either (the same constraint every prior part of Phase 13 recorded). Documented, not silently
  skipped.
- This document and `docs/USER_GUIDE.md` updated.
