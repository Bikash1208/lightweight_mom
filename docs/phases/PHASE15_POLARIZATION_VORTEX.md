# Phase 15 — per-beam polarization, pencil vs vortex beams, a vector pattern model, multi-feed

Brief: `docs/briefs/PHASE15_POLARIZATION_VORTEX_MULTIFEED.md`. Methodology source: Ghosh, Ghosh,
Saha, IEEE OJAP 2026 (`docs/references/Ghosh_OJAP2026_multipolarized_HMS.pdf`), Table 1 and eqs.
1-17; Minatti et al., IEEE TAP 2016, for the general polarization-control framing.

Every monopole-fed HMSA design this app produces radiates a **vortex (doughnut) beam** in HFSS
although the app's pre-Phase-15 pattern model showed a pencil beam. Cause: the monopole's TM0
surface wave has a *radial* tangential field, so the aperture current is radial,
`J = rho_hat(x,y) * J_s(x,y)`. Decomposed into circular polarization, a radial current carries an
`e^{+-j*phi'}` term — an orbital-angular-momentum (OAM) vortex at the beam direction — which the
old *scalar* pattern model (a bare 2-D Fourier transform, no polarization notion at all) could not
see. This phase adds (A) a per-beam polarization-compensation term in the hologram so the app can
request a genuine pencil beam of a chosen polarization (or deliberately a vortex of a chosen OAM
charge); (B) a **vector** aperture/far-field model, via this project's own layered Green's
function, so the app's own pattern plot finally shows what HFSS shows; (C) multi-feed ("multipole")
support; (D) the app UI for all of the above.

Every sign/handedness derivation below was **pinned by test against this project's own convention**
(`docs/theory/CONVENTIONS.md`, `src/polarization.py`), not copied from the paper's own convention
(the brief's explicit instruction, since the two are not guaranteed to agree) — two of them were
caught wrong during development by exactly that discipline (see Part B).

## Part A — per-beam polarization compensation (`src/hologram/polarization_compensation.py`)

**Derivation.** The field a small patch of current radiates into a chosen far-field polarization
`p_hat` is proportional to `conj(p_hat) . rho_hat(x,y)` (`rho_hat` = the aperture current's own
local radial direction about the feed). For the whole aperture to add constructively into `p_hat`
at a beam direction, each cell needs an extra hologram phase canceling that projection's own phase
before the ordinary object-wave steering phase is applied:

```
phi_ext,b(x,y) = exp(-j*arg(conj(p_hat_b) . rho_hat(x,y))) * exp(+j*l_b*phi'(x,y))
```

(`l_b` = a beam's own requested OAM charge, 0 for "Pencil"). This ONE formula, evaluated for any
`p_hat_b`, reproduces every special case the brief and paper describe, purely algebraically:
Horizontal (`p_hat = x_hat`) gives `conj(p_hat).rho_hat = cos(phi')`, a real number whose `arg` is 0
or pi by sign — the "0/pi flip across the line through the feed perpendicular to x"; Vertical gives
`sin(phi')`; Slant ±45° gives `cos(phi' -+ 45deg)`; a circular `p_hat` gives a spiral
`e^{+-j*phi'}/sqrt(2)`.

**Handedness.** Which of `polarization.E_R`/`E_L` is "RHCP" depends on propagation direction
(`CONVENTIONS.md`'s `exp(+j*omega*t)`). An object (desired-beam) wave radiates toward `+z`, so
`p_hat_b` for "rhcp"/"lhcp" is picked by calling `polarization.handedness(label, "+z")` and
choosing whichever of `E_R`/`E_L` matches — never hardcoded. Result: "rhcp" -> `E_L`, "lhcp" -> `E_R`
(the project's own convention; the paper's own sign need not, and does not have to, agree).

**Sector multiplexing** (`src/hologram/sector_multiplexing.py`, paper eq. 17): a nearest-beam-
azimuth Voronoi partition of the aperture about the feed — reduces to equal quadrants for 4 beams
90° apart, well-defined for any beam spacing. Wired into both hologram formulations via
`multiplexing="sectors"`.

**Feed position**: taken from the uploaded reference field's own NaN/singularity location by
default (`fields.reference_field_grid.feed_position_from_nan_mask`, the centroid of the masked
cells), or an explicit override (`HmsaSynthesizeRequest.feed_x_mm/feed_y_mm`, API only — no Step 2
UI control, a documented scope reduction).

**Compatibility invariant (tested, `tests/test_hologram_polarization_compensation.py`,
`tests/test_hologram_hmsa_synthesis.py`)**: `desired_polarization="natural"` (the `Beam` dataclass's
own default) makes `phi_ext,b` identically 1 for every beam, so `combine_beams_with_polarization`
reduces EXACTLY to the pre-Phase-15 `combine_beams_scalar` sum — existing designs left at the new
defaults reproduce bit-identically.

## Part B — a vector aperture/far-field model (`src/pattern/vector_aperture_pattern.py`)

**Method.** Vectorize the scalar aperture current Formulation A/B already compute
(`hmsa_synthesis`'s `aperture_field`, already documented as `~ J_s`): `J = rho_hat(x,y) * J_s(x,y)`,
optionally tapered by an exponential leakage decay `e^{-alpha*rho}` (the brief's Part B, off by
default — this project does not itself derive alpha from unit-cell data; a caller-supplied
estimate only). Far field via reciprocity: for each `(theta, phi)`, evaluate this project's own
`green_function.LayeredGreenFunction.spectral_dyad(kx, ky)` (validated elsewhere against an
independent field-matching solver, `tests/test_green_function_fields.py`) at the CONTINUOUS point
`kx, ky = k0*sin(theta)*(cos(phi), sin(phi))` — not a free-space dipole with an image — giving the
tangential spectral response `(Ex, Ey) = G(kx,ky) @ J~(kx,ky)`. Apply the project's existing
`cos(theta)` obliquity factor (`pattern.aperture_pattern`'s own, cited-Balanis, PHASE13D-derived —
the standard stationary-phase Jacobian of a spectral-domain far-field extraction, reused, not
re-derived), then project Cartesian -> spherical (Ludwig-3): `E_theta = Ex*cos(theta)*cos(phi) +
Ey*cos(theta)*sin(phi)`, `E_phi = -Ex*sin(phi) + Ey*cos(phi)`.

**Two sign bugs found and fixed during development** (both by checking a known-polarization beam's
radiated field at exact boresight against `polarization.handedness`, not by inspection):

1. `circular_components` initially hand-rolled `E_R = (E_theta + j*E_phi)/sqrt(2)` by analogy with
   `polarization.E_R`'s own *basis-vector* formula. But the correct DECOMPOSITION coefficient uses
   the CONJUGATE-basis formula `a_R = (V_x - j*V_y)/sqrt(2)` (`polarization.amplitudes_to_circular`'s
   own, already-tested convention) — the opposite sign. Fixed by calling
   `polarization.amplitudes_to_circular` directly instead of re-deriving the formula.
2. `copol_label` initially assumed the far field ends up in the SAME `E_R`/`E_L` slot as the beam's
   own `p_hat_b` (Part A). It does not: compensation makes `conj(p_hat_b) . J` coherent, and
   `conj(p_hat_b) . J` is exactly `J`'s decomposition coefficient ALONG `p_hat_b` itself (not its
   conjugate) — so, after fixing bug 1, the far field's own `p_hat_b`-labeled coefficient IS what's
   dominant, with no inversion. (Bug 1 and a since-corrected draft of bug 2's fix had briefly
   cancelled each other out, which is exactly why this needed a numeric check at development time,
   not just algebra on paper.)

**`measure_vortex_winding`**: samples a beam's own co-pol component on a small ring around its
direction (in `(u,v) = (sin(theta)cos(phi), sin(theta)sin(phi))` space, flat and non-singular even
at boresight), reports the net phase winding in cycles (standard discrete winding-number sum of
wrapped consecutive phase differences around the closed ring) and the null depth relative to the
pattern's own peak — "vortex vs pencil" as a measurement, per the brief.

### Validation (the brief's 5 cases, `tests/test_pattern_vector_aperture_pattern.py`, 25 GHz,
RT/duroid-5880-like substrate, grounded)

1. **Natural, single beam θ=20°, φ=0°** -> both circular components show a deep null (< -30 dB
   measured, ~ -53 dB) with `|winding| = 1` — the vortex HFSS shows, now measured by the model too.
2. **RHCP/LHCP/HP/VP compensation, same beam** -> pencil beam, peak at the requested direction,
   XPD >= 15 dB (measured: RHCP/LHCP 21.4 dB, HP 341 dB, VP 62 dB — HP/VP exceed CP because the
   medium's TE/TM anisotropy at oblique incidence mixes the two CIRCULAR components somewhat, a
   real, documented limitation of a scalar/isotropic unit cell, not a bug; see "Limits" below).
3. **An extra OAM charge** -> deep null with matching `|winding|`, measured at broadside (theta=0,
   where the aperture current is exactly axisymmetric, the unconfounded case) with
   `oam_charge=-2` (Natural): both circular components show `|winding| = 2` with deep nulls. The
   winding's SIGN for an arbitrary (polarization, oam_charge, steering angle) combination is a
   genuinely more complex function of all three than a simple additive rule (verified during
   development: off-axis, some combinations of compensation + oam_charge partially cancel rather
   than add) — this project reports what it MEASURES for a given design (per the brief's own
   "measured, not assumed"), rather than asserting a universal sign formula it cannot independently
   verify against HFSS for this option specifically.
4. **Scalar model unaffected** -> `pattern.aperture_pattern.compute_pattern` and
   `phase8e.beam_reconstruction.far_field_pattern` were not modified at all in Part B (only new
   files/functions were added); confirmed directly.
5. **Design-1 qualitative reproduction** (25 GHz, RT/duroid 5880-like, 88x88 cells, X̄=165.43Ω,
   M=0.45, beams HP(-35°,0°), LHCP(-18°,0°), RHCP(25°,0°), VP(45°,0°), synthetic TM0 reference from
   `dispersion.analytic_tm0.tm0_beta`) -> all four beams measured peaking near their own requested
   direction and polarization, with a pencil-like (not deep-null) response there. No gain-agreement
   claim: this project's own dispersion table (`dispersion.sheet_model`) is a different,
   independently-documented approximation from the paper's own gap-parameter-patch/HFSS-eigenmode
   unit cell.

## Part C — multi-feed ("multipole") support (`src/hologram/multi_feed.py`)

`FeedSpec`: position, complex amplitude, and either an uploaded reference field or `None` for an
analytic TM0 cylindrical wave about that feed (`phase8b.reference_wave.monopole_reference_field`
reused UNMODIFIED — its own `rho = sqrt(x^2+y^2)` is already "distance from wherever (0,0) of the
passed-in grid is," so shifting coordinates is the entire generalization needed).

- **Superposition** (paper eq. 5): sum of every feed's field, each resampled to a common grid first
  (complex interpolation, real/imaginary parts separately — never magnitude/phase, which would
  mishandle the ±180° wrap). Known literature caveat (not independently measured here, no number
  invented): summing independent cylindrical waves typically costs some aperture efficiency vs a
  single feed; Part B's vector model is what actually measures the real difference for a given
  design (`tests/test_hologram_multi_feed.py`'s own superposition-vs-segmentation regression test
  confirms the two strategies genuinely produce different designs).
- **Segmentation**: each feed owns its own nearest-distance aperture region (Voronoi partition);
  the combined field is piecewise, not summed.

Both return a plain `ReferenceFieldGrid`, so `hmsa_synthesis.synthesize_formulation_a/b` need NO
changes to consume either one. `feed_xy_m_for_beams` builds the per-beam feed-position sequence
Part A's `combine_beams_with_polarization` already accepts (shipped generally enough in Part A
specifically so Part C needed no further changes there).

**Export**: `export.array_export.write_hfss_step`/`write_hfss_readme`/`export_array` gain an
optional `feeds_mm` parameter — a small marker disk per feed in a "Feeds" STEP body, plus a
coordinates note in the README (the brief: "a small marker sheet or a note... with coordinates" —
this does both). `None` (default) reproduces every pre-Phase-15 export unaffected.

## Part C (gap-fill) — multi-port sub-apertures (`src/hologram/ports.py`, `multi_port_synthesis.py`)

The brief's own Part C asked for more than the single-shared-hologram multi-feed support above:
Mishra, Barad, Sarkar's own methodology (IEEE TVT 2026, `docs/references/
Mishra_TVT2026_multiport_1bit_OCMA.pdf`) is five coaxial monopoles, each feeding its OWN rectangular
sub-aperture (not several feeds illuminating one shared board) — ports 1-4 are 5λ×5λ single-beam
sub-apertures, port 5 a 10λ×10λ quad-beam one.

**Design: N independent single-aperture syntheses, not new synthesis machinery.** A
`ports.PortRegion` is just an axis-aligned rectangle (global coordinates) plus its own feed
position. `multi_port_synthesis.synthesize_ports` crops a native grid to each region, builds that
port's own reference field (`hologram.multi_feed.superposed_reference_field` on a single-feed list
— a port's own field IS exactly that function's "one feed" case, reused unmodified), and calls the
SAME `hologram.hmsa_synthesis.synthesize_formulation_a/b` every single-aperture design already uses,
with that port's own feed and beam list. Part A's own `feed_xy_m` parameter already makes one
synthesis call radial about an arbitrary feed — a port needs no new hologram formula at all.

**Documented scope reduction**: every port of one multi-port design shares the SAME formulation and
dispersion table/phase library (`MULTI_PORT_SHARED_METHOD`) — the brief's "each with ... its own
method" is not implemented as a genuinely independent per-port formulation choice in this version.
Each port DOES get its own feed, beam list, and polarization compensation (every other part of
Part A/B is fully per-port).

**Layout presets** (`ports.layout_preset_2x2_plus_1_large`, `layout_preset_1xn`, `custom_layout`):
Mishra's own 2×2-plus-one-large-port layout tiles exactly (no gap) when `gap_m=0`; `gap_m` shrinks
every port's rectangle by `gap_m/2` on every edge (a uniform margin, not only at shared seams — the
simpler of two possible rules, documented as such).

**Isolation** (`isolation_estimate`): a first-order, explicitly-NOT-full-wave estimate — the
Euclidean edge-to-edge gap between every pair of ports' rectangles, in units of the board's own
surface-wave wavelength `lambda_sw = 2*pi/k_swg`, flagged below `min_lambda_sw` (default 1.0) with a
guard-row suggestion (how much MORE gap would reach the threshold). `PORT_ISOLATION_NOT_FULL_WAVE`
is attached whenever this runs, and the HFSS README's new "Multi-port isolation" section says
plainly that a real isolation number needs the full S-parameter matrix of one lumped/coax port per
monopole — this project reports no isolation number of its own.

**Aperture efficiency per port** (`pattern.vector_aperture_pattern.aperture_efficiency_estimate`):
`directivity_estimate_dbi / (4*pi*area/lambda^2)`, the standard maximum directivity of a uniformly-
illuminated aperture of that port's own physical area (Balanis) — reported per port, from that
port's own measured `compute_vector_pattern` directivity, never assumed.

**Export**: `export.array_export.export_multi_port_array` (a `write_hfss_step`/`write_layout_files`/
a new rectangle-safe `write_layout_dxf_rect` orchestration, parallel to `export_array` but without
its square-board assumption — Mishra's own 5-port board is 4×2 small-ports wide, not square) builds
ONE combined board from every port's own solved cells, snapped onto one shared corner-indexed
lattice (`multi_port_synthesis.combined_board_lw_grid_m`); since each port's own grid is centred
independently within its own rectangle (this project's existing single-aperture convention), two
neighbouring ports do not necessarily share one exact lattice phase at their shared seam — nearest-
cell snapping absorbs that sub-period (<= one period) mismatch for the combined drawing only, a
documented limitation, not a claim of a seamless physical layout at every boundary. The README's new
"Multi-port isolation" section lists one port per monopole and names the S-matrix to check.
`export.array_export.write_layout_files`'s JSON also had a latent `total_cells = n*n` bug fixed here
(silently wrong for any non-square grid, which no pre-Phase-15 caller could ever trigger since every
one of them was already square — a true no-op fix for every existing export).

**Validation** (`tests/test_phase15_mishra_validation.py`, `@pytest.mark.slow`): Mishra's own 5-port
layout, qualitatively, on a small synthetic board (this project has no measurement of the paper's
own unit cell, so no gain-agreement claim is made) — beams (15°,0°), (30°,90°), (45°,180°),
(60°,270°) from ports 0-3 and all four from port 4, LHCP requested throughout: every port's own
vector-model pattern peaks within a qualitative tolerance of its commanded direction, with LHCP
measured (not assumed) as the dominant circular component, and the quad-beam port shows elevated
power at all four commanded directions.

## Part D — the app

Step 2's beam table gains Polarization and Beam type columns (UI default for a NEW row: Horizontal/
Pencil — a UI convenience default, independent of the backend's own bit-identical-preserving
"natural"/Pencil default, the same split used throughout this phase). Advanced gains Multiplexing,
Feeds (with per-feed `.fld` upload or analytic fallback, and a Superposition/Segmentation choice),
Leakage decay, and Pattern model (Vector default/Scalar legacy). Step 3 gains co-pol/cross-pol
u-v maps, a per-beam results table (direction, polarization/beam-type requested, XPD, null depth,
winding, a plain-language note), and a phase-around-beam-axis plot per beam.

`PatternRequest.pattern_model` defaults to `"scalar"` at the API/schema level (so a request that
omits it, including every pre-Phase-15 caller, is unaffected byte-for-byte) — the app's own JS
explicitly requests `"vector"` by default, matching Part A's `Beam.desired_polarization` split.

**Verified live** (`docs/USER_GUIDE.md`'s Step 2/3, screenshots `docs/screenshots/phase15_*.png`):
a real walk-through in a real browser against the live app (chromium via Playwright, installed and
removed again for this verification only — not a project dependency) — uploaded the example field,
set a 2-beam design (RHCP/Pencil and Horizontal/Vortex charge 1), exercised Multiplexing/Feeds/
Leakage decay/Pattern model in Advanced, ran the design end to end, and confirmed the per-beam
table/plots: beam 1 (RHCP) measured a clean 19.8 dB-XPD pencil beam (~0 dB null, ~0 winding); beam
2 (Horizontal + a DELIBERATELY requested vortex, charge 1) measured a genuine -30.2 dB null with
winding +1.0 — the requested vortex, independently measured, not assumed. A real bug was found and
fixed during this walkthrough: switching the Feeds/multi-feed controls re-rendered the beam table
from stale `state.beams` without first reading back the user's own edits, silently discarding them
(`index.html`'s `add-feed-btn`/`feed-count-control` handlers now call `readBeamsFromTable()` first).

## Limits (documented, not silently glossed over)

- **First-order model throughout**: `J_s ~ E_ref*(X_achieved/X̄ - 1)` (Formulation A's own
  documented leaky-wave-perturbation approximation, PHASE14 Part C) is what Part B vectorizes — not
  a rigorous full-wave radiation calculation.
- **Scalar, isotropic unit cells cap achievable XPD.** Part B's own Case 2 measurements (21.4 dB CP
  vs 341/62 dB linear) show the layered medium's TE/TM anisotropy at oblique incidence mixes the two
  circular components somewhat — a real physical limit of an isotropic-cell surface, not a synthesis
  bug. Tensor/anisotropic unit cells (Fong 2010's own extension; Minatti et al. 2016) would reduce
  this further but are not implemented here (a documented next step, not attempted).
- **Multi-feed pattern PREDICTION** (Part B/D) uses only the first feed's own position for the
  aperture current's radial direction — exact for single-feed designs (still the common case); a
  full per-region ρ̂ model for the far-field prediction of a genuinely segmented multi-feed design
  is not implemented (synthesis itself, Part A/C, already assigns each beam its own correct feed).
- **OAM-charge winding sign** is reported as measured per design, not asserted by a universal
  formula (Part B validation case 3's own note) — a real, not merely conservative, limitation: this
  project has no independent HFSS ground truth for the `oam_charge` option specifically to pin a
  universal sign convention against.
- **No explicit feed-position-override UI** in Step 2 (API-only, `feed_x_mm`/`feed_y_mm`).

## Done when

- [x] Part A: per-beam polarization compensation + sector multiplexing, derived and tested, bit-
      identical regression for existing designs.
- [x] Part B: vector pattern model via this project's own layered Green's function; all 5 brief
      validation cases pass; scalar model unaffected.
- [x] Part C: multi-feed superposition/segmentation; feed markers in Export.
- [x] Part D: app UI for all of the above; verified live in a real browser; this doc + USER_GUIDE
      updated; fresh screenshots.
- [x] Fast suite green after every part (`pytest -m "not slow" -q`).
