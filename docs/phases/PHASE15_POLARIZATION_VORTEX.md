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

## Part E — the Coded (1-bit / N-bit) design method (`src/hologram/coded_synthesis.py`)

Mishra, Barad, Sarkar's own methodology (IEEE TVT 2026) is a THIRD design method next to the
eigenmode/dispersion (Formulation A) and reflection-phase (Formulation B) methods this app already
has: a continuous per-cell phase target, quantized to 2 (1-bit) or `2**N` (N-bit) discrete
reflective states.

**The phase formula needs no new code.** Mishra's eq. 1, `psi = -k0(x*sin(th0)*cos(ph0) +
y*sin(th0)*sin(ph0)) + beta_sw*rho + chi*atan2(y,x)` (`chi = +-1` for LH/RH CP — the same spiral as
Part A; `chi = 0` with a half-plane flip for LP), decomposes exactly into quantities this phase
already computes: the steering term is `arg(single_beam_object_field(...))` (identical formula);
`beta_sw*rho` is `-arg(monopole_reference_field(...))`, so adding it is exactly Formulation B's own
`required_phase = arg(E_obj) - arg(E_ref)`; the CP spiral/LP flip term is exactly
`arg(polarization_compensation_term(..., oam_charge=0))` — Part A's own phi_ext. So calling
`hmsa_synthesis.formulation_b_control_map` on a SINGLE beam already returns `psi_m(x,y)` — the
brief's own "unify the code so both methods use one phi_ext implementation" falls out for free
(`coded_synthesis.per_beam_phase_maps_deg`).

**The 1-bit/N-bit combination rule is one formula, not two.** Mishra eq. 3,
`f = sum_m w_m*cos(psi_m)`, eq. 2, `B = 1 if f>0 else 0`; the brief's N-bit generalization,
quantize `arg(sum_m w_m*e^{j*psi_m})` into `2**N` levels. Since
`f = sum_m w_m*cos(psi_m) = Re(sum_m w_m*e^{j*psi_m})`, and nearest-level quantization of an angle
into `{0, 180}` by sign is exactly `sign(Re(.))`, the 1-bit rule IS the `N=1` case of the general
N-bit quantization — `coded_synthesis.coded_phase_map_deg` implements it once (derived, then
confirmed by test against the brief's own `f>0` formula directly, not merely asserted).

**Unit cells, 1-bit**: a `two_state_cells.TwoStateChoice` pairs two real, measured candidate cells
(never assumed 180° apart) — (a) `best_separated_pair` brute-force searches any pool of candidates
(e.g. `candidates_from_square_library`'s own square-patch-library wrapper) for the pair closest to
180°, reporting the ACHIEVED separation; (b) a genuinely new unit-cell geometry,
`geometry.AnnularSlotPatch` (a rectangular `a`×`b` patch with a concentric annular slot `[r, R]`),
duck-typing `Geometry`'s own interface (`.polarity`, `.feature_area`, `.corners()`, `.in_feature()`)
so `UnitCell`/`RooftopBasis`/the MoM solve need no changes at all — validated with a real solve
(`tests/test_geometry_annular_slot_patch.py`, `@pytest.mark.slow`): a lossless grounded reflector
(`|S11| ~ 1`) whose phase measurably differs from a plain patch of the same outer footprint (e.g.
−49.9° vs −43.2° at one tested size — the slot does something physically real, not asserted).
Wired through the full sweep pipeline (`sweep.py`'s `PARAMETER_NAMES`/`base_parameters`/
`geometry_from_parameters` generalized to dispatch on geometry kind; `config.py`'s `geometry:`
section gains a `patch_annular_slot` type) exactly like the square patch, including a real
`library_sweep.build_annular_slot_library_cells` R-sweep (one real solve per outer radius) so a
user can search for a good pairing; (c) `two_state_cells.two_state_from_values` accepts an
HFSS-measured two-state table directly. `separation_warning` reports (never silently accepts) a
pair far from 180°.

**Geometry mapping**: 1-bit designs map their quantized `{0,1}` level map DIRECTLY onto the chosen
two states (no library search — these ARE the two available states); N-bit (`N>1`) designs instead
match the quantized target phase against a CONTINUOUS phase library via the SAME
`map_control_to_geometry_library`/`nearest_response_match` the continuous Formulation B path
already uses — an N-bit design is reachability-reported exactly like a continuous one.

**Pattern and quantization loss**: the same vector model (Part B) — `aperture_field` uses each
cell's REAL achieved phase (the chosen state's own measured value, not the ideal 0°/180° target),
exactly like Formulation B already does. `coded_synthesis.continuous_reference_aperture_field`
builds the hypothetical zero-quantization-error field from the same design's own
`continuous_phase_deg`, so comparing `compute_vector_pattern` of both measures quantization loss
directly. The quantization loss is defined in the "Part E checks" subsection below; the directivity
difference alone is NOT the loss (see there).

**Export**: `export.array_export` generalizes its cell representation to accept EITHER a plain
`(length_m, width_m)` tuple or a `{"kind": "annular_slot", "a_m", "b_m", "R_m", "r_m"}` dict,
per cell, in `write_hfss_step` (two cadquery faces per annular cell — the outer patch-minus-a-disk
ring, plus the inner disk re-inserted as its own island, exactly matching `AnnularSlotPatch.
in_feature`'s own definition), `write_layout_dxf`/`write_layout_dxf_rect` (native DXF `CIRCLE`
entities for the two slot radii), and `write_layout_files` (JSON always, CSV only gains the new
`kind`/`R_m`/`r_m` columns when the grid actually has an annular-slot cell — a plain-rectangle-only
grid still gets the EXACT original 4-column CSV, byte-for-byte). `app.plots.patch_layout_png` draws
both cell types on one board (the brief: "the patch-layout PNG must show both cell types").

### Part E checks: surface wave, reference plane, quantization definition, beam checks

Added after the first Part E commit, as the reviewed requirements for a 1-bit pair before it is trusted:

- **Reference plane.** The auto-built libraries store S11 referenced 5 mm ABOVE the patch
  (`top_air_thickness=5e-3`). The surface-wave model needs the patch plane, so
  `two_state_cells.patch_plane_gamma` de-embeds with `S11_patch = S11_top*exp(+2j*k0*5mm)` (the rule
  of `sparameters.Sparameters.at_planes`). Every report states both: the patch-plane phase and the
  library-reference phase differ by `2*k0*5mm` (127 deg at 10.54 GHz) — a common offset that does not
  change the pattern but does change any absolute-phase check.
- **Surface wave per state.** `dispersion.sheet_model.find_beta_from_gamma` (the Phase 14 model,
  reused unmodified in its equations, no new MoM solve) gives each state's beta and X_s, or reports
  that it carries no TM0 root. Validated: the bare grounded slab's own reflection returns the
  analytic `tm0_beta` (`tests/test_phase15_part_e_checks.py`). Measured on the walk-through
  library (L = 1…7 mm at 10.54 GHz): every state carries the wave, beta 225.0…236.8 rad/m, X_s
  73…146 ohm.
- **Pair rule.** The chosen pair is the largest separation among pairs where BOTH states carry the
  wave AND `|beta_a - beta_b| / mean(beta) <= 5%` (`BETA_SPREAD_LIMIT_FRACTION`). The 5% limit is a
  project choice, not from the paper; the report shows the measured spread against it. Otherwise
  the job fails with `TWO_STATE_NO_VALID_PAIR`. This gate applies to `beta_mode: "state"` only. The
  paper and model modes (below) use one common beta for the whole surface, so the states' own betas
  need not agree.
- **Phase used in psi.** The 1-bit control map uses `psi = arg(E_obj) + beta_bar*rho + polarization
  term` with `beta_bar` = the mean of the two states' betas (Mishra's `beta_sw`, the surface wave
  the impedance surface carries), not the feed-measured field's phase. The aperture amplitude still
  comes from the real reference field.
- **Patch-plane resonance and magnitude warnings.** The open-circuit patch resonance is at 0 deg at the
  patch plane (Huang & Encinar). A state within 30 deg of 0 deg is refused (`NEAR_PATCH_PLANE_RESONANCE`),
  not flagged. Magnitudes differing by more than 1 dB are flagged.
- **+-5% robustness.** Only the chosen pair is re-solved at 0.95f and 1.05f (the library holds one
  frequency). Walk-through pair (L=5 mm, L=7 mm): separation 147.3 deg at 10.013 GHz, 176.3 deg at
  10.54 GHz, 165.0 deg at 11.067 GHz.
- **Edge phase from the beta mismatch.** `(|beta_a - beta_b|/2) * rho_max`, flagged above 30 deg.
  Walk-through: 60 deg at the aperture corners, so that design is flagged.
- **Cell size and lambda_sw.** Reported as `lambda_sw = 2*pi/beta_bar`, `lambda_sw/cell`, and the
  cell pitch. Walk-through: cell 8.00 mm, lambda_sw 27.11 mm = 3.39 cells, lambda0 28.44 mm.
- **Quantization loss, defined explicitly.** The continuous design `E_ref*exp(j*psi_cont)` is
  compared with the ideal quantized levels `E_ref*exp(j*level)` (0/180 deg for 1 bit). Same aperture
  amplitude, so equal input power. Three numbers:
  1. the peak far-field power ratio at each design's own main beam (the main-beam gain loss) — the
     reported "quantization loss";
  2. the aperture coupling `-10*log10(|<E_q,E_c>|^2/(|E_q|^2*|E_c|^2))`, computed without any
     pattern grid;
  3. the state error (ideal vs achieved levels), also an aperture coupling, so a constant phase
     offset common to both states is not counted.
  The textbook 1-bit value is `-10*log10(sinc^2(pi/2)) = 3.92 dB` (`coded_synthesis.
  textbook_quantization_loss_db`). Walk-through measurements: peak 3.89 dB, coupling 3.96 dB, textbook
  3.92 dB. They agree.
- **Why the directivity difference is not the loss.** An earlier version compared directivities, which
  gave about 0.5 dB and did not match the textbook. The cause: directivity is normalised by the power
  a design RADIATES. A 1-bit aperture with 8 mm (0.28 lambda) cells keeps only about 45% of its
  aperture power in the visible spectrum; the rest is evanescent. Its visible-hemisphere integral
  therefore drops by about 3.5 dB, and that cancels most of the peak drop. On the walk-through design
  the radiated power change is -3.49 dB and the directivity change only 0.40 dB. The directivity
  number is still reported, labelled as reference only.
- **Secondary beam and mirror.** The strongest local maximum outside the main lobe's first-null
  radius (`first_null_half_width_deg`, `asin(1.22*lambda/D)`, widened by `1/cos(theta0)`), and the
  level at the mirror direction `(theta0, phi0+180 deg)` from the nearest grid sample (the phi grid
  is 360/73 spaced, so the sample can be up to 2.5 deg off). Walk-through: strongest secondary at
  -4.9 dB at (12 deg, 0 deg), 8 deg from the main beam at (20 deg, 0 deg), so a real shoulder of the
  1-bit beam, not an artifact of the search; mirror beam at -22.0 dB.

**OAM sign rule (measured, submission).** With a monopole feed, a clean circular vortex needs the
OAM charge sign opposite to the handedness in this code's convention: RHCP needs l<0, LHCP needs
l>0. Measured on the walk-through board at 10.54 GHz: RHCP l=-1 gives winding -1 and an on-axis
null of -31 dB (slow test, real library); RHCP l=-2 gives winding -2 and -32 dB; RHCP l=+1 equals
the natural (uncompensated) beam; RHCP l=+2 puts its power into the LHCP pencil component. The app
flags the wrong sign in the Step 2 beam table (`hologram.polarization_compensation.oam_sign_mismatch`).
Handedness names come from `polarization.handedness` (IEEE, +z): the basis entry `L` is RHCP and
`R` is LHCP in this code's convention.

**Status of the Coded method: Experimental.** The pair screen (`two_state_cells.screen_gamma`) refuses a
state that has no TM0 root, lies within 30 deg of the open-circuit patch-plane resonance (0 deg), or lies
outside the sheet model's range (stop band `beta*P/pi >= 0.9`, or period `>= lambda0/2`). Phase is built
with a path-integrated beta: `phi(x,y) = integral of beta(s) ds` along the radial ray from the feed, which
reduces to `beta*rho` for a uniform beta. The walk-through pair shown above predates the screen and should
not be relied on.

**The inductive-at-patch-plane check was dropped.** For a lossless cell `Z = j*eta0*cot(theta/2)`, so two
states 180 deg apart always have reactances of opposite sign. The check would reject every valid pair, so
it is not part of the screen.

**Limits added for Part E checks:** the 5% beta limit and the 30 deg thresholds are project choices;
the resonance convention depends on the reference plane; the +-5% check re-solves only the chosen
pair; uploaded two-state values are taken to be at the patch plane; 1-bit with multiple feeds is
refused (one feed defines the surface-wave reference); the Coded ports path takes one shared two-state
source for every port.

### Paper reproduction at 26 GHz (Mishra, Barad, Sarkar, IEEE TVT 2026)

Script: `examples/make_paper_board_26ghz.py`. Outputs: `docs/submission/paper26/` (`comparison.md`,
`metrics.json`, figures, `step/`). Run time about 25 s after the cell solves.

**Cell 1 is read from the paper's Figure 2, not from its text.** The text calls cell 1 "full metal", but
Fig. 2 draws it in the substrate colour, with copper only in cell 2, and the full-antenna drawing shows
the same bare regions. We therefore model cell 1 as an EMPTY cell: bare RT/Duroid 5880 (1.575 mm) with no
metal. Cell 2 is a 2.17 mm square patch with an annular slot (r = 0.46 mm, R = 0.75 mm). The period is
2.57 mm (lambda0/4.5 at 26 GHz); the paper's table gives it as `a = 2.57 mm`. The cell solve uses the app's
fast mesh (32 x 32). The empty cell's S11 is the bare-stack reflection (zero currents), which the MoM solve
reduces to with no metal.

Cell pair, patch-plane reference (`LIBRARY_REFERENCE_OFFSET_M` = 5 mm de-embedded):

| f (GHz) | empty phase | slotted phase | separation | `|S11|` difference | beta empty (rad/m) | beta slotted (rad/m) |
|---|---|---|---|---|---|---|
| 25 | +56.4 deg | -130.2 deg | 173.5 deg | 0.000 dB | 580.2 | 666.5 |
| 26 | +49.0 deg | -135.2 deg | 175.9 deg | 0.002 dB | 607.9 | 702.3 |
| 27 | +41.3 deg | -139.5 deg | 179.2 deg | 0.004 dB | 636.1 | 738.1 |

The empty cell's beta equals `tm0_beta` at 26 GHz (607.92 rad/m), the bare-slab check. The separation is
within 180 +- 20 deg at all three frequencies, so the paper's own pair is used. Both states have
`|S11|` about 0 dB: a lossless reflectarray cell reflects almost all power by construction.

**Beta modes** (both reported; the coded phase uses one of them per run):

- Paper mode, Mishra eq. (1): `beta_sw = k0*sqrt((eps_r+1)/2) = 689.27 rad/m` at 26 GHz.
- Model mode: beta of the mixed 1-bit surface, `find_beta_from_sheet_admittance` applied to the mean
  of the two states' sheet admittances (equal-area parallel loading). The mean is an assumption of this
  model, not a paper result. It gives 662.07 rad/m at 26 GHz.
- Aperture-edge phase difference between the two modes, `(beta_paper - beta_model) * rho` at the 28.8 mm
  edge of a 5 lambda0 sub-aperture: 44.9 deg (89.9 deg at the 57.7 mm edge of the large port).

**Board.** 5 lambda0 sub-apertures (57.65 mm) with no guard gap, as the paper's 20 lambda0 x 10 lambda0
board. At 2.57 mm pitch the board is 230.6 x 115.3 mm, 90 x 45 cells. The paper's 92 x 46 cells would
need a pitch of about 2.51 mm, so the cell count differs from the paper's by 2 x 1; the pitch is the
paper's stated period and is kept.

**Beams, LHCP.** Each port's beams are LHCP (the paper's polarisation), with the port's own feed, the
analytic TM0 monopole, and two equalisation iterations for port 5. The table of pointing errors,
LHCP/RHCP levels and directivity is in `docs/submission/paper26/comparison.md`. Summary: single-beam
ports point within 0 to 5 deg of the request; port 5 (paper mode) is 6 to 11 deg off, and in model mode
within 3.2 deg. The paper gives target directions, not measured beam angles, in its text (its measured
patterns are in Fig. 5 only), so no measured angle is compared.

**What is not compared, and why.** The paper's realised LHCP gains (14.8 to 16.8 dBic) include losses
the model does not have, so they are set against our lossless directivity (18 to 21 dBi), not against a
realised gain. The paper's isolation (below -30 dB) is a full-wave measurement; this model gives only the
first-order edge-gap flags, so isolation is checked in HFSS with the exported STEP (`docs/submission/
paper26/step/`). The paper's `|S11|` below -10 dB over 25 to 27 GHz is a feed-match quantity; a
reflectarray cell's own `|S11|` is about 0 dB, so the two are not comparable.

**Limits.** The cell 1 reading (empty, not solid copper) is from the figure; if the figure is misread the
pair and every number here change. The model-mode beta depends on the equal-area averaging assumption.
Results are first-order vector-model results (no mutual coupling between sub-apertures, no conductor or
mismatch losses) and need HFSS confirmation.

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
- **Multi-port designs share one formulation/library across every port** (`MULTI_PORT_SHARED_METHOD`)
  — not a genuinely independent per-port method choice, though each port DOES get its own feed, beam
  list, and polarization compensation.
- **The combined multi-port board drawing snaps to the nearest solved cell** at each shared seam
  (`multi_port_synthesis.combined_board_lw_grid_m`) since each port's own grid is centred
  independently within its own rectangle — a bounded (<= one period) visualization/export
  imprecision at port boundaries, not a claim of a seamless physical layout there.
- **Port isolation and inter-port S-parameters are not computed** — only a first-order edge-distance
  heuristic (`PORT_ISOLATION_NOT_FULL_WAVE`); the HFSS README says plainly what full-wave check is
  needed instead.
- **N-bit (N>1) coded designs need an existing continuous phase library** (e.g. the square-patch
  one) to match their quantized target against — only 1-bit designs use the dedicated two-state
  cell search (`two_state_cells`); no N-bit-specific multi-level library-building tool is provided.
- **The annular-slot geometry has no eigenmode/dispersion-table counterpart** (`sheet_model.py`/
  `mom_scaffold.py`'s `square_patch_unit_cell` were not extended) — Mishra's own coded method uses
  reflection phase (Formulation B's own library family), which is the only path implemented for it.

## Done when

- [x] Part A: per-beam polarization compensation + sector multiplexing, derived and tested, bit-
      identical regression for existing designs.
- [x] Part B: vector pattern model via this project's own layered Green's function; all 5 brief
      validation cases pass; scalar model unaffected.
- [x] Part C: multi-feed superposition/segmentation; feed markers in Export.
- [x] Part D: app UI for all of the above; verified live in a real browser; this doc + USER_GUIDE
      updated; fresh screenshots.
- [x] Fast suite green after every part (`pytest -m "not slow" -q`).
