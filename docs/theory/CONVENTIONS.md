# Numerical and sign conventions

Single source of truth. If code and this file disagree, that is a bug; fix one
of them and add a test. Nothing below may change silently between modules.

## Units
- Internally **SI only**: metres, hertz, seconds, radians, ohms.
- mm / GHz / degrees appear only at the input/output boundary (`src/units.py`,
  `*_deg` config keys). Plain numbers in `config.yaml` are SI. A unit suffix
  (`"8 mm"`, `"10 GHz"`) converts on load; suffixes are case-sensitive.
- Constants (`src/constants.py`): `C0` exact; `MU0` CODATA 2018
  (1.25663706212e-6 H/m); `EPS0 = 1/(MU0*C0^2)`; `ETA0 = MU0*C0`.

## Time convention
- **`exp(+jωt)`** everywhere. A phasor `F` represents `Re{F·exp(+jωt)}`.
- Consequences (all tested):
  - Lossy dielectric: `ε_c = ε_r (1 − j·tanδ)`, so `Im(ε_c) ≤ 0`, `Im(k) ≤ 0`.
  - A wave travelling toward `+z` is `exp(−jkz)`; toward `−z` is `exp(+jkz)`.
  - Plane wave: `exp(−j k·r)`.
- Tools using `exp(−iωt)` (most physics texts, some solvers) report `Im(ε) > 0`
  and conjugated S-parameters. Conjugate when importing/comparing.

## Coordinates and stack (Phase 1)
- Unit cell in the x–y plane, origin at the cell centre,
  `x ∈ [−Px/2, Px/2]`, `y ∈ [−Py/2, Py/2]`. Rectangular lattice.
- Metal pattern lies in the plane `z = 0` (zero thickness, PEC).
- Stack, top to bottom: **air (z > 0, incident side)** | pattern at z = 0 |
  **substrate (−h < z < 0)** | **air (z < −h)**. No ground plane (default), no superstrate.
  **`backing`** (Phase 12A; `UnitCell.backing`, `LayeredGreenFunction.backing`, config key
  `substrate.backing`): `"air"` (default, every pre-Phase-12A config/call site, unchanged) is the
  stack above; `"pec"` replaces the bottom air region with a **PEC ground plane at z = −h**
  (air | pattern | substrate | ground). A grounded stack has **no bottom port**: `S12`, `S21`,
  `S22` are undefined (not zero — reading them raises), `ReferencePlanes.bottom_air_thickness` is
  `None`, and `reference_planes.bottom_air_thickness` must be **absent** from the config (strict:
  present and grounded is an error). See GREEN_FUNCTION.md Sec. 12A and
  PHASE12A_GROUND_PLANE.md.
- Rotation of a feature is counter-clockwise (from +x toward +y) about its centre.

## Incident wave
- Travels toward `−z`; `k_inc = (kx0, ky0, −kz0)`.
- `θ` = angle of incidence from the surface normal; `φ` = azimuth of the in-plane
  part of `k_inc` from +x toward +y:
  `kx0 = k0 sinθ cosφ`, `ky0 = k0 sinθ sinφ`, `kz0 = k0 cosθ`.
- Polarization = direction of the incident E-field (`x` or `y`), must be ⟂ `k_inc`.
- **TE / TM (Phase 20 R1)**, relative to the plane of incidence (azimuth `φ`): `p = (cosφ, sinφ)`, `q = ẑ × p = (−sinφ, cosφ)`.
  TE: `E ∥ q` (⟂ plane of incidence). TM: tangential `E ∥ p` (E in the plane of incidence; `E = (cosθ cosφ, cosθ sinφ, sinθ)`).
  At `θ = 0, φ = 0`: **TE = y, TM = x**; in general the pair is the x/y pair rotated by `φ`. S-parameters at oblique
  incidence are given in this basis, index 0 = TE, 1 = TM (`ObliqueSParameters`, S_PARAMETERS.md Sec. 11).
  Circular-polarization labels are defined in DATABASE.md Sec. 7.

## Floquet modes and the kz branch
- Quasi-periodicity: `E(x+Px, y) = E(x, y)·exp(−j·kx0·Px)` (same for y).
- Harmonics `exp(−j(kxm·x + kyn·y))`,
  `kxm = kx0 + 2πm/Px`, `kyn = ky0 + 2πn/Py`.
- Spatial transform pair implied by these harmonics:
  `f(x,y) = Σ f̃_mn exp(−j(kxm x + kyn y))`, `f̃_mn = 1/(Px·Py) ∬ f exp(+j(kxm x + kyn y)) dA`.
- `kz = sqrt(k² − kxm² − kyn²)` with **`Im(kz) ≤ 0`** and `Re(kz) ≥ 0`
  (propagating lossless: real positive; evanescent lossless: `−j·|·|`).
  Scattered fields go as `exp(−j·kz·z)` for `z > 0` and `exp(+j·kz·z)` for `z < −h`,
  so both decay when evanescent.

## Green's function (Phase 2; details and derivation in GREEN_FUNCTION.md)
- Electric surface current `J` [A/m] lies in `z = 0` (air/substrate interface).
- **`Ẽ_t = G · J̃`** with `G` in ohms. `G` is the *physical field* response, so
  `G = −Z` (free space: `E = −(η0/2) J`). Sign is fixed here; do not redefine it.
- TE/TM directions: `p = (kx,ky)/kt` (TM), `q = ẑ×p = (−ky,kx)/kt` (TE);
  `V, I` of the transmission lines are defined by `E_t = V ê`, `H_t = I ẑ×ê`;
  the sheet is a shunt source with `I(0⁺) − I(0⁻) = −J_e`.
- `G` contains no `1/(Px·Py)`; that belongs to the Floquet coefficients.
- At `kt = 0`, `Z_TE = Z_TM` and `G = −Z·I₂`.

## MoM (Phase 3; details and derivation in MOM_FORMULATION.md)
- Unknown: PEC surface current on the metal only (slots = metal complement, no magnetic current), as x/y
  **rooftop** coefficients `I` [A/m] (value at the mesh edge). Ordering: all x-directed, then all y-directed.
- **`Z I = V`**, `Z_ij = −(1/A) Σ_mn conj(Φ_i) G_{d_i d_j} Φ_j`, `V_i = conj(Φ_i(k_00)) E_stack,d_i`.
  `Z` is an impedance (Hermitian part PSD), so the sign of `G` is carried into `Z`. Galerkin inner product has
  **no conjugation** (`Z` complex-symmetric at normal incidence).
- Rooftop transform `Φ(k) = ∬ f e^{+jk·r} dA`; seam (periodic) edges are referenced to the copy at `−Px/2`, `−Py/2`.
- Incident field: unit amplitude (1 V/m) with phase reference **at the pattern plane z = 0**; `E_stack` is the tangential
  field of the *unpatterned* stack at z = 0 (Airy factor `1+Γ`).
- Mesh cell is metal iff its centroid is metal. Floquet truncation must satisfy `2·order + 1 ≥ max(cells_x, cells_y)`.
- **Reference planes** (explicit in `config.yaml`, used from Phase 4, not in `Z`/`V`):
  `z_top = +top_air_thickness`, `z_bottom = −(h + bottom_air_thickness)`.
- **Mesh alignment** (Phase 12B; details in `discretization.py`): mesh lines sit at `−P/2 + k·P/N`
  for integer `k = 0..N` (`N` = `cells_x`/`cells_y`, `P` the period on that axis). A feature edge
  that does **not** sit exactly on one is staircased (nonzero `metal_area_error`); one that does
  has zero staircase error there. Two geometries that discretize to the *identical* metal mask on
  a given mesh are mesh-identical: they solve to bit-identical `Z`, `V` and hence `S` at every
  frequency, regardless of how different their continuous geometry is
  (`discretization.metal_mask_key`/`duplicate_geometry_groups`; `sweep.SweepPlan.mesh_duplicate_groups`).
  This is a numerical-discretization fact, not a new physical assumption, and never changes the
  centroid-is-metal classification rule above.

## Waves, ports and S-parameters (Phase 4; details in S_PARAMETERS.md)
- `exp(+jωt)`: `exp(−jkz)` travels +z, `exp(+jkz)` travels −z; delay = negative phase.
- **Ports** are the two reference planes: port 1 at `z_top = +top_air_thickness`, port 2 at
  `z_bottom = −(h + bottom_air_thickness)`, both in air; only the specular (0,0) order at normal incidence is a port mode.
- Wave amplitudes are tangential-E amplitudes referenced to the plane: port 1 incident `a1 e^{+jk0(z−z_top)}` (−z),
  reflected `b1 e^{−jk0(z−z_top)}` (+z); port 2 incident `a2 e^{−jk0(z−z_bot)}` (+z), reflected `b2 e^{+jk0(z−z_bot)}` (−z).
- `b = S a`; `S_ij` is 2×2, **`S_ij[received, incident]`**, index 0 = x, 1 = y, the *same global x̂, ŷ at both ports*.
  Co-pol = diagonal; cross-pol `S21_yx` = y out for x in. Reciprocal ⇒ `S = Sᵀ` (4×4); lossless ⇒ `SᴴS = I`.
- Incident amplitude is *measured* from E and H at the plane (up/down split with η0), not assumed.
- Moving planes: `S11′ = S11 e^{−2jk0 D1}`, `S22′ = S22 e^{+2jk0 D2}`, `S21′ = S12′ = S21 e^{+jk0(D2−D1)}`,
  `D1 = z_top′−z_top`, `D2 = z_bottom′−z_bottom`.
- Fields at any z (`sheet_fields`, `stack_fields`): `E_t` continuous everywhere; `H_t` jumps by `ẑ×ΔH = J` at the sheet
  (`sheet_side` selects the limit at exactly z = 0). Port-2 incident amplitude is 1 at `z = −h`, port-1 at `z = 0`.
- **Grounded stack** (Phase 12A, `backing="pec"`): only port 1 exists. `SParameters.grounded = True`; `s11` is the same 2×2
  block as always, `s12`/`s21`/`s22` are `None` and any attempt to read them (`.component`, `.matrix()`, ...) raises.
  `E_t = H_t = 0` for `z < −h` (nothing exists below the ground); `port="bottom"` raises on `stack_fields`/
  `stack_tangential_field`.

## Circular polarization, phase and databases (Phase 5; details in DATABASE.md)
- **Circular basis** (global x–y vectors, same at both ports and for incident and scattered waves — R/L labels the
  transverse vector, not the propagation direction): `e_R = (e_x + j e_y)/√2`, `e_L = (e_x − j e_y)/√2`;
  `S_circ = Tᴴ S_lin T`, `T = [e_R e_L]` unitary; entries `S<out><in>_<received><incident>` (`S21_LR`).
- **Handedness under `exp(+jωt)`**: `E(t) = Re{e e^{jωt}}`; `e_R` gives `(cos ωt, −sin ωt)/√2`, i.e. rotates **clockwise seen from +z**
  (x toward −y). Wave handedness by the right-hand rule (thumb along propagation): a wave travelling **−z** (port-1 incident,
  transmitted into port 2) with `e_R` is **RHCP**; travelling **+z** (port-1 reflected, port-2 incident) `e_R` is **LHCP**.
  So transmission `RR` preserves handedness while reflection `RR` reverses it (a PEC sheet has `S11_RR = −1`).
  The optics naming is the opposite and is not used.
- Co-polar = same label received/incident (`xx, yy` | `RR, LL`); cross-polar = the rest.
- Phase is stored as `arg(S)` in degrees on **(−180, 180]** (negative real axis = +180); a delay is negative phase. Unwrapping is
  valid only for < 180° true change per sample; anchored to the first sample.
- Frequency grids and parameter ranges include both end points; the span must be an integer number of steps.
- A database record is the linear-basis 4×4 S-matrix at the stored reference planes; files are deterministic (no timestamps).

## Not yet fixed (decide in later phases, then document here)
- Which handedness naming (IEEE right-hand rule as used here, or the optics naming) the holographic synthesis and any
  measurement/simulation reference will use; the basis vectors themselves are fixed above.
- (Oblique-incidence port definitions: fixed in Phase 20 R1, see above and S_PARAMETERS.md Sec. 11.)
