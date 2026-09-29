# Periodic layered-medium spectral dyadic Green's function

Implementation: [src/green_function.py](src/green_function.py) (`LayeredGreenFunction`).
Conventions: [CONVENTIONS.md](CONVENTIONS.md) (`exp(+jωt)`, SI, `Im(kz) ≤ 0`).

## 0. Scope, provenance, and honesty about sources

**What it is.** For a unit-cell stack **air (z>0) | substrate (−h<z<0) | air (z<−h)**
and an electric surface current `J` in the plane `z = 0` (the air/substrate
interface, where the PEC pattern lives), the function returns the 2×2 matrix `G`
that maps each Floquet harmonic of `J` to the tangential electric field on that
plane: `Ẽ_t = G · J̃`.

**What it is not.** No MoM matrix, basis functions, S-parameters or geometry sweep (those came in Phases 3–5).
(Phase 4 added tangential `E`/`H` at arbitrary `z`, see Section 12 and S_PARAMETERS.md; the core
of this document, the spectral dyad at `z = 0`, is unchanged.)

**Provenance.** The formulas below are *derived here from Maxwell's equations*
(Sections 2–4), not copied from a text. The method (spectral-domain,
TE/TM-decomposed transmission-line / immittance formulation) is the standard one
in the references of Section 10. **I did not have those texts available while
writing this, so I have not checked the equations against them line by line.**
Instead they are verified numerically against an independent implementation
that uses neither the transmission-line model nor these closed forms
(Section 9): full 3-D plane-wave field matching, agreeing to ~1e-15, plus an
energy-conservation check. Please cross-check against your preferred reference
before relying on it for publication.

## 1. Definitions

All quantities SI, time factor `exp(+jωt)`, fields `∝ exp(−j(kx x + ky y))`.

- `k0 = ω/c0`, `η0 = μ0 c0`.
- Substrate: `ε_c = ε_r(1 − j tanδ)`, `k1 = k0 √ε_c` (principal root: `Im k1 ≤ 0`), `μ = μ0`,
  thickness `h`. Both air regions are vacuum.
- `kt² = kx² + ky²` (real `kx, ky`).
- Normal wavenumbers, both on the documented branch `Im ≤ 0`, `Re ≥ 0`
  (`periodic_bc.normal_wavenumber`):

  ```
  a = kz(air)       = sqrt(k0² − kt²)
  b = kz(substrate) = sqrt(k1² − kt²)
  ```

- In-plane unit vectors for each `(kx, ky)`:

  ```
  p = (kx, ky)/kt          (TM direction: along k_t)
  q = ẑ × p = (−ky, kx)/kt (TE direction)
  ```

- Floquet series for a current periodic up to the Bloch phase of the incident
  wave (CONVENTIONS.md): `J(x,y) = Σ_mn J̃_mn exp(−j(kxm x + kyn y))`,
  `kxm = kx0 + 2πm/Px`, `kyn = ky0 + 2πn/Py`, and likewise
  `E_t(x,y,0) = Σ_mn Ẽ_mn exp(−j(kxm x + kyn y))`.
- **Definition of the Green's dyad:**

  ```
  Ẽ_mn = G(kxm, kyn) · J̃_mn            G in ohms
  ```

  `G` carries **no** `1/(Px Py)`: that factor belongs to the Floquet
  coefficients `J̃_mn = (1/(PxPy)) ∬ J exp(+j k·r) dA` and will appear in the MoM
  assembly (Phase 3).

## 2. TE/TM decomposition and the transmission-line equations

Source-free region, `∇×E = −jωμH`, `∇×H = jωεE`. Write, for one `(kx, ky)`,
`E_t = V(z) ê`, `H_t = I(z) ẑ×ê` (common factor `exp(−j k_t·ρ)` dropped) with

| mode | `ê` | other field components |
|---|---|---|
| TE^z (`E_z = 0`) | `q` | `H_z = kt V/(ωμ)` |
| TM^z (`H_z = 0`) | `p` | `E_z = −kt I/(ωε)` |

Using `(∇×E)_t = ẑ×∂_z E_t + ∇_t E_z × ẑ`, `ẑ×q = −p`, `ẑ×p = q`,
`∇_t → −j kt p`, and eliminating `H_z` (TE) or `E_z` (TM) with the
z-components of Maxwell's equations gives, for **both** modes,

```
dV/dz = −j kz Z_c I ,     dI/dz = −j kz Y_c V ,     Y_c = 1/Z_c ,     kz² = k² − kt²
```

with characteristic impedances

```
Z_TE = ωμ/kz        Z_TM = kz/(ωε)
```

(`ε = ε0 ε_c` in the substrate, `ε0` in air; `μ = μ0` everywhere). Forward
(+z) waves `V ∝ exp(−j kz z)` carry `I = V/Z_c` and positive z-directed Poynting
flux (checked for both modes). These are the standard TL equations under
`exp(+jωt)`.

**Degeneracy.** At `kt = 0`, `kz = k` and `Z_TE = ωμ/k = η = k/(ωε) = Z_TM`, in every
layer. Away from `kt = 0` they differ by `O(kt²)`.

## 3. The current sheet as a shunt source; equivalent impedance

Boundary condition at the sheet (standard, e.g. Balanis, reference 8):
`ẑ × [H(0⁺) − H(0⁻)] = J`, `E_t` continuous. With `J = J_p p + J_q q` and the
field definitions above, `ẑ × ΔH_t = −ΔI ê`, so for each mode

```
I(0⁺) − I(0⁻) = −J_e ,        J_e = J_p (TM) or J_q (TE) ,  V continuous.
```

Looking up from `z = 0` the line is semi-infinite air (impedance `Z_up = Z_c,air`).
Looking down it is the substrate line of length `h` terminated by semi-infinite
air. With `Γ_L = (Z_air − Z_sub)/(Z_air + Z_sub)`, the textbook input impedance
(`exp(+jωt)`, complex `kz` allowed) is

```
Z_down = Z_sub (Z_air + j Z_sub tan(b h)) / (Z_sub + j Z_air tan(b h))
       = Z_sub (1 + Γ_L u)/(1 − Γ_L u) ,       u = exp(−2j b h) .
```

Since `I(0⁺) = V/Z_up` and `I(0⁻) = −V/Z_down`, the jump condition gives

```
V(0) = − Z_eq J_e ,        Z_eq = Z_up Z_down / (Z_up + Z_down) .
```

Sign check: free space, `Z_up = Z_down = η`, gives `E = −(η/2) J`, the textbook
field of a uniform current sheet. So, per mode, `V_TE = −Z_TE^eq J_q`,
`V_TM = −Z_TM^eq J_p`.

## 4. Closed forms actually implemented

The direct `tan(bh)` expression is fine numerically only for small `|Im b h|`
(overflow otherwise) and has removable `0/0` or `1/0` points at `a = 0`
(`kt = k0`) and `b = 0` (`kt = k1`). Substituting the exponential form and
writing `1 − u = 2bσ` with

```
σ = (1 − exp(−2j b h)) / (2 b)        (σ → j h  as b → 0;  computed with expm1)
```

the algebra (first do `Z_eq` in `1 ± u` form, then divide numerator and
denominator by `2b`) gives

```
Z_TE = ω μ0 · [1 + σ (a − b)] / [2a + σ (a − b)²]

Z_TM = a/(ω ε0) · [a − b σ (a − b/ε_c)] / [2a + ε_c σ (a − b/ε_c)²]
```

Equivalent "raw" forms (useful for the limits below), with `A = a`, `B = b/ε_c`:

```
Z_TE = ωμ0  [(a+b) − u(a−b)]        / [(a+b)² − u(a−b)²]
Z_TM = AB/(ωε0) · [(A+B) + u(A−B)]  / [(A+B)² − u(A−B)²]
```

Limits (all tested):

- `h → 0` (`u = 1`): `Z_TE = ωμ0/(2a)`, `Z_TM = a/(2ωε0)` (free space, Z/2).
- `h → ∞`, lossy (`u → 0`): `Z_TE = ωμ0/(a+b)`, `Z_TM = (AB/(A+B))/(ωε0)` (two
  half-spaces in parallel).
- `ε_c = 1`: `b = a`, `B = A`, independent of `h`: free space, for any thickness.
- `kt = 0`: `a = k0`, `b = n k0`, `B = k0/n` and `Z_TE = Z_TM` (both reduce to the
  normal-incidence line `η0 ∥ Z_in`).
- `kt → ∞`: `u → 0`, `Z_TE → jωμ0/(2kt)`, `Z_TM → −j kt / ((ε_c+1) ω ε0)`.
- `kt = k0` (`a = 0`): `Z_TM = 0`, `Z_TE` finite. `Z` is linear in `a`, i.e. has a
  square-root branch point in `kt` (Rayleigh anomaly).
- `b = 0`: `Z_TE = ωμ0 (1 + j a h)/(a (2 + j a h))`, finite.

**Stability.** On the branch `Im b ≤ 0`, `|u| = exp(2 Im(b) h) ≤ 1`, so nothing
overflows, and `σ` is regular, so no cancellation at `b → 0`.

## 5. The dyad

In `(x, y)` components, with `p pᵀ = (1/kt²)[[kx², kx ky],[kx ky, ky²]]`,
`q qᵀ = (1/kt²)[[ky², −kx ky],[−kx ky, kx²]]`:

```
G = −( Z_TM p pᵀ + Z_TE q qᵀ )

G_xx = −(Z_TM kx² + Z_TE ky²)/kt²
G_yy = −(Z_TM ky² + Z_TE kx²)/kt²
G_xy = G_yx = −(Z_TM − Z_TE) kx ky/kt²
```

Physical meaning: `−G` is the impedance dyad; `Re(−J̃ᴴ G J̃)/2` is the power the
current delivers to the fields (Section 9).

### `kt = 0` handling

`p, q` are undefined at `kt = 0`, but `Z_TE = Z_TM ≡ Z`, so
`p pᵀ + q qᵀ = I₂` gives `G = −Z I₂` **independent of any direction choice**.
The kt → 0 limit is direction-independent, with error `O(kt²)`
(`Z_TM − Z_TE = O(kt²)`). The code uses `G = −½(Z_TE + Z_TM) I₂` (exactly
isotropic, `G_xy = 0` exactly) when `|kt| ≤ 10⁻¹² k0`; the neglected difference
there is `~10⁻²⁴` relative. Above the threshold the general formula is used
directly (no cancellation: `Z_TE, Z_TM` are computed separately).
`te_tm_unit_vectors` returns the arbitrary pair `p = x̂, q = ŷ` when asked to
apply a floor.

The Floquet harmonic `(m,n) = (0,0)` at normal incidence is exactly this case.

## 6. Floquet harmonics

`LayeredGreenFunction.floquet_dyads(lattice, kx0, ky0, N)` builds
`kxm = kx0 + 2πm/Px`, `kyn = ky0 + 2πn/Py` for `|m|,|n| ≤ N`
(`RectangularLattice.mode_grid`, `indexing='ij'`, element `[m+N, n+N]`) and
returns `(kx, ky, G)` with `G` of shape `(2N+1, 2N+1, 2, 2)`. Only `kx0, ky0` enter
the incident-wave dependence; the stack itself is independent of it.

## 7. Properties (each is a test)

- **Reciprocity/symmetry:** `G` is complex-symmetric (`G_xy = G_yx`, *not*
  Hermitian when lossy) and even, `G(−k) = G(k)`; together `G(−k)ᵀ = G(k)`
  (Lorentz reciprocity). Mirror `kx → −kx` flips `G_xy` only; axis exchange swaps
  `G_xx ↔ G_yy`; rotation covariance `G(Rk) = R G(k) Rᵀ` (isotropic stack).
- **Passivity:** `Re Z_TE, Re Z_TM ≥ 0` for lossy substrates (`> 0` except
  `Z_TM = 0` at `kt = k0`).
- **Lossless:** `Re Z = 0` for `kt > k0` (reactive), `Re Z > 0` for `kt < k0`
  (radiation). Evanescent TE is inductive (`Im Z_TE > 0`), evanescent TM
  capacitive (`Im Z_TM < 0`) under `exp(+jωt)`.
- **Poles:** for a lossless substrate the denominators vanish at the guided-mode
  (surface-wave) wavenumbers `k0 < kt < k1` of the ungrounded slab. They are
  real-axis poles; `G` is finite except exactly at them and there is **no**
  special-casing. Near-hits by a Floquet harmonic give large values (physical
  resonance); for `tanδ > 0` they move off the real axis.
- **Large `kt`:** `|Z_TE| ~ 1/kt`, `|Z_TM| ~ kt`.

## 8. Implementation notes

`equivalent_impedances` computes `a, b` via `normal_wavenumber` (branch enforced),
`σ` via `-expm1(-2j b h)/(2b)` (with `σ = jh` at `b = 0`), then the two closed
forms of Section 4. `spectral_dyad` broadcasts over any `(kx, ky)` shape and
returns `shape + (2, 2)`. No global state; inputs are validated (`frequency > 0`,
finite; substrate via `Substrate`).

## 9. Verification (all in `tests/`)

| Check | File | Tests |
|---|---|---|
| Dyad vs. **independent** 3-D plane-wave field matching (`H = k×E/ωμ0`, BCs; no TL model), 4 substrates (lossless, lossy, high-ε, very lossy), kt from 0 to 20 k0, both current polarizations, complex `J`, 5 thicknesses; agreement ~1e-15 | `test_green_function_reference.py` | 22 |
| Energy balance `−½Re(J*·E) = P_up + P_down + P_diss` (dissipation by numerical quadrature of `\|E\|²`); outgoing radiation only for `kt<k0`; deliberately wrong `kz` branch is detected | same | (incl. above) |
| Closed forms vs. textbook `tan(bh)` TL formula (lossless, `tanδ` up to 0.3); free-space, `h→0`, `h→∞`, `kt=0`, `kt→∞`, `kt=k0`, `kt=k1` limits | `test_green_function_properties.py` | 83 |
| TE/TM degeneracy at `kt=0`; kt→0 direction-independence and O(kt²) approach; threshold seam | same | |
| Symmetry, evenness, mirror, axis exchange, rotation covariance, Floquet-grid symmetries, indexing | same | |
| `kz` branch (`Im ≤ 0`), `\|u\| ≤ 1`, no overflow up to `kt = 1e8`, reactance signs | same | |
| Floquet-order convergence of a bilinear form with a smooth rooftop×pulse current: partial sums, Cauchy differences shrink with ratio → 4 (tail `~N⁻²`), Richardson extrapolation, lossless `Re` exact after the propagating order, oblique Bloch phase, x/y symmetry | `test_green_function_convergence.py` | 8 |

Mutation checks (done in a scratch copy, not committed): dropping the overall
sign of `G`, swapping TE/TM, dropping `ε_c` in `Z_TM`, flipping the sign of the
exponent in `σ`, disabling the `kt=0` branch, using `ε0` for `μ0`, and flipping
the sign of `G_xy` are each caught by at least one test.

## 10. References (method and conventions; equations here are derived, see Section 0)

1. T. Itoh, "Spectral domain immitance approach for dispersion characteristics of
   generalized printed transmission lines," *IEEE Trans. Microw. Theory Tech.*,
   vol. MTT-28, no. 7, pp. 733–736, 1980. (TE/TM immittance formulation.)
2. K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in
   integral equation formulations," *IEEE Trans. Antennas Propag.*, vol. 45,
   no. 3, pp. 508–519, 1997. (Transmission-line Green's functions for layered media.)
3. R. Mittra, C. H. Chan, and T. Cwik, "Techniques for analyzing frequency
   selective surfaces — a review," *Proc. IEEE*, vol. 76, no. 12,
   pp. 1593–1615, 1988. (Floquet-harmonic spectral-domain MoM for periodic surfaces.)
4. B. A. Munk, *Frequency Selective Surfaces: Theory and Design*, Wiley, 2000.
5. N. Marcuvitz and J. Schwinger, "On the representation of the electric and
   magnetic fields produced by currents and discontinuities in waveguides. I,"
   *J. Appl. Phys.*, vol. 22, no. 6, pp. 806–819, 1951. (TE/TM transmission-line representation.)
6. J. A. Kong, *Electromagnetic Wave Theory*, Wiley, 1986; W. C. Chew, *Waves and
   Fields in Inhomogeneous Media*, IEEE Press, 1995. (Layered-medium TE/TM decomposition.)
7. D. M. Pozar, *Microwave Engineering*, Wiley. (Input impedance of a terminated
   line, `exp(+jωt)`.)
8. C. A. Balanis, *Advanced Engineering Electromagnetics*, Wiley. (Boundary
   condition `n̂×(H₁−H₂) = J_s`; radiation of an infinite current sheet, `E = −ηJ/2`.)

## 11. Hand-off status (updated after Phase 4)

Resolved in Phase 3 (MOM_FORMULATION.md): the MoM equation and its normalization/sign, the Galerkin weighting
and the `1/(PxPy)` factor, the excitation vector, the slot decision (metal-complement currents; no magnetic-source dyad).
Resolved in Phase 4 (S_PARAMETERS.md): `V(z)`, `I(z)` and hence `E_t`, `H_t` at arbitrary `z` for the sheet and for the unpatterned stack.

Still open:
- Acceleration of the algebraically converging Floquet sums (measured `~N⁻²`).
- Guided-mode poles near the real axis for low-loss substrates.

## 12. Fields at arbitrary z (Phase 4)

`LayeredGreenFunction.sheet_fields(kx, ky, J̃, z, sheet_side)` and `.stack_fields(kx, ky, e_inc, z, port)` return the tangential
`(E_t, H_t)` of one Floquet harmonic radiated by the sheet, respectively of the unpatterned stack for a unit incident wave from the
top or bottom. Derivation (substrate solution `V = W[e^{jbz} + ρ e^{−jb(2h+z)}]`, air regions, sheet jump `I(0⁺) − I(0⁻) = −J_e`,
mirror symmetry for the bottom port), stability (all exponents ≤ 0), the `kt = 0` mean-of-modes rule and the `kz = 0` floor are in
S_PARAMETERS.md Sec. 3. At `z = 0` `E_t = G J̃` exactly; verified against the independent field-matching solver to ~1e-14.

## 12A. Optional PEC ground plane (Phase 12A)

`LayeredGreenFunction(substrate, frequency, backing="pec")` (default `backing="air"`, every prior call site unchanged) replaces
the semi-infinite air region below the substrate with a short circuit at `z = −h`: **air (z>0) | pattern at z=0 | substrate
(−h<z<0) | PEC ground (z=−h)**. PHASE12A_GROUND_PLANE.md gives the full numerical verification; this section is the derivation.

**Bottom-face load reflection.** Every field formula in Sections 3-4 and 12 (`_mode_constants`'s `rho`, `Z_down`'s `Γ_L`) is
written in terms of the bottom-face load's reflection coefficient for a wave inside the substrate. For a semi-infinite air load
that is `ρ = (Z_air − Z_sub)/(Z_air + Z_sub)` (Section 3); for a short circuit (`Z_load = 0`) it is exactly `Γ_L = −1`,
independent of mode or frequency. Substituting `Γ_L = −1` into `Z_down = Z_sub(1 + Γ_L u)/(1 − Γ_L u)` gives
`Z_down = Z_sub(1 − u)/(1 + u)`.

**Equivalent impedance seen by the sheet.** `Z_eq = Z_air Z_down/(Z_air + Z_down)` as always (Section 3); substituting the
grounded `Z_down` and simplifying with the same `σ = (1 − u)/(2b) → jh` substitution as Section 4 (so `1 − u = 2bσ`,
`1 + u = 2(1 − bσ)`) gives, after cancelling a common factor of `b`:

```
Z_TE = ωμ0 · σ / (1 + σ(a − b))

Z_TM = A·B·b·σ / (ωε0 · (A − bσ(A − B))),        A = a, B = b/ε_c
```

finite everywhere on the `Im b ≤ 0` branch (including `b = 0`, where `σ = jh`), verified numerically against the direct
`Z_eq = 1/(1/Z_air + 1/(j Z_sub tan(bh)))` textbook form to ~1e-14 across `kt` from 0 through the substrate Rayleigh point, both
polarizations, lossy and lossless (`PHASE12A_GROUND_PLANE.md`). Limits: `h → 0` gives `σ → jh → 0`, so `Z_eq → 0` (the ground
short-circuits the pattern plane directly); `b = 0` gives the finite `Z_TE = ωμ0 σ/(1 + σa)`, `Z_TM = 0` (same as the ungrounded
`kt = k1` limit, since the numerator carries a factor of `B = b/ε_c → 0`); a quarter-wave-shorted-line resonance (`bh = π/2`
lossless, `u = −1`) makes `Z_TE` diverge — a genuine physical resonance (the classic shorted-quarter-wave-line "looks open"
result), **not specially handled**, exactly like the ungrounded slab's own guided-mode poles (Section 7).

**Substrate standing wave and the field below the ground.** Because every substrate/below-ground field formula in Section 12
(`_sheet_mode_transfer`, `_stack_mode_transfer`) is already written purely in terms of the bottom-face reflection coefficient
(called `rho` there), substituting `rho = −1` for `backing="pec"` reproduces the correct grounded fields with **no other code
change**: at `z = −h`, the substrate standing wave's two exponential terms coincide (`e^{jb(-h)} = e^{-jb(2h + (-h))}`), so the
`(1 + rho)` factor that multiplies them is exactly zero there — `E_t(−h) = 0`, the PEC boundary condition, falls out of the
existing algebra rather than being imposed separately. The same `(1 + rho)` factor also multiplies every "z < −h" field term
in both `_sheet_mode_transfer` and `_stack_mode_transfer`, so it is *identically* zero there too: **there is no field below the
ground**, again without a separate "z < −h" code path.

**Top-of-stack reflection Γ (needed only for `stack_fields`/the unpatterned-stack excitation).** The existing formula
`Γ = r(1 − u)/(1 − r²u)` (MOM_FORMULATION.md Sec. 4) assumes the bottom load is air too (`rho = −r`, the *same* magnitude as the
top interface's own reflection `r`, since both interfaces are between the same two media) — it cannot be reused as-is for an
arbitrary bottom load. Re-deriving `Γ = (Z_down − Z_air)/(Z_down + Z_air)` for a general bottom-face load reflection `rho`
(independent of `r`) gives

```
Γ         = (r + rho·u) / (1 + rho·u·r)
1 + Γ     = (1 + r)(1 + rho·u) / (1 + rho·u·r)
```

which reduces exactly to the existing formula when `rho = −r` (ungrounded), and to `Γ = (r − u)/(1 − ru)`,
`1 + Γ = (1+r)(1−u)/(1−ru)` when `rho = −1` (grounded) — verified `|Γ| = 1` for a lossless grounded slab across 10-15 GHz,
`|Γ| < 1` lossy, and `Γ → −1` as `h → 0`, to the solver's own floating-point noise floor.

**Surface-wave pole.** A grounded dielectric slab supports a TM0 surface wave with **no cutoff** (it exists for any `h > 0` at
every frequency, unlike the ungrounded slab's guided modes, Section 7), with transverse wavenumber in `k0 < kt < Re(k_substrate)`.
`sparameters.surface_wave_proximity` checks (as an a priori interval check, in the same spirit as the existing grating-order
check, not a pole root-solve) whether any Floquet harmonic used in the impedance-matrix assembly falls inside or near that
interval; verified for the project's own 4.8 mm and 8 mm periods across 10-15 GHz at normal incidence to have at least one full
`k0` of margin (PHASE12A_GROUND_PLANE.md). `compute_s_parameters` issues a `UserWarning` when it does (like the existing
grating-order warning; `warn_surface_wave=False` suppresses it, same convention as `warn_grating=False`), and `runner.run_sweep`
records the result once per grounded cell, at the sweep's own highest frequency, as `CellEntry.validation["surface_wave"]`
(`None` for an ungrounded cell) -- not folded into the strict, fixed-key-set per-record `Record.diagnostics` schema, so that
schema (and every file already written under it) stays exactly as it was before Phase 12A (DATABASE.md Sec. 6A).

**Bottom port.** There is nothing below the ground, so `port="bottom"` raises `ValueError` on `stack_fields` (and hence
`stack_tangential_field`, `scattering_from_fields`) for `backing="pec"`.
