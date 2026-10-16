# S-parameter extraction (Phase 4)

Implementation: [src/sparameters.py](src/sparameters.py) (extraction, data model),
[src/green_function.py](src/green_function.py) (`sheet_fields`, `stack_fields`),
[src/excitation.py](src/excitation.py) (`ReferencePlanes`, `stack_tangential_field`, now with `port`),
[src/solver.py](src/solver.py) (`solve` with several right-hand sides).
Conventions: [CONVENTIONS.md](CONVENTIONS.md). Upstream: [GREEN_FUNCTION.md](GREEN_FUNCTION.md), [MOM_FORMULATION.md](MOM_FORMULATION.md).

**Pipeline:** solved current `I` → mean current `J̃₀₀` → tangential **E and H at the reference planes** (unpatterned-stack field
+ sheet field) → up/down wave decomposition → `S11, S21` (and `S12, S22`), 2×2 in polarization.
Not included: oblique incidence, TE/TM ports, multi-order (grating) ports. (Sweeps, the unit-cell database and the circular basis derived
from these linear-basis S-parameters are Phase 5: [DATABASE.md](DATABASE.md).)

## 0. Provenance and honesty about sources

As in the earlier phases, the equations are **derived here** (Sections 3–5) from the transmission-line (TL) model of
GREEN_FUNCTION.md and are standard in form (multilayer/Airy transmission-line theory, spectral-domain MoM; references in
Section 10). I did not have those texts available to check line by line. They are verified instead by independent means (Section 8):
a 3-D plane-wave field-matching solver that uses no TL model, closed-form Fresnel/Airy slab results, a full PEC sheet, energy
conservation, reciprocity, and mutation checks. Please cross-check against your preferred reference before publishing.

## 1. Conventions (all with time factor `exp(+jωt)`)

Coordinates as in CONVENTIONS.md: pattern plane `z = 0`, substrate `−h < z < 0`, air above and below.

**Reference planes** (from `config.yaml`, `ReferencePlanes`; explicit and never inferred):

```
port 1 = plane z = z_top    = +top_air_thickness                (in the upper air region)
port 2 = plane z = z_bottom = −(h + bottom_air_thickness)       (in the lower air region)
```

**Port modes.** Only the specular `(0,0)` plane-wave order at **normal incidence** is a port mode. Both ports are in air, so the
tangential-E amplitude ratio equals the power-wave amplitude ratio (`|b|²/|a|²` = power ratio, no impedance factors).
Field components are taken along the **same global axes `x̂, ŷ` at both ports** (no axis flip on the far side).

**Wave amplitudes** (2-vectors `(x, y)` of tangential E, in V/m, *referenced to the plane*):

```
port 1:  incident   E = a1 · exp(+j k0 (z − z_top))      travels −z
         reflected  E = b1 · exp(−j k0 (z − z_top))      travels +z
port 2:  incident   E = a2 · exp(−j k0 (z − z_bot))      travels +z
         reflected  E = b2 · exp(+j k0 (z − z_bot))      travels −z
```

Under `exp(+jωt)` a wave `exp(−jkz)` travels toward +z, `exp(+jkz)` toward −z, and propagation delay is a *negative* phase.

**Definition:**

```
[b1]   [S11 S12] [a1]
[b2] = [S21 S22] [a2]            each S_ij is 2×2, indexed  S_ij[received, incident],  index 0 = x, 1 = y
```

- `S11`: reflection at the top plane; `S21`: transmission, incident on the top plane, received at the bottom plane;
  `S12`, `S22` the same for incidence from below.
- **Co-polarized**: `[x,x]` and `[y,y]`. **Cross-polarized**: `S_ij[y,x]` = y-polarized output for x-polarized input (named
  `"S21_yx"` in `SParameters.component`), and `[x,y]`.
- The complete 4×4 matrix over (port1 x, port1 y, port2 x, port2 y) is `SParameters.matrix()`.

Consequences used as tests: for a reciprocal cell `S = Sᵀ` (4×4), i.e. `S12 = S21ᵀ`, `S11 = S11ᵀ` (so `S11[x,y] = S11[y,x]`),
`S22 = S22ᵀ`; for a lossless cell whose only propagating order is `(0,0)`, `SᴴS = I`; with loss `I − SᴴS ⪰ 0`.

**Incident amplitude is measured, not assumed:** `a` in `S = b/a` is obtained from the fields at the plane (Section 5), so a wrong
phase reference or a leaking radiation condition shows up as a wrong `a` or a non-zero "incoming" residual
(`SParameters.incident_amplitudes`, `.incoming_residual`).

## 2. Internal normalization of the two excitations

The MoM excitation is a unit wave; where "unit" is measured is a convention that cancels in `S = b/a` because `a` is re-measured at the plane:

- port 1: incident tangential E has amplitude 1 at `z = 0` (travelling −z);
- port 2: incident tangential E has amplitude 1 at `z = −h` (travelling +z).

So `a1 = e^{+jk0 z_top}` and `a2 = e^{+jk0 d_bottom}` (with `d_bottom = bottom_air_thickness`); tests assert this.

## 3. Tangential E and H at arbitrary z (`LayeredGreenFunction`)

Notation of GREEN_FUNCTION.md: per `(kx, ky)` and TE/TM mode, `E_t = V(z) ê`, `H_t = I(z) ẑ×ê`; `a, b` = kz in air / substrate;
`u = e^{−2jbh}`; `Y = 1/Z_c` the characteristic admittances (TE: `a/(ωμ0)`, `b/(ωμ0)`; TM: `ωε0/a`, `ωε0 ε_c/b`); and

```
ρ = (Z_air − Z_sub)/(Z_air + Z_sub) = { TE: (b − a)/(b + a) ;  TM: (a − B)/(a + B), B = b/ε_c } = −r
```

`ρ` is the reflection coefficient, for a wave *inside the substrate travelling toward the bottom face*, at that face (load = air).

### 3.1 Substrate solution (both sources)

A homogeneous line with an air load at `z = −h` has `V = W [ e^{+jbz} + ρ e^{−jb(2h+z)} ]`, `I = Y_sub W [ ρ e^{−jb(2h+z)} − e^{+jbz} ]`
(`W` a constant). Both exponentials have modulus ≤ 1 on the branch `Im b ≤ 0` for `−h ≤ z ≤ 0`, so nothing overflows.

### 3.2 Field radiated by the current sheet (`sheet_fields`), per unit `J_e`

```
V0 = V(0) = −Z_eq          (Z_eq: Phase 2 closed forms, parallel of the up- and down-looking impedances)

z > 0:      V = V0 e^{−jaz},                        I = Y_air V
−h ≤ z < 0: V = V0 [e^{jbz} + ρ e^{−jb(2h+z)}]/(1+ρu),   I = Y_sub V0 [ρ e^{−jb(2h+z)} − e^{jbz}]/(1+ρu)
z < −h:     V = V0 (1+ρ) e^{−jbh}/(1+ρu) · e^{+ja(z+h)},  I = −Y_air V
```

Checks (all tested): `V(0) = G J` (the Phase 2 dyad); `E_t` continuous at `z = 0` and `z = −h`, `H_t` continuous at `z = −h`;
the sheet condition `ẑ×[H(0⁺) − H(0⁻)] = J`, i.e. `I(0⁺) − I(0⁻) = −J_e`, with `I(0⁻) = −V0/Z_down`,
`Z_down = Z_sub (1+ρu)/(1−ρu)`. At `z = 0` exactly `H` is discontinuous; `sheet_side='above'|'below'` selects the limit.

### 3.3 Field of the unpatterned stack (`stack_fields`), per unit incident mode voltage

Top port (incident at `z = 0`, travelling −z), with `Γ = r(1−u)/(1−r²u)` and `1+Γ = (1+r)(1−ru)/(1−r²u)`:

```
z ≥ 0:      V = e^{+jaz} + Γ e^{−jaz},           I = Y_air (Γ e^{−jaz} − e^{+jaz})
−h ≤ z < 0: V = W_s [e^{jbz} + ρ e^{−jb(2h+z)}],   I = Y_sub W_s [ρ e^{−jb(2h+z)} − e^{jbz}],   W_s = (1+Γ)/(1+ρu)
z < −h:     V = W_s (1+ρ) e^{−jbh} e^{+ja(z+h)},   I = −Y_air V
```

The transmitted amplitude at the bottom face is `W_s(1+ρ)e^{−jbh} = (1 − r²) e^{−jbh}/(1 − r²u)` (the Airy slab transmission).
**Bottom port** (incident at `z = −h`, travelling +z): the unpatterned stack is symmetric under `z → −h − z`, so
`V_bottom(z) = V_top(−h − z)` and `I_bottom(z) = −I_top(−h − z)` (`H` flips sign, `E` does not).
The incident wave must propagate (`|k_t| < k0`; checked).

### 3.4 Assembly into tangential vectors

With `p = k_t/|k_t|`, `q = ẑ×p` and an in-plane input `w` (the current `J̃` or the incident tangential E):

```
E_t = V_TM (w·p) p + V_TE (w·q) q,          H_t = −I_TE (w·q) p + I_TM (w·p) q
```

At `kt = 0` TE = TM and `p, q` are arbitrary; the code replaces the two modes by their mean so `E_t = V w`, `H_t = I ẑ×w` exactly.
An exactly zero `kz` (Rayleigh point, or lossless `kt = Re k_sub`) is replaced by `10⁻¹² k0` (`I = V/Z_c` would be 0/0); the fields are
continuous through it, and near `kz_air = 0` they vary linearly in `kz_air` (square-root branch point in `kt`), which the tests check.

## 4. Excitation for both ports

`stack_tangential_field(green, excitation, port)` = `stack_fields(...)` at `z = 0` (E is continuous there). Top port:
`(1+Γ) V⁺`, as in MOM_FORMULATION.md Sec. 4. Bottom port: the slab transmission `(1−r²)e^{−jbh}/(1−r²u)` (normal incidence,
`= T` of the bare slab, referenced to the two faces). `Z` does not depend on the excitation, so `compute_s_parameters` assembles `Z` once and
solves four right-hand sides (2 ports × 2 polarizations) with one factorization.

## 5. Extraction

For each (port, polarization), with `J̃₀₀ = (1/A) Σ I_b Φ_b(0,0)` (`RooftopBasis.current_coefficients`):

```
E, H = stack_fields(0, 0, e_pol, z_planes, port) + sheet_fields(0, 0, J̃₀₀, z_planes)      at z_planes = (z_top, z_bottom)
E_up   = (E − η0 ẑ×H)/2,      E_down = (E + η0 ẑ×H)/2,      ẑ×H = (−H_y, H_x)
```

(Up-going waves have `H = ẑ×E/η0`, down-going `H = −ẑ×E/η0`.) Then

```
port-1 excitation:  a = E_down(z_top)[pol]     b1 = E_up(z_top),    b2 = E_down(z_bottom)     →  S11[:,pol] = b1/a,  S21[:,pol] = b2/a
port-2 excitation:  a = E_up(z_bottom)[pol]    b1 = E_up(z_top),    b2 = E_down(z_bottom)     →  S12[:,pol] = b1/a,  S22[:,pol] = b2/a
```

Only the sheet-radiated field and the stack field are added, so nothing incoming other than the specified wave may remain:
`E_down(z_top)[other pol]`, `E_up(z_bottom)` (port 1) and the mirror set (port 2) are returned as `incoming_residual`
(`~1e-16` in all tests; it would flag a wrong radiation condition or reference plane).

Only the `(0,0)` harmonic enters: evanescent harmonics do not carry power to a plane, and `compute_s_parameters` warns when other orders
propagate in air (`period ≥ λ0`), since then `S` is specular-only and energy is not conserved in the specular channel.

## 5A. Grounded stack: one port only (Phase 12A)

`backing="pec"` (GREEN_FUNCTION.md Sec. 12A) removes port 2 entirely: nothing exists below the ground, so there is no bottom
reference plane and no wave to measure there. `compute_s_parameters`/`scattering_from_fields` then solve only the 2 (not 4)
right-hand sides (port-1 x and y polarization) and return an `SParameters` with `grounded=True`: `s11` is the same 2×2 block
as always (Section 1's definition, indexed the same way), but `s12`, `s21`, `s22` are `None` — **undefined, not zero** — and
`.component()`, `.matrix()`, `.to_record()` for anything other than S11 raise `ValueError` rather than silently returning 0
(this matters physically: a lossless grounded stack genuinely has `S21 ≡ 0`, but representing "no port 2 exists" as a
measured zero would make that a silent assumption baked into the type rather than a fact the solver actually established).
`incident_amplitudes` and `mean_currents` drop their port axis (shape `(1, 2)` / `(1, 2, 2)`). `power_out_per_incident()`
still works, generalized to sum over the 2 available rows instead of 4 (for a lossless grounded stack this is exactly 1:
everything that leaves the pattern plane returns through port 1, since there is nowhere else for it to go).
`ReferencePlanes.bottom_air_thickness = None` marks a plane set as grounded (`.grounded` property); `_check_planes` and
`at_planes` require this to match the stack's own `backing` (mixing a grounded stack with two-port planes, or vice versa,
raises rather than silently ignoring the mismatch).

## 6. Reference-plane phase convention (explicit)

For planes `(z_top, z_bottom)` and `(z_top′, z_bottom′)`, with `D1 = z_top′ − z_top` and `D2 = z_bottom′ − z_bottom`
(outward motion is `D1 > 0`, `D2 < 0`), the (0,0) waves are propagated through air only:

```
S11′ = S11 exp(−2j k0 D1),      S22′ = S22 exp(+2j k0 D2),      S21′ = S12′ = S21 exp(+j k0 (D2 − D1))
a1′ = a1 e^{+jk0 D1},  a2′ = a2 e^{−jk0 D2}                                 (implemented in `SParameters.at_planes`)
```

Special cases: planes on the substrate faces (`top = bottom = 0` thicknesses) give the "face-referenced" S; a full PEC sheet gives
`S11 = −e^{−2jk0 z_top}`; a slab-free empty stack gives `S21 = e^{−jk0 (z_top − z_bottom)}` (pure delay, negative phase).
`|S|` and all conservation/reciprocity relations are invariant under this change; the tests confirm that `at_planes` equals a direct
recomputation at the new planes.

## 7. Analytic results used for validation

- **Bare slab (no metal, `J = 0`)**, `r = (1−n)/(1+n)`, `n = √ε_c`, `u = e^{−2jk0 n h}`, faces as planes:
  `S11 = S22 = r(1−u)/(1−r²u)`, `S21 = S12 = (1−r²)e^{−jk0 n h}/(1−r²u)`; half-wave slab `S11 = 0`, `S21 = −1`; lossless
  `|S11|² + |S21|² = 1`.
- **Full PEC sheet**: the total tangential E on `z = 0` vanishes, so the reflection referenced to `z = 0` is `−1`:
  `S11 = −e^{−2jk0 z_top} I₂`, `S21 = S12 = 0`; `S22` is the reflection of a short-circuited substrate line seen from below:
  `S22 = (Z_in − 1)/(Z_in + 1) · e^{−2jk0 d_bottom}` with `Z_in = j tan(k0 n h)/n` (lossless; complex `n` with loss).

## 8. Verification (tests)

| Claim | File |
|---|---|
| Sheet and stack fields (E **and** H) vs independent field matching in all three regions, both sides of the sheet, both ports, TE/TM, oblique, lossy: agreement ~1e-14; source jump condition; continuity at both interfaces; Poynting-flux power balance (incl. dissipation by quadrature); decay/overflow safety; Rayleigh and substrate-cutoff points; mirror symmetry of the ports | `test_green_function_fields.py` |
| Bare slab = Fresnel/Airy closed form (lossless/lossy, 4 plane placements, all four S blocks) and = independent field-matching decomposition; half-wave slab; free space | `test_sparameters.py` |
| Full PEC sheet: `S11 = −e^{−2jk0 z_top}`, `S21 = 0` (any mesh/order), `S22` shorted-line analytic | same |
| Energy conservation `SᴴS = I` (6 lossless cells incl. rotated dipole, slots, spanning wire; 2 plane sets), lossy passivity | same |
| Reciprocity `S = Sᵀ` incl. non-zero cross-pol (rotated dipole/slot), lossless and lossy | same |
| Mirror-symmetric cell ⇒ zero cross-pol; 90° rotation covariance `S_B = R S_A Rᵀ`; x↔y symmetry | same |
| Measured incident amplitudes; radiation-condition residual; `at_planes` = recomputation; `exp(+jωt)`: phase slope, positive group delay | same |
| Mesh convergence (first order, ratio 1.89), Floquet convergence (~M⁻²), vanishing metal → bare slab | same |
| Bottom-port excitation vs field matching and slab transmission; multi-RHS solve | `test_mom_excitation.py`, `test_mom_solution.py` |

Measured: slab vs Airy `~1e-16`; PEC sheet `S11 = −1` and `|S21| ~ 1e-15`; lossless `‖SᴴS − I‖∞ ~ 4e-16`; `‖S − Sᵀ‖∞ ~ 2e-16`
with `|S21[x,y]| ≈ 0.09` (rotated dipole; not vacuous). Mutation checks (scratch copy): flipping the sign of the wave decomposition, of
the plane phases, of the `H` sign, of the bottom-port mirror or current sign, of `ρ`/`r`, swapping ports or planes are each caught.

## 9. Limitations and hand-off (Phase 5 has since been implemented: DATABASE.md)

- Normal incidence, x/y linear basis (`compute_s_parameters`); oblique incidence is Sec. 11 (TE/TM, power-wave normalised).
- Only the specular order; grating lobes are not ports (warning issued).
- The circular basis (`e_R = (x+jy)/√2`, `e_L = (x−jy)/√2`) and its handedness under `exp(+jωt)` are fixed in DATABASE.md Sec. 7; how the holographic synthesis uses them is still open.
- Geometry and frequency sweeps, and the unit-cell database, are implemented in Phase 5 (DATABASE.md); they call `compute_s_parameters` unchanged.
- Validation against an external full-wave solver (HFSS/CST/openEMS) with `exp(−jωt)` conjugation is still open.

## 11. Oblique incidence (Phase 20 R1; additive, `compute_s_parameters_oblique`, `ObliqueSParameters`)

**What already carried the transverse wavenumber.** `LayeredGreenFunction` (TE/TM equivalent impedances at every `kt`), `stack_fields` and
`sheet_fields` (per-mode transmission lines for any `(kx, ky)`), `assemble_impedance_matrix` and `assemble_excitation_vector`
(`kx0, ky0` in the Floquet sums and in `Φ(k0)`) and `RooftopBasis.current_coefficients` were already correct at `kt < k0` (audit:
`tests/test_oblique_audit.py`, against textbook transmission-line impedances, Fresnel/Airy reflection, rooftop quadrature and
`Z(k)ᵀ = Z(−k)` reciprocity). Only `sparameters.py` assumed normal incidence (`0.0, 0.0` into assembly, excitation, currents and the (0,0)
field decomposition) and `config.py` rejected `theta_deg ≠ 0`.

**Basis and reference planes.** `p = (cosφ, sinφ)`, `q = ẑ × p`; mode 0 = TE (`E ∥ q`), mode 1 = TM (tangential `E ∥ p`);
`θ = 0, φ = 0`: TE = y, TM = x. Ports are the same planes as in Sec. 1/6 (`ReferencePlanes`), only the (0,0) order, with `kz0 = k0 cosθ`;
re-referencing multiplies by `exp(∓2j kz0 D)` (`ObliqueSParameters.at_planes`). Time factor `exp(+jωt)` as everywhere.

**Wave amplitudes.** `E_t = V_TM p + V_TE q`, `H_t = −I_TE p + I_TM q` (as `LayeredGreenFunction._apply_modes`); characteristic
admittances `Y_TE = kz0/(ωμ0)`, `Y_TM = ωε0/kz0`; `V_up/down = (V ± I/Y)/2` per mode, measured from the total (stack + sheet) field at
each plane (so the incident amplitude used in `S = b/a` is the measured one). The returned S are **power-wave normalised**,
`S_ij = sqrt(Y_i/Y_j) · V_i/V_j`: co-pol terms equal the tangential-E ratio, a TE↔TM term carries `sqrt(Y_TM/Y_TE) = 1/cosθ`, and for a
lossless cell with only the (0,0) order propagating `Σ|S|² = 1` (grounded: over the one port's two modes; ungrounded: over both ports).
Cross-polarisation is output (full 2×2 per block). A higher propagating Floquet order at `(θ, φ)` raises the same `UserWarning` and sets
`higher_orders_propagate` in the object and in `to_record()` / `sweep_s_parameters_oblique` records.

**Checks.** Bare grounded and ungrounded slab vs the analytic TE/TM transmission-line result at 0°, 30°, 60° (≤ 1e-10); energy
conservation with a patch at oblique angles (≤ 2e-9); `(θ, φ) ↔ (θ, φ+180°)` equality for a centrosymmetric cell; TE/TM at `φ = 0`
equal to those at `φ = 90°` for a square patch on a square lattice; `θ = 0` equals the x/y solver with TE = y, TM = x; and the
bit-identical gate: normal-incidence results equal the Phase 20 base-commit values (`tests/data/r1_theta0_reference.json`).

## 10. References (method; equations here are derived, see Section 0)

1. D. M. Pozar, *Microwave Engineering*, Wiley (transmission-line input impedance, S-parameter and reciprocity/unitarity relations).
2. M. Born and E. Wolf, *Principles of Optics*, Cambridge Univ. Press (Fresnel coefficients, Airy summation for a plane-parallel slab).
3. J. A. Kong, *Electromagnetic Wave Theory*, Wiley, 1986; W. C. Chew, *Waves and Fields in Inhomogeneous Media*, IEEE Press, 1995 (layered-medium TE/TM
   transmission-line fields).
4. K. A. Michalski and J. R. Mosig, IEEE Trans. Antennas Propag., vol. 45, no. 3, pp. 508–519, 1997 (layered-media Green's functions).
5. R. Mittra, C. H. Chan, and T. Cwik, Proc. IEEE, vol. 76, no. 12, pp. 1593–1615, 1988; B. A. Munk, *Frequency Selective Surfaces*, Wiley, 2000
   (reflection/transmission from Floquet-harmonic MoM for periodic surfaces).
