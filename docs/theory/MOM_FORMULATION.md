# MoM discretization and matrix assembly (Phase 3)

Implementation: [src/basis_functions.py](src/basis_functions.py) (mesh, rooftops, transforms),
[src/mom_matrix.py](src/mom_matrix.py) (Z, V), [src/excitation.py](src/excitation.py)
(unpatterned-stack field, reference planes), [src/solver.py](src/solver.py)
(`assemble_system`, `build_system`, `solve`).
Conventions: [CONVENTIONS.md](CONVENTIONS.md). Kernel: [GREEN_FUNCTION.md](GREEN_FUNCTION.md).

**Pipeline:** geometry → mesh / rooftop basis → excitation vector `V` → impedance matrix `Z` → `Z·I = V`.
Not included here: S-parameters (Phase 4, S_PARAMETERS.md), sweeps and the database (Phase 5, DATABASE.md), magnetic-current (aperture) formulation.

## 0. Provenance

As in Phase 2, the equations are **derived here** (Sections 3–6), not copied from a text; the
method (spectral-domain Galerkin MoM with rooftop functions for periodic printed structures) is
standard, see Section 12. I did not have those references available to check line by line. They
are instead verified by independent means (Section 10): exact analytic cases, Parseval identities
with exactly solvable mock kernels, an independent spatial-domain evaluation, and an independent
plane-wave field-matching solver for the excitation.

## 1. Problem and decisions

- One rectangular unit cell `[-Px/2, Px/2] × [-Py/2, Py/2]`, stack air | substrate | air, PEC pattern
  of zero thickness on the air/substrate interface `z = 0`.
- Unknown: the tangential **electric surface current** `J(x, y)` [A/m] on the metal `M ⊂ z = 0`.
- **Slots are metal-complement currents**: for `type: slot` the metal `M` is everything in the cell
  *except* the slot (`UnitCell.is_metal` returns the complement). No magnetic-current formulation exists.
- Basis: x- and y-directed rooftops on a uniform rectangular mesh (Section 5). Galerkin testing.
- Kernel: the existing `LayeredGreenFunction` (`Ẽ_t = G J̃`), used unchanged.

## 2. Notation and conventions (all inherited, none changed)

`exp(+jωt)`, SI. `A = Px·Py`. `k_mn = (kx_m, ky_n)`, `kx_m = kx0 + 2πm/Px`, `ky_n = ky0 + 2πn/Py`,
`(kx0, ky0)` the incident transverse wavevector (0 at normal incidence).

- Floquet series: `J(ρ) = Σ_mn J̃_mn e^{−j k_mn·ρ}`, `J̃_mn = (1/A) ∬_cell J e^{+j k_mn·ρ} dA`.
- Green's dyad: `Ẽ_t,mn = G(k_mn) J̃_mn` (ohm; no `1/A` inside `G`).
- Plane-wave transform of a function on the whole plane: `F(k) = ∬ f(ρ) e^{+j k·ρ} dA`.
- Inner product in testing: `⟨a, b⟩ = ∬ a·b dA`, **no complex conjugation** (Galerkin with the real
  rooftop functions; this is what makes `Z` complex-symmetric, not Hermitian).
- Incident field: unit amplitude `|E_inc| = 1 V/m`, phase referenced to the pattern plane `z = 0`.

## 3. Integral equation

The total tangential field is `E_t^stack + E_t^scat`, where `E_t^stack` is the field the stack would
have at `z = 0` **without** the metal (Section 4) and `E_t^scat` is radiated by `J` in the stack. PEC:

```
E_t^stack(ρ) + E_t^scat(ρ) = 0        for ρ in the metal M          (every cell, by periodicity)

E_t^scat(ρ) = Σ_mn G(k_mn) J̃_mn e^{−j k_mn·ρ}
```

## 4. Excitation: field of the unpatterned stack at z = 0

Use the transmission-line (TL) model of GREEN_FUNCTION.md Sec. 2 for each TE/TM mode
(`p = k_t/|k_t|`, `q = ẑ×p`; TE has `E_t ∥ q`, TM has `E_t ∥ p`). An incident wave travelling
toward −z with mode voltage `V⁺` (its tangential-E amplitude along `q` or `p`) is partly reflected by
the substrate slab. With `a, b` = kz in air / substrate, `ε_c` the complex permittivity,
`u = exp(−2j b h)` (`|u| ≤ 1`):

```
r_TE = (a − b)/(a + b)                        (Z_TE ∝ 1/kz)
r_TM = (B − A)/(B + A),  A = a, B = b/ε_c     (Z_TM ∝ kz/ε)
```

`r` is the reflection coefficient of the air→substrate interface for the TL voltage; the substrate
bottom face has `−r`. Summing the multiple reflections (Airy):

```
Γ = r (1 − u)/(1 − r² u)                  reflection coefficient of the slab seen from air
1 + Γ = (1 + r)(1 − r u)/(1 − r² u)       ratio  V(0)/V⁺  of the total tangential E at z = 0
```

The tangential field at `z = 0` for a unit incident field `E_inc` (3-vector, ⟂ k_inc) with tangential
part `E_t,inc = (E_x, E_y)`:

```
E_t^stack = (1+Γ_TM) (E_t,inc·p) p + (1+Γ_TE) (E_t,inc·q) q
```

At `kt = 0` TE = TM and `p, q` are arbitrary: `E_t^stack = (1+Γ)·E_t,inc`. The incident phase is
`e^{−j k_00·ρ}`, so `E_t^stack(ρ) = E_t^stack e^{−j(kx0 x + ky0 y)}`.

(`stack_tangential_field`; verified to ~1e-16 against an independent 3-D field-matching solve,
including TM and lossy stacks; and against thin-film analytics: half-wave slab `1+Γ = 1`, quarter-wave
slab `1+Γ = 2/(1+ε_r)`.)

## 5. Mesh and rooftop basis

### Mesh
`nx × ny` equal rectangles, `dx = Px/nx`, `dy = Py/ny`, rectangle `(i, j)` centred at
`(−Px/2 + (i+½)dx, −Py/2 + (j+½)dy)`. A rectangle is **metal iff its centroid is metal**
(`is_metal`). Feature edges that lie on mesh lines are represented exactly (choose `nx` such that
`dx` divides the feature dimensions); other edges are staircased.

### Unknowns
One rooftop per mesh edge whose **two** neighbouring rectangles are both metal:

- x-directed on the edge at `x = x_e` between rectangles `(i−1, j)` and `(i, j)` (indices wrap mod `nx`);
  `y_c` is the centre of row `j`.
- y-directed on the edge at `y = y_e` between `(i, j−1)` and `(i, j)` (wrap mod `ny`).

Ordering: all x-directed (row-major in `[i, j]`), then all y-directed.

```
f_x(x, y) = tri((x − x_e)/dx) · pulse((y − y_c)/dy) · x̂
f_y(x, y) = pulse((x − x_c)/dx) · tri((y − y_e)/dy) · ŷ
tri(t) = max(0, 1 − |t|),   pulse(t) = 1 for −½ ≤ t < ½ else 0
```

Peak value 1 at the edge; the coefficient `I_b` is the current density [A/m] there.
Support: two mesh rectangles. The component normal to each mesh line is continuous
(divergence-conforming), the tangential component is a pulse; `∇·f = ∓1/dx` piecewise constant on the
two halves, so each rooftop carries zero net charge.

### PEC boundary condition on `J`
A rooftop exists only between two metal rectangles, so **the current component normal to a
metal/non-metal boundary vanishes identically** for any coefficients (the nearest rooftop's tri function
reaches zero exactly on that boundary). No current exists outside the metal.

### Periodic wrap and connected structures
Edge indices wrap, so the edge at `x = −Px/2` and the one at `x = +Px/2` are **one** unknown (referenced
to the copy at `−Px/2`, likewise `−Py/2` in y). This is essential for wires spanning the cell and for
slot complements (metal touching the cell boundary): the current continues into the neighbouring cell.
Its Bloch extension over the plane is

```
f^B(ρ) = Σ_{p,q ∈ Z} f(ρ − p Px x̂ − q Py ŷ) · e^{−j(kx0 p Px + ky0 q Py)}     ⇒   f^B(ρ + Px x̂) = e^{−j kx0 Px} f^B(ρ).
```

### Closed-form transform
With `np.sinc(z) = sin(πz)/(πz)` and `Φ_b(k) ≡ F[f_b](k)` for the *unwrapped* rooftop centred at `(x_c, y_c)`:

```
Φ_x(kx, ky) = e^{+j(kx x_c + ky y_c)} · dx sinc²(kx dx/2π) · dy sinc(ky dy/2π)
Φ_y(kx, ky) = e^{+j(kx x_c + ky y_c)} · dx sinc(kx dx/2π)   · dy sinc²(ky dy/2π)
```

**Key fact:** for `k = k_mn`, the transform over one cell of the whole Bloch extension equals
`Φ_b(k_mn)` of the single unwrapped rooftop, because the lattice sum of phases `e^{+j k_mn·(pP)}` cancels the
Bloch phases `e^{−j(kx0 pPx + …)}` (`e^{j 2π m p} = 1`). Hence wrapped and interior rooftops are treated
identically and **no explicit wrap phases appear** in Z or V.

Floquet coefficients of a solved current: `J̃_mn = (1/A) Σ_b I_b Φ_b(k_mn)` (`current_coefficients`).

## 6. Galerkin equations

Test the PEC equation of Section 3 with each **unwrapped** rooftop `f_i` (real; its support lies in the
metal of the periodic plane, so `E_tot = 0` there):

```
Σ_j I_j ⟨f_i, E^scat[f_j]⟩ + ⟨f_i, E^stack⟩ = 0.
```

Insert the Floquet sum. For a Bloch-periodic field `E = Σ Ẽ_mn e^{−jk_mn·ρ}`,
`⟨f_i, E⟩ = Σ_mn Φ_i(−k_mn)·Ẽ_mn`, and for real `f`, `Φ_i(−k) = conj(Φ_i(k))` (real k). With
`J̃_mn = (1/A) Σ_j I_j Φ_j(k_mn)` this gives, writing `d_b ∈ {x, y}` for the direction of rooftop `b`:

```
Z I = V

Z_ij = −(1/A) · Σ_{|m|,|n| ≤ M}  conj(Φ_i(k_mn)) · G_{d_i d_j}(k_mn) · Φ_j(k_mn)          [ohm·m²]

V_i  =  conj(Φ_i(k_00)) · (E_t^stack)_{d_i}                                               [V·m]
```

`I` in A/m. The **minus sign** is `G = −Z_eq` (the physical E of a current is `−Z_eq J`), which makes `Z`
an impedance with positive-semidefinite Hermitian part. `M` is `floquet_order`.
Normalization check: for a uniform current (all `J̃_mn = 0` except `(0,0)`) `Z` reduces to
`−dx dy G(0)` per rooftop row, giving `Z_eq = η0/2` for a free-space sheet.

Assembly (`assemble_impedance_matrix`): the four blocks `xx, xy, yx, yy` are computed independently as
`(conj(Φ_a) ⊙ G_ab) Φ_bᵀ` over chunks of harmonics (memory bounded); nothing is symmetrized by hand.
Cost `O(N_b² (2M+1)²)`, `N_b` unknowns.

### Bloch phase reference
The unknown for a seam edge is referenced to the copy inside the cell, so at oblique incidence the
current coefficient carries the Bloch phase of that copy. This is a fixed, consistent gauge; it changes nothing
physical and no explicit phase is needed (see key fact above).

## 7. Structural properties of Z (all tested)

- **Normal incidence:** `Z = Zᵀ` (complex-symmetric, not Hermitian when lossy). Reason: `G` symmetric and even, and the
  harmonic set is symmetric under `k → −k`.
- **Oblique:** `Z(k₀)ᵀ = Z(−k₀)` (Lorentz reciprocity for a Bloch problem); `Z(k₀)` alone is not symmetric.
- **Passivity:** `½(Z + Zᴴ) = Re Z` is positive semidefinite (`½ Re(Iᴴ Z I)` = power delivered to the fields);
  strictly positive definite with loss. Lossless with one propagating order: `rank(Re Z) = 2` (the two polarizations
  of the (0,0) harmonic).
- **Rank rule:** rooftop centres repeat with the mesh, so harmonics `m` and `m + nx` (`n` and `n + ny`) are
  indistinguishable there; the truncated set must contain every residue class:
  **`2·floquet_order + 1 ≥ max(nx, ny)`**, otherwise `Z` is rank-deficient (`MeshSettings` rejects it). This is the
  *minimum*; accuracy needs a larger order.
- **Uniform current:** on a fully metal cell, `Σ_j Z_ij = −dx dy G_xx(0)` for every mesh and order (`sinc²` vanishes at all
  aliases `m = ±nx, ±2nx, …`).

## 8. Analytic case: full-cell PEC sheet

With metal everywhere, `E_tot(z=0) = 0`, so the air side has a PEC mirror, the substrate side is unexcited, and
`J = 2 H_tan,inc`: for normal incidence `J = 2E₀/η0 = 5.3088 mA/m`, **independent of the substrate**.
In the model: `J = E_t^stack/(Z_eq)` with `E_t^stack = (1+Γ)V⁺` and `Z_eq = Z_air ∥ Z_down`,
`1+Γ = 2Z_down/(Z_down+Z_air)`, so `J = 2V⁺/Z_air`. For oblique incidence:
TE `J = 2cosθ/η0` (`Z_TE,air = η0/cosθ`, `V⁺ = E₀`); TM `J = 2/η0` (`V⁺ = E₀cosθ`, `Z_TM,air = η0 cosθ`).
The uniform current is exactly representable, so this is reproduced to `~1e-16` for any mesh, any order satisfying the
rank rule, and any substrate (lossless, lossy, `ε_r` 1–10.2). It validates sign, `1/A` normalization, the excitation/Green
consistency (TE and TM), and the wrap logic. (Oblique TE/TM: the rooftop interpolation of a phase-varying current is only
approximate; observed error `2e-5 → 2e-6 → 3e-7` for `nx = 8, 16, 32`.)

## 9. Discretization behaviour (measured, not proven)

For a 6 × 2 mm rectangle in an 8 × 8 mm cell (RT5880, 10 GHz), aligned meshes `nx = 8, 16, 32`, order `4·nx`:

- Zeroth Floquet coefficient of the solved current `J̃_00,x`: successive differences ratio ≈ 1.9 (≈ first order,
  as expected for rooftops at PEC edges).
- Floquet order at fixed mesh: Cauchy differences shrink by ≈ 4 per doubling (tail `~ M⁻²`, as in GREEN_FUNCTION.md).
- **PEC condition is weak, not pointwise:** the tangential total E on metal falls ≈ 1/N when the mesh is refined
  (interior mean `0.131 → 0.065` of `|E_inc|` for `nx = 16 → 32`; smooth-window average `0.170 → 0.084`), while
  pointwise values in the two end cells of the strip stay `~2|E_inc|` (edge singularity of the true current).
  Do not interpret a pointwise residual there as a bug.
- Cost: `nx = 32` (352 unknowns, `M = 128`) ≈ 0.7 s; dense `Z`, dense LU.

Practical rules: align `nx, ny` with the feature dimensions; use `floquet_order` of order `2–4 × max(nx, ny)`; increase both
and check `J̃_00` convergence.

## 10. Verification map (tests)

| Claim | File |
|---|---|
| Support, peak, linear ramps, continuity of normal component, tangential jump, charge neutrality; periodic (Bloch) continuity across the seam; partition of unity; zero normal current at PEC edges; no unknown across metal/gap; slot = complement; transform vs quadrature; transform of Bloch extension = unwrapped transform | `test_basis_functions.py` |
| Parseval identity with mock kernels (`G = I` ⇒ `Z = −` analytic rooftop Gram matrix; `G = [[0,1],[1,0]]` ⇒ exact x–y overlap), definition-by-brute-force, chunk invariance, uniform-current identity, free-space `η0/2`, symmetry/reciprocity (normal, oblique), PSD, rank 2, independent spatial-domain evaluation, Floquet convergence | `test_mom_matrix.py` |
| Stack field vs independent field matching (TE, TM, normal, oblique, lossy), thin-film analytics, `V` normalization/phase/projection by quadrature, reference planes independent of `Z`, `V` | `test_mom_excitation.py` |
| Full PEC sheet `2/η0` (normal, TE, TM, lossy, any mesh/order), PEC residual convergence, no current in gaps, mesh/Floquet convergence (incl. slot), mirror parity, 90° rotation covariance, `x↔y`, linear-system residual, config route | `test_mom_solution.py` |

Mutation checks (scratch copy, not committed): flipping the sign of `Z`, dropping `1/A`, mis-wiring the `xy` block,
dropping the conjugate in `V`, wrong periodic neighbour, wrong Bloch-phase sign, `sinc²` on the wrong axis, halved rooftop
width, flipped transform phase, wrong TE/TM/Airy excitation terms, and removing the rank rule are each caught.

## 11. Reference planes (explicit; used from Phase 4)

`reference_planes` in `config.yaml` gives the thicknesses of the modelled air regions; the planes sit at their outer
boundaries: `z_top = +top_air_thickness`, `z_bottom = −(h + bottom_air_thickness)` (`h` = substrate thickness).
They do **not** enter `Z` or `V` (the Green's function treats both air regions as semi-infinite; verified bit-identical).
They only define where reflected/transmitted phases are referenced when S-parameters and de-embedding are added.
The incident wave's phase reference stays at `z = 0` (Section 2).

## 12. Limitations and hand-off to Phase 4

- x/y rooftops on a uniform rectangular mesh only; non-aligned or rotated features are staircased (first-order geometry error).
- Dense `Z`; `O(N_b²(2M+1)²)` assembly. For larger problems exploit the two-level Toeplitz/aliasing structure (not done).
- Slots only via metal-complement currents (large unknown count); no aperture (magnetic-current) formulation.
- **S-parameters (done in Phase 4, S_PARAMETERS.md):** `V(z)` outside `z = 0` is now available
  (`sheet_fields`, `stack_fields`); the bottom-port excitation is `stack_tangential_field(..., port="bottom")`.
- Finite metal conductivity and metal thickness: not modelled.

## 13. References (method; equations here are derived, see Section 0)

1. R. F. Harrington, *Field Computation by Moment Methods*, Macmillan, 1968. (Galerkin MoM.)
2. A. W. Glisson and D. R. Wilton, "Simple and efficient numerical methods for problems of electromagnetic radiation
   and scattering from surfaces," *IEEE Trans. Antennas Propag.*, vol. 28, no. 5, pp. 593–603, 1980. (Divergence-conforming
   surface-current basis functions.)
3. S. M. Rao, D. R. Wilton, and A. W. Glisson, "Electromagnetic scattering by surfaces of arbitrary shape,"
   *IEEE Trans. Antennas Propag.*, vol. 30, no. 3, pp. 409–418, 1982. (Basis-function conventions.)
4. J. R. Mosig, "Integral equation technique," in *Numerical Techniques for Microwave and Millimeter-Wave Passive
   Structures*, T. Itoh, Ed., Wiley, 1989. (Spectral-domain MoM for printed structures.)
5. R. Mittra, C. H. Chan, and T. Cwik, "Techniques for analyzing frequency selective surfaces — a review," *Proc. IEEE*,
   vol. 76, no. 12, pp. 1593–1615, 1988. (Floquet-harmonic MoM for periodic surfaces; convergence/relative convergence.)
6. B. A. Munk, *Frequency Selective Surfaces: Theory and Design*, Wiley, 2000.
7. Plus the references of GREEN_FUNCTION.md (layered-medium TE/TM transmission-line formulation).
