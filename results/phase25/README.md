# Phase 25 — Butler-fed sector scanner (model study)

Brief: `docs/briefs/PHASE25_BUTLER_SECTOR_SCANNER.md`. Everything here is **MODEL_ONLY, NOT_HFSS_VALIDATED**.
Numbers are on the held-out 0.5° × 1° grid unless a table says otherwise.

## A. Forward model for an arbitrary excitation (`forward_validation.json`)

Code: `src/hologram/butler_forward.py`; tests: `tests/test_hologram_butler_forward.py`; driver:
`scripts/research/phase25_forward_validation.py` (`python … <out.json>`; `--mom-only` reruns Part II only).

Model: `E(d) = Σ_m a_m · G(d)·FT[ρ̂_m · E_m · Γ_n](d)` — each cell radiates `Γ_state × E_inc`, the feeds' radial
currents add linearly in the complex excitation `a`. First order, local periodicity, **no feed–feed coupling**.

### Part I — pattern path, hex4 (side 150 mm, pointing-corrected, 5698 cells), designed excitations

| Check | Result |
|---|---|
| (a) design's own incident field (`incident='stored'`) vs `sll_board.evaluate` / `_FarField` | complex far field identical (max rel. diff 0.0 on 4000 hemisphere points); every metric identical. This is a plumbing check (same operations by construction). |
| (b) incident field rebuilt analytically at the cell centres (`incident='analytic'`) | directivity ≤ 0.005 dB, pointing 0.000°, SLL ≤ 0.07 dB, XPD ≤ 0.03 dB. Aperture field: median 0.009°, p99 0.3°, max 7.8° (6 cells within 4 mm of a feed), mean stored modulus 0.989 (bilinear interpolation of the unit phasor). |
| (c) keep the ρ^-1/2 surface-wave taper (`amplitude_model='physical'`) | directivity −0.3…−0.5 dB, **global SLL worse by 2.4…3.0 dB single-port (2.8 dB all-ports; worst all-ports y-cut 5.2 dB)**, pointing ≤ 0.14°. The existing synthesis is phase-only; the taper is a real, unvalidated sensitivity. |

Existing numbers reproduced (single port, held-out): 25.84 dBi, pointing 0.10°, global co-pol SLL −13.7 dB, XPD 17.8 dB.

### Part II — what the repo's full-array MoM can check

The only multi-cell MoM in the repo is the **periodic N×N supercell under a normal-incidence plane wave**
(`phase12c.grounded_supercell`, `phase12d.floquet_power`). There is no finite board fed by a surface wave, so the
monopole's incident field (and feed–feed coupling) is `NOT_AVAILABLE_FROM_CURRENT_DATA`. What it does test is the
local-periodicity step "a cell radiates its isolated reflection phasor": for a 5×5 block (D = 1.12 λ, the smallest
supercell with propagating orders besides the specular one) of the 1-bit pair (empty / square), grounded RT/duroid 5880
1.575 mm, 10.54 GHz, the forward model's Floquet-order powers (Fourier sum of isolated-cell phasors from the SAME mesh
and harmonic cut-off, renormalised to the MoM's radiated total) are compared with the full supercell solve.
Known answer (uniform board): 0.0 pp at every mesh.

Worst per-order power-fraction error, percentage points (designed block of port 1 / block of port 3 / seeded random):

| mesh (cells / 6.4 mm) | represented square side, gap | isolated square phase | designed | port-3 block | random | forward model radiates (designed) |
|---|---|---|---|---|---|---|
| 12 | 5.33 mm, 1.07 mm | +75.7° | 2.9 | 10.1 | 4.1 | 0.79 |
| 16 | 5.60 mm, 0.80 mm | +35.9° | 11.9 | 30.4 | 15.3 | 0.54 |
| 20 | 5.76 mm, 0.64 mm | −2.7° | 32.6 | 57.2 | 42.4 | 0.32 |
| 24 | 5.87 mm, 0.53 mm (≈ design 5.85 mm) | −33.9° | 42.1 | 64.5 | 68.6 | 0.24 |

(The rooftop mesh represents a patch as an even number of mesh cells; meshes 6 and 8 fill the whole 6.4 mm cell and are
not used. The empty cell is +138.14°. The isolated square's phase moves ~−0.2°/µm (5.33→5.87 mm: +75.7°→−33.9°), so the pair is strongly resonant.)

**Named state: `LOCAL_PERIODICITY_POOR_AT_DESIGN_GAP_PLANE_WAVE_CHECK`.** At the design gap the isolated-cell model
puts 54–81 % of the incident power into evanescent orders (the MoM cannot: it radiates 99 %), i.e. the cells' in-context
reflection differs strongly from the isolated one. The check is plane-wave, normal incidence, 5×5, one mesh family —
it does not prove the surface-wave-fed board is equally wrong, and it does not show the pattern model's *relative*
beam shapes are wrong; it shows the cell-wise Γ model is not supported at this patch gap by the repo's only full-array
solver. HFSS of the 1-bit pair as one cell (STATUS open step 2) and of a small block is the way to settle it.

---

## Follow-ups (2026-10-13), requested after the stop at A

### 1. Is it a mesh artefact? (`block_locality.json`, `scripts/research/phase25_block_locality.py`)

Small square | empty supercells of the design pair, the 5.867 mm square (= 22 cells at mesh 24 = 44 cells at mesh 48,
gap 0.533 mm), normal-incidence plane wave, full periodic solves (the 1×2 stripe is a 2×2 supercell with identical
rows). Isolated square: −32.4° at mesh 24, −40.6° at mesh 48 (Phase 22's Richardson value is −41.4°); empty +138.1°.
Specular |S11_xx| (all that propagates for a 12.8 mm period), MoM vs the isolated-cell model:

| block | mesh 24: MoM / model | mesh 48: MoM / model | worst Floquet-order complex error (mesh 24 → 48) |
|---|---|---|---|
| 1×2 stripes | 0.999 / 0.083 | 0.999 / 0.011 | 1.06 → 1.04 |
| 2×2 checkerboard | 0.999 / 0.083 | 0.999 / 0.011 | 1.04 → 1.03 |
| 2×2, three squares + one empty | 0.996 / 0.502 | 0.996 / 0.497 | 1.93 → 1.91 |

Halving the mesh moves the isolated phases by ~8° and does not change the disagreement: **not a mesh artefact.**
Caveat on reading it: a period below λ radiates everything into the specular order, so |S11| ≈ 1 by energy
conservation for any pattern, and the isolated-cell model (|mean Γ| < 1, the rest in evanescent orders) is *forced*
to disagree whenever the two states are far apart. These blocks show the isolated model cannot describe the near
field of adjacent states; the 5×5 blocks above (propagating orders) are the more representative test, and the
quantity HFSS can check on the small blocks is the **phase** of S11 (package below).

### 2. Pairs with a minimum gap (`gap_pair_study.json`, `pair_locality.json`)

Phase 22 `auto_states`, 1-bit (empty + the square closest to empty + 180°), 6.4 mm cell, substrate as above:

| min gap | 10.33 GHz | 10.54 GHz | 10.75 GHz | d(phase)/d(side) at 10.54 | mesh uncertainty |
|---|---|---|---|---|---|
| 0.2 mm (current) | 5.90 mm, 178.6° | 5.85 mm, 179.6° | 5.80 mm, 178.5° | 265°/mm | 5.0° |
| 1.0 mm | 5.40 mm, 72.7° | 5.40 mm, 81.0° | 5.40 mm, 90.1° | 126°/mm | 1.3° |
| 1.1 mm | 5.20 mm, 53.5° | 5.20 mm, 59.4° | 5.20 mm, 65.8° | 85°/mm | 1.3° |
| (8.0 mm cell, gap ≥ 1.0) | 6.95 mm, 175.8° | 6.85 mm, 178.7° | 6.80 mm, 177.2° | 174°/mm | 4.0° |
| (7.0 mm cell, gap ≥ 1.0) | 5.87 mm, 89.9° | 5.87 mm, 100.5° | 5.87 mm, 111.8° | 149°/mm | 1.7° |

(separation between the empty cell and the square; the 1.0 mm, 1.1 mm and 7.0 mm rows are `UNREACHABLE_LEVEL`: the
square is capped at P − gap before the phase reaches 180° away.) **At the 6.4 mm period no square with a ≥ 1.0 mm gap gets past 90°.**

Locality error of each candidate in the same plane-wave check (worst Floquet-order power error, percentage points;
designed block / port-3 block / seeded random):

| pair (isolated separation at the checked mesh) | gap | D | error |
|---|---|---|---|
| 6.4 mm, 5.33 mm square (68°) | 1.07 mm | 1.12 λ | 3.7 / 12.2 / 5.1 |
| 6.4 mm, 5.60 mm square (102°) | 0.80 mm | 1.12 λ | 11.9 / 30.4 / 15.3 |
| 6.4 mm, 5.76 mm square (141°) | 0.64 mm | 1.12 λ | 32.6 / 57.2 / 42.4 |
| 6.4 mm, 5.87 mm square (172°, the current pair) | 0.53 mm | 1.12 λ | 42.1 / 64.5 / 68.6 |
| 8.0 mm, 6.86 mm square (172°) | 1.14 mm | 1.12 λ (4×4) | 29.3 / 10.9 / 55.3 |

**The error tracks the separation, not the gap or the period.** A ~180° contrast needs the patch near resonance
(hundreds of degrees per mm), and a resonant patch is detuned by whatever sits next to it. The 8 mm cell fixes the
gap and the separation but not the locality.

**Recommendation: no pair satisfies "separation ≥ ~150° with small locality error" with a resonant square at this
substrate and frequency** (6.4 mm: separation and locality trade directly; 8.0 mm: 172° but 11–55 pp). The two
candidates the package carries are the ends of that trade-off: the 6.4 mm / 5.33 mm square (locality-best: 3.7–12 pp,
≈ 65° separation, ≈ −5.4 dB main-beam amplitude relative to a 180° pair, from `sin(Δ/2)` — arithmetic, not a
simulation) and the 8.0 mm / 6.85 mm square (separation-best, poor locality). Choosing between them, or looking for a
different cell (non-resonant states, other shapes, a thinner or higher-ε substrate), is a design decision that needs
the HFSS check below first. All of this is plane-wave, normal incidence, 4×4/5×5 blocks, model-only.
