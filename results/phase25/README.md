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

---

## 3. The continuous, modest-modulation family (`continuous_family.json`, `scripts/research/phase25_continuous_family.py`)

Question: does the classic route — a continuous square-patch library restricted to a non-resonant size range, with a
smooth sinusoidal-type modulation of moderate depth — give a local-periodicity picture the isolated-cell model can
live with? Same plane-wave 5×5 supercell, mesh 24 (squares are even numbers of 0.267 mm cells), 10.54 GHz,
grounded RT/duroid 5880. Everything is MODEL_ONLY.

### Library range, phase and reactance span

Isolated patch-plane Γ of each even-cell square (mesh 24; mesh 48 agrees to ≤ 3° up to 4.8 mm):

| side | 0 | 2.13 | 3.20 | 3.73 | 4.27 | 4.80 | 5.33 |
|---|---|---|---|---|---|---|---|
| phase | 138.1° | 136.8° | 132.7° | 128.3° | 120.6° | 105.7° | 70.0° |
| normalised sheet reactance X = Im[(1+Γ)/(1−Γ)] | 0.38 | 0.40 | 0.44 | 0.48 | 0.57 | 0.76 | 1.43 |
| local slope | ~1°/mm | | | 12°/mm | 29°/mm | 28→68°/mm | resonance onset |

**Non-resonant range used: empty … 4.80 mm (gap ≥ 1.6 mm): phase 138.1° → 105.7° (span 32°), X 0.38 → 0.76, |Γ| ≥ 0.999.**
Allowing 5.33 mm (gap 1.07 mm, slope ~70–90°/mm, X 1.43) extends the span to 68° and is the edge of the resonance
(the 5.6 mm square is already at 25° and 5.87 mm at −41°). Known answer (uniform 4.80 mm block): 0 error.

### Blocks and depths

Centre phase 112°; depth = peak-to-peak phase swing 15° / 30° / 50° (the 50° case needs the 5.33 mm square).
`x_p1`: m = cos(2πx/D) (modulation along the incident E); `diag_p11`: m = cos(2π(x+y)/D); `hex_crop`: the Fong-type
pattern m = cos(k_sw ρ − k0 sinθ (x cosφ + y sinφ)) of the hex150 port-1 beam (30°, 45°) on a 5×5 crop of the
real design's cells (non-periodic: the supercell has one seam). The exact hex150 neighbourhoods are therefore
reproduced up to that seam; the local gradient range of the real pattern (115–335 rad/m) brackets the two periodic
harmonics (196, 277 rad/m).

| depth | block | worst propagating-order power error | dominant modulated order | MoM \|A\| / model \|A\| | phase MoM − model | specular error |
|---|---|---|---|---|---|---|
| 15° | x_p1 | 2.4 pp | (±1,0) | 0.040 / 0.084 = **0.47** | −1° | 0.008 |
| 15° | diag_p11 | 0.0 (evanescent) | (±1,±1) | 0.035 / 0.084 = 0.41 | +34° | 0.010 |
| 15° | hex_crop | 0.3 pp | (0,1) / (±1,0) | 0.051 / 0.052 = **0.98** / 0.45 | +20° / −5° | 0.008 |
| 30° | x_p1 | 3.9 pp | (±1,0) | 0.055 / 0.111 = 0.49 | −2° | 0.012 |
| 30° | diag_p11 | 0.0 | (±1,±1) | 0.048 / 0.111 = 0.43 | +32° | 0.015 |
| 30° | hex_crop | 0.6 pp | (0,1) / (±1,0) | 0.074 / 0.072 = 1.02 / 0.47 | +18° / −7° | 0.012 |
| 50° | x_p1 | 5.8 pp | (±1,0) | 0.073 / 0.139 = 0.52 | −3° | 0.019 |
| 50° | diag_p11 | 0.0 | (±1,±1) | 0.064 / 0.139 = 0.46 | +30° | 0.025 |
| 50° | hex_crop | 1.8 pp | (0,1) / (±1,0) | 0.148 / 0.152 = 0.97 / 0.44 | +23° / −4° | 0.039 |

Reading it:
* **Power-level locality is no longer a problem**: ≤ 6 pp in every propagating order (the 1-bit pairs: 42–69 pp at the
  current gap, 4–12 pp only at ~65° separation), specular order to ≤ 4 % in amplitude, and it degrades gently with depth.
* **Amplitude is not reproduced where the modulation runs along the incident E**: the isolated-cell model
  over-predicts the first-order amplitude by ≈ 2× (ratio 0.41–0.52, depth-independent), the order along y
  (⊥ E) to 2–3 %. Phase is within a few degrees for the E-parallel order but 20–34° off for the y and diagonal ones.
  The incident field here is a plane wave polarised along x; a surface wave's E is radial, so a board sees a
  position-dependent mixture of the two cases, and a polarisation-averaged correction is not something this check
  can supply.

### Expected beam-amplitude cost versus the ideal 1-bit 180° pair

Ideal sinusoidal phase modulation of peak-to-peak swing Δ puts J1(Δ/2) into each first order; an ideal 0/180°
square wave puts 2/π = 0.637. Cost = 20 log10(J1(Δ/2) / 0.637):

| depth Δ | ideal J1 | cost vs ideal 1-bit 180° | measured in the MoM (E ∥ modulation / E ⊥ modulation) |
|---|---|---|---|
| 15° | 0.065 | −19.8 dB | −24.1 / −21.9 dB |
| 30° | 0.130 | −13.8 dB | −21.3 / −18.7 dB |
| 50° | 0.212 | −9.5 dB | −18.8 / −12.7 dB |

For comparison the 65° locality-best 1-bit pair gives 2 sin(Δ/2)/π ⇒ −5.1 dB (arithmetic). A multiplexed hologram
that carries four beams shares the same total swing, costing roughly another 6–12 dB per beam. Not captured here
and not obtainable from this repo: the radiated *fraction* of the surface-wave power. A weak modulation leaks little
per wavelength, and a 150 mm aperture is only ~5 surface wavelengths, so a modest-depth design may leave most of the
guided power unradiated (the model's `aperture_field` has unit amplitude everywhere and no leakage term). Absolute
gain and efficiency of this family are `NOT_AVAILABLE_FROM_CURRENT_DATA`.

### Recommendation for B/C

**Use the continuous modest-modulation family if B/C is run at all, and only as a relative study; do not use
1-bit pairs.** It is the only family here for which the cell-wise forward model has support (power errors ≤ 6 pp,
specular ≤ 4 %). Conditions:
1. Treat the forward model's modulated-order amplitudes as upper bounds along E (×~0.5 measured), report beam
   *shapes*, spacing, crosstalk and steering trajectories, and label every absolute gain/efficiency number
   `NOT_AVAILABLE_FROM_CURRENT_DATA`.
2. Keep the library to the non-resonant range (≤ 4.80 mm; ≤ 5.33 mm only at the deepest depth) and re-run this
   plane-wave check on the actual synthesised board's neighbourhoods before quoting anything.
3. The honest outcome may be "the sector scanner is not worth building with a passive modest-modulation board":
   multiplexing four beams at 10–25 dB below a 1-bit pair, with unknown leakage, will probably be radiating a
   small fraction of the feed power. The HFSS locality package (current pair first) and one HFSS run of a
   modest-depth modulated block would settle whether the continuous family's 0.5× amplitude factor is real.
