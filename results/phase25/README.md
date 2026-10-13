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
