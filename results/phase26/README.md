# Phase 26 — leakage-controlled tensor-impedance CP hexagon (model study)

Brief: `docs/briefs/PHASE26_LEAKY_TENSOR_CP_HEXAGON.md`. Everything on the model side is **MODEL_ONLY,
NOT_HFSS_VALIDATED**; the HFSS column is the user's run of `hfss_exports/four_port_hex150_gap0`
(`hfss_results/hex150_gap0/`, port 3 alone, 10.5 GHz).

## A. Leakage model, re-prediction of the existing hex150 port 3 (`hex150_leakage.json`)

Code: `src/dispersion/leakage_model.py` (closed form + Floquet harmonic balance), `src/hologram/leaky_forward.py`
(ray accounting, illumination taper), tests `tests/test_dispersion_leakage_model.py`,
`tests/test_hologram_leaky_forward.py`; drivers `scripts/research/phase26_hex150_leakage.py` and
`scripts/research/phase26_context_contrast.py`. Additive: no existing module changed, the Phase 23 rectangle
snapshot gate is bit-identical.

### The model
A TM surface wave of constant β on the grounded slab, loaded by a two-state patch sheet, leaks into the n = −1
Floquet harmonic that the hologram points at the beam. First order in the modulation (sheet-modulation /
penetrable-sheet definition; the slab is the non-local background):

`α = ½ |c₁ η₀κ/k₀|² [cos²ζ Re Y_up,TM/|Y_TM|² + sin²ζ Re Y_up,TE/|Y_TE|²] / (2 w)`, `κ = √(β²−k₀²)`,
`c₁ = |Y₁−Y₀| sin(πd)/π` (duty `d` = the board's own patch fraction, 50.5 %), `ζ = φ_beam − ψ` the in-plane angle
between the ray and the radiated harmonic, `w = |∂B/∂β| (η₀κ/k₀)²/4` the surface-wave power per unit width and
|H₀|² (air + slab). Checked against an independent truncated Floquet harmonic balance (TM+TE, every harmonic
with its own slab admittance): **≤ 0.5 % for weak modulation, 3 % at ΔY/Y ≈ 0.4**, α ∝ M² to 3 %.
The converted fraction of the launched surface-wave power along a ray is `1 − exp(−2 α ρ_edge(ψ))`, the
incident field of every cell is multiplied by `ρ^−1/2 exp(−α ρ)`.

### What the model says (port 3, 1414 cells, 714 patches, states re-solved at 10.5 GHz: empty +138.3°, 5.85 mm square −39.0°)

| Case (named state) | contrast ratio κ | α (n=−1) | converted η | D_aperture | gain, referred to accepted power, **launch fraction = 1** |
|---|---|---|---|---|---|
| `PHASE23_CELLWISE_UNIT_AMPLITUDE` (the existing model: all power through a phase-only aperture) | — | — | 100 % (assumed) | 25.7 dBi | 25.7 dBi |
| `LEAKY_ISOLATED_CELL_CONTRAST` | 1.00 | 0.93–1.20 Np/m | 25.1 % | 25.3 dBi | **19.3 dBi** |
| `LEAKY_MOM_CONTEXT_CONTRAST_LOW` | 0.86 (1×2 stripes, 2-cell period) | 0.51–0.62 | 14.4 % | 25.3 dBi | **16.9 dBi** |
| `LEAKY_MOM_CONTEXT_CONTRAST_HIGH` | 1.04 (3-cell-period stripes) | 1.10–1.43 | 29.1 % | 25.3 dBi | **19.9 dBi** |
| HFSS port 3, realized gain RHCP (peak) | | | | | **19.85 dBi** (directivity 19.97, total efficiency 97.7 %) |

(κ = in-context ÷ isolated contrast of the sheet admittance, from the plane-wave supercell MoM: `context_contrast.json`
and `results/phase25/block_locality.json`. 1×2 stripes 0.87/0.86 (mesh 24/48), 2×2 checkerboard 0.97/0.97, 3-cell
stripes 1.03 (duty 1/3), 1.04 (duty 2/3); 4-cell stripes 1.52 are excluded — the 25.6 mm period is 0.9 λ, not a
homogenisable supercell. The hologram's own periods are 18.6–54 mm.)

**The result is a bracket, not a calibration.** The leakage rate scales as κ², so the MoM evidence (κ 0.86–1.04)
spans 16.9–19.9 dBi; ±0.1 in κ is ±1.5–2 dB (`gain_upper_bound_vs_contrast_ratio`). HFSS lies inside the bracket,
at the high end. Three things the model does not settle:

1. **Launch fraction.** The monopole's share of power entering the surface wave (the rest radiates directly) is
   not computed — there is no pin-launching model in this repo — so every gain above is an UPPER BOUND:
   0.75 → 18.0, 0.5 → 16.3, 0.35 → 14.7 dBi at κ = 1. The HFSS pattern has a roughly azimuth-uniform component at
   θ 60–90° (0.022–0.049 of the power per 45° sector, vs 0.2 in the beam sector at 20–70°), which is what a
   directly radiating λ/4 pin would look like; this is a hypothesis the data do not exclude.
2. **First-order accuracy.** The Floquet balance (total over all radiating harmonics, stripe assumption) gives
   α between 0.58 and 1.87 Np/m against the closed form's 0.93–1.20 at this depth (ΔY/|Y| ≈ 1): ±40–55 % in
   α, ±2 dB in η. Its total converted fraction is 29.3 % (κ = 1).
3. **Where the power goes.** The model puts η·(aperture-pattern fraction) = 14 % within 8° and 22 % within 20° of
   the beam; HFSS has 19.4 % and 34.7 % (RHCP 17.0 % / 20.8 %). The model's aperture pattern is narrower in power
   terms; its peak is right only because the two errors partly cancel. The ≈ 75 % unconverted remainder is
   assumed to radiate outside the beam — it is not modelled (no edge diffraction).

HFSS power budget (`port3_realized_gain_RHCP_LHCP_10p5GHz.csv`, 2° grid, integral 1.05 × accepted power, ±5 %):
≤ 8° 19.4 %, ≤ 20° 34.7 %, θ > 70° 38.0 %, back hemisphere 14.4 %, LHCP 48 % of the radiated power (14 % of the total
within 20° of the beam, 15 % at 20–70° outside that cone); XPD at the RHCP peak 22.3 dB (the model's phase-only aperture: 17.8–18.9 dB).

Gate result: `LEAKAGE_MODEL_BRACKETS_HFSS_NOT_DISCRIMINATING` (central contrast −0.6 dB, pessimistic edge −3.0 dB,
launch fraction unmodelled). The brief's stop criterion (> 2 dB miss) is neither cleanly met nor cleanly avoided, so
Phase 26 stops here for the user's decision before item B.

## A2. Launch fraction of the actual pin (`launch_fraction.json`)

Code: `src/dispersion/pin_launch.py`, tests `tests/test_dispersion_pin_launch.py`, driver
`scripts/research/phase26_launch_fraction.py`. The pin is a vertical current filament through the grounded
slab (the same TM transmission line along z as `LayeredGreenFunction`), optionally ending in a top-hat disk;
the radiated power splits into the TM0 pole residue and the `kt <= k0` space-wave integral.

Checks: no-slab short dipole over ground 0.2 % from `eta0 k0^2 (IL)^2/(6 pi)`; pin + top-hat + ground images
against an independent free-space far-field integral 0.01–0.1 %; TM0 pole = `tm0_beta` to 1e-9; lossy-slab
numerical total vs `P_space + P_sw` within 1 %; thin-slab limit launches < 0.1 % and recovers the no-slab
power; launch fraction grows with slab thickness; converged in segments (12 → 48) and angles to < 1 %.

| Reading of "PinLength 110 mil, radius 12 mil" | top above the surface | launch fraction (TM0 share of radiated power) | R_rad | X (thin-wire TL est.) |
|---|---|---|---|---|
| from the ground plane (`mono_h` convention) | 1.219 mm | **0.592** | 1.7 Ω | −115 Ω |
| above the substrate surface | 2.794 mm | **0.600** | 6.5 Ω | −67 Ω |

The launch fraction is **insensitive to pin height** (0.587 at −0.8 mm, 0.61 at 4.4 mm, 0.69 at 10 mm): a thin
slab (h = 0.055 λ) puts roughly 60 % of a vertical dipole's radiated power into TM0 whatever the pin.
Consistency flag: either reading has R_rad of a few ohms and a large capacitive X, i.e. |S11| of −0.1…−3 dB, while
the HFSS S11 of port 3 is −16.4 dB and the −10 dB band 8.5–11.4 GHz. The model's thin-wire estimate resonates
at 6.4 mm total pin length (R ≈ 23–26 Ω). So the feed in the HFSS project is electrically not a 110 mil bare
pin in this model's sense (a longer effective pin, a gap/annulus, a feed structure) — the launch fraction of
that feed is therefore still an estimate, and `hfss_exports/launcher_calibration/` measures it.

### Item-A bracket re-run with the launch fraction (launch 0.592 / 0.600)

| Case | gain (launch 1, A) | gain with launch | vs HFSS 19.85 dBi |
|---|---|---|---|
| leaky, isolated contrast κ = 1 | 19.3 | **17.0 / 17.1** | −2.8 / −2.8 dB |
| leaky, MoM context low κ = 0.86 | 16.9 | **14.6 / 14.7** | −5.2 / −5.2 dB |
| leaky, MoM context high κ = 1.04 | 19.9 | **17.6 / 17.7** | −2.2 / −2.2 dB |

With a physical launch fraction the model **under**-predicts HFSS by 2.2–5.2 dB. The launch fraction explains
about 2.3 dB of the 5.9 dB between the old phase-only 25.7 dBi and HFSS, the leakage conversion (25 %) about
6.0 dB, and the remaining ≈ 2.5 dB is power that reaches the HFSS beam and the model does not account for. The
HFSS needs a converted fraction of about 0.48 of the launched TM0 power against the model's 0.25 (κ = 1), i.e.
a leakage rate ≈ 2.3× larger (equivalently a sheet-contrast ratio ≈ 1.5; the 4-cell-period supercell MoM gave
1.52 but that cell is not a valid homogenisation). Named state:
`LEAKAGE_PLUS_LAUNCH_MODEL_UNDERPREDICTS_HFSS`. Candidate causes, none excluded by the data: the pin's direct
radiation grazing over the patterned aperture and scattering into the beam (a reflectarray-like path the
surface-wave picture omits); a feed whose TM0 share is not 0.6; the first-order leakage rate (±40 % at this
depth); the hologram's non-stripe local structure.

## A3. Launcher as a design variable (`launcher_sweep.json`)

`scripts/research/phase26_launcher_sweep.py`: 10.54 GHz, pin top −1.2…+10 mm about the substrate surface
(negative = fully inside the substrate), top-hat radius 0–5 mm, plus the thin-wire resonant-tuned heights.
Match is the thin-wire estimate (`MATCH_ESTIMATE_THIN_WIRE_TL`, tens of ohms), only a guide.

| Option | launch fraction | est. |S11| | launched share of accepted power | recovery vs current pin (launch fraction) |
|---|---|---|---|---|
| current pin (1.2 mm above the surface) | 0.594 | −0.1 dB | 0.013 | — |
| bare pin, resonant: 5.0 mm above the surface | 0.617 | −9.6 dB | **0.549** (best) | **+0.17 dB** |
| 10 mm pin, 4 mm top-hat | **0.821** | −1.3 dB | 0.208 | +1.4 dB, only if matched by other means |
| pin fully inside the substrate | 0.587 | worse | — | −0.05 dB |

**The launcher is not the lever**: the best matched option recovers 0.2 dB and even an unmatched ideal one
1.4 dB, against a 2.2–2.8 dB model deficit. A top-hat helps the fraction (0.69 → 0.82 for a 10 mm pin) but
its match is what the model cannot judge; pin fully inside the substrate is no better than the current pin.
Note the share column is large for the 5 mm pin because the *model's* estimate of the current pin is
unmatched while the HFSS one is −16 dB: it is a model-internal ranking, not a prediction of HFSS gain.

## A4. HFSS calibration package

`hfss_exports/launcher_calibration/` (A: current pin, B: 5 mm resonant pin, C: tall pin + top-hat): the bare
grounded board, one feed at port 3's position, a radial line 1 mm above the surface (rho 10–140 mm, phi
225°). README lists the exports requested and the extraction (launch fraction from the E_z decay).

## B. Tensor cell family, tensor sheet extraction, database coverage (**stop point**)

Code: `src/dispersion/tensor_sheet.py` (extraction + tensor transverse resonance),
`leakage_model.leakage_alpha_first_order_tensor`, tests `tests/test_dispersion_tensor_sheet.py`, drivers
`scripts/research/phase26_tensor_database.py`, `phase26_tensor_coverage.py`; data `tensor_database.json`
(140 rows: 20 (L, W) families × rotations 0–180° in 15° steps, mesh 64), `tensor_coverage.json`,
`tensor_coverage.png`. MODEL_ONLY, `TENSOR_SHEET_FROM_NORMAL_INCIDENCE`.

**Extraction.** The cell is a rotated rectangular PEC patch on the 6.4 mm cell (`Geometry.rotation`, grounded
RT/duroid 5880, 10.54 GHz). The MoM's full 2 × 2 normal-incidence reflection at the patch plane gives
`Y_in = Y0 (I − Γ)(I + Γ)⁻¹` and the patch layer's tensor `Y_s = Y_in − Y_slab I`. The dispersion along azimuth φ
is the root of `(Y_kk + A(β))(Y_ss + B(β)) − Y_ks² = 0` (TM / TE background admittances); for a sheet diagonal in
the propagation frame it IS `sheet_model`'s equation, and the code path is the same (bit-identical at φ = 0).
`sheet_model.py` and `mom_scaffold.py` are untouched.

**Validation.**

| Check | Result |
|---|---|
| square patch → isotropic tensor | off-diagonal and anisotropy < 1e-9 (also for a 45°-rotated square) |
| scalar limit vs `sheet_model` | β bit-identical (φ = 0), 1e-12 rotated frame; tensor of a square = the scalar extraction to 1e-9 |
| reciprocity (symmetric tensor) | max `|Y_xy − Y_yx|/|Y_iso|` = 9e-16 over all 140 entries |
| quarter turn ≡ swapped length/width | identical to 1e-9 (mask-exact) |
| principal axis of the rotated cell | max error 1.8°, rms 0.6° over all rotations and families |
| rotation equivariance `Y(a) = R(a) Y(0) Rᵀ(a)` | median 1.7 %, max 7.8 % — **not** exact: the neighbour gap of a rotated rectangle in a square lattice depends on the angle, so the isotropic part rises up to ≈ 7 % towards 45°; the database holds the actual tensors rather than assuming the rotation law |
| staircase scatter of the rotated mask (mesh 32 / 48 / 64 / 96) | δ and Y_iso ±2–3 % at mesh ≥ 48 (up to 8 % at mesh 32); mesh 64 used |
| mixed-polarisation dispersion | principal-axis β are the scalar β of Y_xx / Y_yy; at 45° to the axes the mode is a TM/TE hybrid (TE content 5 % for 4.0 × 2.4 mm) and β lies between |

**The non-resonant window.** Phase 25's locality criterion (rotated bounding box ≤ 4.8 mm, gap ≥ 1.6 mm)
restricts rotation: a rectangle rotates through every angle only if `√(L² + W²) ≤ 4.8 mm` (12 of the 20
families, e.g. 4.0 × 2.4, 4.0 × 1.6, 3.2 × 1.6). Longer patches only exist at 0° / 90°.

**Coverage (what item C can ask the database for).** A constant-anisotropy, rotating-axis tensor surface
`Y = Ȳ I + δ R(α) diag(1, −1) Rᵀ(α)`, `α = ψ/2` (Fong 2010 / Minatti 2016, impedance definition: sheet
admittance of the penetrable patch layer on the grounded slab) needs one (L, W) per `(Ȳ, δ)`; α is then free.

| full-rotation family (L × W, mm) | Ȳ_iso [mS] | δ [mS] | Ȳ, δ spread over rotation | first-order CP conversion over 0.14 m |
|---|---|---|---|---|
| 4.0 × 2.4 | 0.92 | 0.38 | 2.6 %, 3.2 % | 0.07 % |
| 4.0 × 1.6 | 0.60 | 0.38 | 4.8 %, 3.4 % | 0.06 % |
| 4.0 × 0.8 | 0.34 | 0.30 | 2.6 %, 3.1 % | 0.03 % |
| 3.2 × 1.6 | 0.36 | 0.18 | 3.2 %, 3.4 % | 0.01 % |
| axis-aligned only: 4.8 × 3.2 | 1.94 | 0.72 | — | 0.4 % |

**The map is empty where it matters.** The CP leakage constant of a rotating-axis tensor sheet is
`α ≈ 0.003 Np/m` for the best fully rotatable family (0.015 for 4.8 × 3.2, axis-aligned) — against 1.06 Np/m
for the 1-bit pair of the existing board (which already converts only 25 %) and 8 Np/m needed for 90 % over
0.14 m. Converting 50 % / 90 % at the sheet loading of those families needs δ ≈ 13 mS / 23 mS; the
non-resonant window offers δ ≤ 0.4 mS (0.8 mS axis-aligned). Relaxing the window towards resonance
(axis-aligned, gap 0.4 mm: 6.0 × 3.2 mm) reaches δ = 2.6 mS and 14 %; the first-order formula (checked ≤ 3 %
against the Floquet balance in item A at this depth) is more accurate for these smaller modulations.
Leakage rises steeply with the MEAN sheet loading (`Ȳ_iso`): at δ = 2.5 mS, 2.0 % conversion at Ȳ = 0, 14 % at
4 mS, 38 % at 6 mS, 54 % at 8 mS (β/k₀ 1.018 → 1.067), because the wave is then more tightly bound.

Consequences / options (a decision for you before item C):
1. Allow loaded, near-resonant elements (mean sheet 4–8 mS, δ ≈ 2–3 mS): the conversion is reachable but this
   is the regime where Phase 25 found the isolated-cell/locality assumption failing; each element would need
   the supercell MoM treatment (and the extraction here is normal-incidence, isotropic-in-angle).
2. Change the board: a thicker or higher-ε substrate binds TM0 tighter at the same patch loading (the model is
   ready: only `Substrate` changes). This is the cheapest lever; not applied here.
3. Accept a low-conversion CP design (≤ 1 % non-resonant): the gain would sit 20–30 dB below the aperture
   directivity (conversion 0.1–1 %), i.e. not useful.
4. Longer / bigger aperture than 150 mm hexagon sides: conversion ∝ 1 − exp(−2αρ); still needs α ≳ 1 Np/m.

**Known gap (noted, not built).** The launch model assumes a bare slab; on the real board the pin is
surrounded by the patterned surface. Estimate with the sheet model's loaded β around the feed
(`launch_fraction.json: loaded_slab_mean_sheet_mS`): at the hologram's mean sheet (≈ 4 mS) β/k₀ rises
1.018 → 1.033 and the launch fraction 0.592 → 0.664 (0.600 → 0.691 for the other pin reading), i.e. +0.5 dB;
8 mS (strongly loaded cells) 0.74 / 0.80. The keep-out cell under the feed and the pin impedance change are not
modelled; `PinLaunch(..., y_sheet_s=...)` is the additive hook.

## B2.1 Required map (`required_map.json`, `required_map.png`)

`scripts/research/phase26_required_map.py`: for a mean sheet `Y_bar = j y mS` (→ β/k₀) and a rotating-axis tensor
anisotropy δ, first-order CP leakage (`leakage_alpha_first_order_tensor`, beam at θ = 30°) integrated over the
150 mm hexagon's azimuths; launch fraction from `pin_launch` with the loaded slab (current pin, both PinLength
readings); aperture directivity of the illuminated hexagon (continuous phase, ρ^-½ e^(−αρ), cos θ₀) — an
estimate. Target: converted fraction ≥ 85 %.

| mean loading y [mS] | β/k₀ | δ required [mS] | modulation depth δ/\|Y_down\| | launch fraction (from ground / above surface) |
|---|---|---|---|---|
| 4 | 1.034 | 9 | 0.89 | 0.665 / 0.693 |
| 5 | 1.040 | 7 | 0.75 | 0.681 / 0.717 |
| 6 | 1.047 | 6 | 0.70 | 0.696 / 0.742 |
| 7 | 1.056 | 5 | 0.64 | 0.714 / 0.771 |
| **8** | **1.067** | **4** (min) | 0.56 | 0.738 / 0.805 |
| 9 | 1.079 | 5 | 0.77 | 0.770 / 0.844 |
| 10 | 1.093 | 5 | 0.83 | 0.807 / 0.882 |
| 12 | 1.123 | 7 | 1.35 | 0.877 / 0.943 |

**Target region: mean loading 5–11 mS (β/k₀ 1.04–1.11), anisotropy δ ≥ 4–7 mS** (minimum 4 mS at 8 mS /
β/k₀ 1.067; δ steps of 1 mS). Beyond ≈ 9 mS the radiated harmonic's own loading suppresses the leakage, so δ
required rises again. The gain estimate inside the region is 24–26 dBi (e.g. 25.3 dBi at 8 / 4 mS: D_ap 27.3 dBi,
η 0.85, launch 0.74) — at the edge of the 24 dBi target; the launch fraction improves with mean loading
(0.59 bare → 0.74–0.88). Caveat: the modulation depth δ/|Y_down| is 0.55–0.9 in the region, beyond the 0.4 where
the first-order formula was checked to 3 % in item A — the true region will move; the Floquet balance is the
check once cells exist. Compared with the non-resonant database (δ ≤ 0.4 mS, Y ≤ 2.6 mS) the target needs
**10× the anisotropy at 3–4× the mean loading**.

## B2.2–B2.3 Slotted-circle family and Bloch extraction at the operating point (**stop point**)

Code: `src/dispersion/slotted_cells.py` (`SlottedCircularPatch`, `SlottedSquarePatch`, `CellRow`; duck-typed
`Geometry`), `tensor_sheet.bloch_sheet_tensor / bloch_beta / bloch_frequencies`, tests
`tests/test_dispersion_bloch_sheet.py`; drivers `phase26_slotted_family.py`, `phase26_ramp_locality.py`,
`phase26_eigenmode_package.py`; data `slotted_family.json` (mesh 48, 19 cells), `slotted_family_mesh64.json`,
`ramp_locality*.json`. MODEL_ONLY.

**Extraction.** The cell is driven by an impressed unit tangential field in the (0,0) Floquet harmonic at the
Bloch wavenumber `(kx0, ky0)` (any `kt`, including `kt > k0`); the cell-mean current and the total field
`E_tot = E_stack + G_00 J_00` (`G_00` contains the slab pole) give `Y_s(kt) = J E_tot^-1`, and the dispersion zero is
`det Im(Y_s(kt) + A pp^T + B qq^T) = 0`. At `kt = 0` it equals the normal-incidence tensor to 1e-8 (tested); the
solver needed no change for `kt > k0` — the old eigen-search problem of `mom_scaffold` does not arise because
the slab pole sits in `G_00`, not in the solved matrix. Spatial dispersion matters: at `kt = 1.05 k0` even a
bare circle has `Y_longitudinal != Y_transverse` (2.5 vs 1.8 mS for D 5.0 mm).

**How far is the normal-incidence sheet model off for these cells?** `beta(normal) - beta(Bloch)` [k0], mesh 48:
along the slot ≤ 0.005 (negligible); at 45° / 90° to the slot: D 5.6, 4.6 × 1.2: −0.001 / −0.003; D 5.6,
5.0 × 1.6: +0.02 / −0.02; D 5.8, 5.0 × 1.6: −0.09 / −0.02; D 5.8, 5.0 × 2.0: −0.13 / −0.04; D 6.0, 5.0 × 2.0:
−0.14 / −0.05; D 6.0, 5.4 × 1.6: −0.20 / −0.08; slotted squares up to −0.30 / −0.13. A bare circle: ≤ 0.003. **The normal-incidence model under-estimates the loading
of exactly the cells that are in the target region, by up to ~15 % of β.**

**Mesh sensitivity (the dominant numerical caveat).** A slot that ends 0.2–0.3 mm from the circle edge leaves a
copper bridge of about one mesh cell: Y_yy of D 5.8, 5.0 × 2.0 is 12.4 (mesh 48), 11.1 (64), and at mesh 32 the
bridge is lost entirely (Y_xx 0.7 vs 13.6 mS). Bridges ≥ 0.4 mm are stable (D 6.0, 4.8 × 2.0: Y_yy 10.6 / 11.0 at
mesh 48 / 64). The joint cell-and-wave rotation check of the circular outline shows 5–7 % staircase scatter
(exact at 90°).

**Coverage vs the required map** (operating point kt = 1.05 k0; Ȳ = (Y_∥ + Y_⊥)/2 across/along the slot,
δ = half the difference; required δ interpolated from `required_map.json`; first-order conversion over the
hexagon):

| Cell (circle D, slot L × W) | bridge | mesh | Ȳ [mS] | δ [mS] | δ required | conversion est. |
|---|---|---|---|---|---|---|
| 6.0, 4.8 × 1.6 | 0.5 mm | 64 | 8.1 | 2.7 | 4.5 | 0.58 |
| 5.8, 4.8 × 2.0 | 0.3 mm | 64 | 7.3 | 2.8 | 4.5 | 0.60 |
| 6.0, 4.8 × 2.0 | 0.4 mm | 64 | 8.8 | 3.3 | 5.0 | 0.72 |
| 6.0, 5.0 × 1.6 | 0.39 mm* | 64 / 48 | 9.0 / 9.1 | 3.6 / 4.0 | 5.0 | 0.76 / 0.82 |
| 5.8, 5.0 × 2.0 | 0.22 mm* | 64 / 48 | 8.6 / 9.4 | 4.1 / 5.0 | 5.0 | 0.86 / 0.92 |
| 6.0, 5.0 × 2.0 | 0.33 mm* | 48 | 10.6 | 5.5 | 5.6 | 0.88 |
| 6.0, 5.4 × 1.6 | 0.19 mm* | 48 | 12.8 | 7.8 | 7.8 | 0.89 |
| square 5.4, 5.0 × 1.6 | — | 48 | 14.6 | 9.3 | 10.2 | 0.84 (β 1.04 → 1.32) |

(* numerically fragile: mesh 48 → 64 moves Y by 10–25 %, and a 0.2–0.3 mm bridge is a fabrication risk.)
**Result:** with mesh-stable bridges the family reaches 65–85 % of the required anisotropy (conversion ≈ 0.6–0.75
first order, short of 0.85); the 0.85+ cells are the thin-bridge ones. Slotted squares reach the δ but with a
β contrast of 1.04 → 1.32 across the slot (one cell goes past resonance, no bound root at 90°): not usable.
Operating-point caveat: the leakage coupling uses Y at the radiated harmonic's kt ≈ 0.5 k0 as well; the
formulas here use the surface-wave kt only (a refinement for item C).

**Large directional dispersion.** For the in-region cells β/k₀ is 1.04 along the slot and 1.10–1.20 across
it: the local phase velocity of the surface wave varies by 8–15 % with the slot angle. The reference-wave
phase of a hologram must be solved with the local tensor dispersion (eikonal), not a fixed β — and the
launch region (feed at the centre, slot angle varying around it) is strongly anisotropic.

## B2.4 Locality of continuous modulation (`ramp_locality*.json`)

Supercell MoM (one periodic solve, Bloch `kx0 = β`, own Floquet orders) vs the adiabatic sum of the local
cells' β along the ramp; wave along the ramp; slot-angle ramps are exactly periodic (0 → 180°). Gate: ≤ 15°.

| Ramp (mesh 24 unless noted) | local β/k₀ range | super β/k₀ | mean local | accumulated phase error | gate |
|---|---|---|---|---|---|
| angle 0→180°, 8 cells, D 5.6, 4.4 × 1.2 | 1.031–1.040 | 1.0363 | 1.0365 | −0.1° | PASS |
| angle 0→180°, 8 cells, D 6.0, 4.6 × 1.6 | 1.041–1.105 | 1.0640 | 1.0754 | **−7.4°** | PASS |
| angle 0→180°, 8 cells, **D 6.0, 4.8 × 2.0** (target-like) | 1.041–1.170 | **1.0593** | 1.1286 | **−44.9°** | **FAIL** |
| size 5.2→5.8→5.2, 10 cells | 1.028–1.038 | 1.0321 | 1.0322 | −0.1° | PASS |
| size 5.6→6.0→5.6, 10 cells | 1.031–1.041 | 1.0372 | 1.0373 | −0.1° | PASS |
| size 5.2→6.1→5.2, 10 cells | 1.025–1.044 | 1.0330 | 1.0333 | −0.3° | PASS |

(Second-start roots agree to 1e-8; D 6.2 mm has no bound root at mesh 24 — the 0.2 mm neighbour gap is below
one mesh cell. A ramp of the in-region D 5.8, 5.0 × 2.0 cell at mesh 32 is invalid: its slot bridge is not
resolved there, see above.)

**Reading.** Weak-to-moderate modulation (size ramps, Δβ ≤ 0.07 k₀) is local to < 8°. A strongly anisotropic
rotating slot — the regime that meets the conversion target (β varying 1.04 → 1.17 over the cycle) — is **not
local**: the supercell's β is 0.07 k₀ (6 %) below the average of the local β's, −45° accumulated over 8 cells
(51 mm ≈ 1.9 guided wavelengths: the slot rotates through 180° inside two guided wavelengths, which is the
rate the hologram needs). An effective-medium predictor — the area-averaged local tensor — does much better
(β 1.085 vs 1.059, +17° — just outside the gate) and is the natural "local model" for item C, but it is not
adiabatic ray optics. Phase 25's abrupt 1-bit failure is not repeated for gentle modulation; it is for the
modulation depth the conversion target requires.

## B2.5 HFSS eigenmode ground truth (`hfss_exports/slotted_cell_eigenmode/`)

Five master/slave unit-cell scripts (D 5.8 / 6.0 mm × slot 0° / 45° + bare D 5.8), two Bloch settings each
(β = 1.04 k₀ along x: 84.2° per cell; on the diagonal 59.5° + 59.5°); predictions written before any run
(`predictions.json`, README table): the bound TM-like mode at 10.17–10.57 GHz (bare circle 10.566 / 10.556 GHz;
D 5.8 slot 0°: 10.574 along x, 10.27–10.32 on the diagonal; slot 45° swaps them; D 6.0: 10.535 / 10.19 and
10.21 / 10.515), a TE-like neighbour 0.04–0.3 GHz above on the slotted cells; mesh 48 vs 64 agree to ≤ 0.05 GHz.
The scripts have not been run in AEDT (boundary property names may need a manual assignment; README).

## B2.4b Corrected locality gate: Floquet balance of the modulated tensor sheet (`ramp_floquet.json`, `database_crosscheck.json`) — **gate FAILS for the target cell**

Labels: `SOLVER_INTERNAL`, `MODEL_ONLY`, `NOT_HFSS_VALIDATED` (the HFSS eigenmode package was not run, by decision).

**Why B2.4's reference was wrong.** The Floquet constant of a modulated surface is not the mean of the local
constants. The reference is now the harmonic balance of the staircase sheet built from the LOCAL cell tensors
(`src/dispersion/tensor_floquet.py`, the tensor twin of `leakage_model.floquet_complex_beta`, reusing its
harmonic blocks unmodified; unit tests: uniform-tensor limit = `find_beta_tensor`, scalar two-level staircase =
the scalar balance to 1e-8, Fourier series reconstruct the profile). Cell tensors are the Bloch tensors
`Y_s(kt)` at one wavenumber for every harmonic. Gate: |β_Floquet − β_supercell| ≤ 0.01 k₀.

| Ramp (8 cells, slot 0→157.5°, mesh 24) | β_super/k₀ | local mean | effective medium | **Floquet (tensors at β_super)** | Floquet (own β_i) | gate |
|---|---|---|---|---|---|---|
| D 5.6, 4.4×1.2 (weak) | 1.0363 | 1.0365 | 1.0359 | 1.0363 (err −0.0000) | 1.0363 | PASS |
| D 6.0, 4.6×1.6 (moderate) | 1.0640 | 1.0754 | 1.0643 | 1.0667 (err +0.0028) | 1.0669 | PASS |
| **D 6.0, 4.8×2.0 (primary / target)** | **1.0593** | 1.1286 | 1.0850 | **1.0763 (err +0.0170 k₀)** | 1.0773 | **FAIL** |

* The Floquet balance is a real improvement on the local mean (error 0.069 k₀ → 0.017 k₀) and on the
  effective medium (0.026 → 0.017), and is within 0.003 k₀ for moderate anisotropy. For the primary cell it
  still over-predicts β by 0.017 k₀ — over 8 cells (51.2 mm) that is **+11° accumulated phase**, inside the old
  15° phase gate but outside the 0.01 k₀ gate set for this step. n_max 6 vs 10 agree to 1e-4 k₀; tensors at β_super
  vs each cell's own β_i differ by 0.001 k₀ (not the cause).
* **α agreement** (supercell α is one Newton step of the complex MoM tensor determinant — radiation into n = −1
  and dielectric loss, an estimate, not a converged complex root): target 4.1 Np/m (Floquet) vs 2.8 (supercell);
  moderate 2.2 vs 2.7; weak 0.07 vs 0.09. Same ordering, ±50 % at the target — not an agreement to build on.
* Not tested here, the likely residual: (i) the same tensor for every harmonic although the n = ±1 harmonics
  see kt = β ± K (spatial dispersion, large for this cell family); (ii) a cell's tensor changes when its
  neighbours differ (mutual coupling is in the supercell MoM, not in a per-cell tensor).

## B2.4c Database Bloch β vs uniform-supercell solve (`database_crosscheck.json`) — PASS (`SOLVER_INTERNAL`)

D 6.0, 4.8×2.0, slot 0° / 45° / 90°, propagation 0° and 45°, 4 identical cells in one MoM (orders 4m × m):
|Δβ| ≤ 1e-9 k₀ in all six cases (gate 0.005 k₀). This is a self-consistency check (a uniform supercell must
reproduce the cell: Floquet-order truncation, staircase masks, root selection), **not** physical validation, and
it does not test the modulated case. Along the diagonal slot 0° and slot 90° coincide (1.1082 k₀) as the mirror
symmetry requires; β/k₀ for slot 0° is 1.041 along x, 1.138 for slot 90° along x.

## Stop condition
Step 1 of the B2.4b instruction fails for the primary cell (0.017 k₀ > 0.01 k₀); per the instruction, item C is
not started. Not yet done: step 3 (mesh-stable restriction recomputation, which only restates the primary cell).

## B2.4d Harmonic-dependent tensors (`ramp_floquet_dispersive*.json`) — gate still FAILS for the primary cell

Coupling of harmonics n <- q uses the cell Bloch tensors at kt = beta_super + (n+q)K/2 (`floquet_complex_beta_dispersive`,
n_max 3; n_max 4 identical to 1e-4). Same three 8-cell slot-angle ramps, mesh 24:

| Ramp | beta_super/k0 | constant-tensor Floquet | **harmonic-dependent** | error [k0] | alpha [Np/m] (super est.) | gate 0.01 k0 |
|---|---|---|---|---|---|---|
| D 6.0, 4.8x2.0 (primary) | 1.0593 | 1.0769 | **1.0720** | **+0.0127** | 2.58 (2.77) | **FAIL** |
| D 6.0, 4.6x1.6 | 1.0640 | 1.0668 | 1.0633 | -0.0006 | 2.26 (2.68) | PASS |
| D 5.6, 4.4x1.2 | 1.0363 | 1.0363 | 1.0363 | -0.0001 | 0.08 (0.09) | PASS |

The dispersion between harmonics removes about a quarter of the primary-cell residual (0.017 -> 0.0127 k0, and
makes alpha agree to 7 %) but not enough: the rest is cell-environment dependence (a cell's tensor changes with
its neighbours), which a per-cell database cannot carry. Per the instruction, the Floquet-of-local-tensors model is
**not** used for item C; the SOLVER_SUPERCELL table below replaces it.

## B2.4e SOLVER_SUPERCELL table (`supercell_table_*.json`)

Slot-angle-rotating row (slot angle 180 i/N deg, i = 0..N-1; modulation vector parallel to propagation,
gamma = 0), circle D 6.0 mm, mesh 24/cell, propagation period N P, K = 2 pi / (N P). The complex root
beta - j alpha is found by analytic continuation of f(beta) = det(Y_s(beta) + diag(A, B)) (7 real-axis MoM
points, degree-6 polynomial; the degree-5 root agrees to <= 6e-5 k0 in all 20 cells and the polynomial's own
value at Re(root) matches the MoM to 1e-7 relative). The real-axis root of Im det differs from Re(root) by up
to 0.01 k0 where leakage is strong — it is not the complex root. alpha_leak = alpha - alpha of the uniform
row of the slot-0 cell (dielectric loss; 0.03-0.2 Np/m).

| slot L x W [mm] | uniform beta/k0 | N=3 (K 1.48 k0) beta / alpha_leak | N=4 (1.11) | N=5 (0.89) | N=6 (0.74) | N=8 (0.56) |
|---|---|---|---|---|---|---|
| 4.8 x 2.0 | 1.0305 | 1.0679 / 4.20 | 1.0637 / 2.40 | 1.0660 / 2.97 | 1.0728 / 2.08 | 1.0686 / 2.80 |
| 4.8 x 1.6 | 1.0405 | 1.0626 / 4.26 | 1.0592 / 2.51 | 1.0608 / 2.82 | 1.0611 / 3.00 | 1.0632 / 3.15 |
| 4.6 x 1.6 | 1.0405 | 1.0620 / 2.99 | 1.0592 / 2.51 | 1.0580 / 1.94 | 1.0608 / 2.67 | 1.0628 / 2.69 |
| 4.4 x 1.2 | 1.0407 | 1.0538 / 0.67 | 1.0509 / 0.47 | 1.0538 / 0.79 | 1.0522 / 0.58 | 1.0516 / 0.42 |

alpha_leak in Np/m. alpha is not monotone in period (a resonance of the row's periodicity with the cells' own
dispersion: dips at N = 4 and 5-6) and depends on depth only weakly above 4.6 x 1.6; the table is not smooth
enough to interpolate across periods with confidence — nearest-grid values are what it supports.
(4.8 x 1.6 and 4.6 x 1.6 coincide at N=4: 2.51; checked, same root.)

## B2.4f Reachable conversion (`table_reach.json`) — TARGET NOT REACHABLE with this family on this board

Best case: every ray of the 150 mm hexagon (beam 30 deg / 225 deg, beta 1.065 k0, K = 0.565-1.565 k0) uses the
largest alpha_leak of the table at its local K. alpha_available (mean) = 3.4 Np/m; 85 % conversion on every ray
needs 7.0 Np/m. **Best-case conversion 0.60 (rays 0.48-0.72), against the 0.85 target** — before the launch
fraction (0.66-0.88 with the loaded slab, 0.59 bare) and before the oblique-modulation (gamma != 0) approximation.
Gain estimate then: aperture directivity ~26 dBi + 10 log(0.60 x ~0.7) ~ 22 dBi, above the HFSS 19.85 dBi of the
1-bit board but well below the 24 dBi requirement and without the >= 85 % conversion. Mesh-stable cells cannot
carry more modulation on 1.575 mm RT5880 (bridge >= 0.35 mm).

## C. Pure-RHCP rotating-slot tensor hologram, single 150 mm hexagon (`tensor_cp_hexagon.json`, `tensor_cp_hexagon_cells.json`)

Labels: SOLVER_SUPERCELL, MODEL_ONLY, NOT_HFSS_VALIDATED, GAMMA_0_COLLINEAR_ONLY (the table has modulation parallel to
propagation only; the hexagon's rays have modulation vectors at angles up to ~60 deg to the wave: |K| is used,
the obliquity is an approximation, not a result). Module `src/hologram/tensor_leaky_cp.py` (tests:
`tests/test_hologram_tensor_leaky_cp.py`: energy conservation of the march, table-capped conversion, decomposition into
an in-phase circular component and an evanescent one), script `scripts/research/phase26_tensor_cp_hexagon.py`.

**Design.** 1430 cells on the 6.4 mm lattice (D 6.0 circle; feed cell and r < 6.5 mm empty), hexagon circumradius 150 mm
(vertex at +/-y, port-3 orientation), feed at the centre, beam theta 30 / phi 225, 10.54 GHz. A radial march (720 rays,
0.25 mm) takes, at every point, the table column of the NEAREST modulation wavenumber K = |beta rho_hat - k_t|
(grid 0.556, 0.741, 0.889, 1.111, 1.481 k0 — no interpolation across periods), the slot (4.4x1.2 / 4.6x1.6 / 4.8x1.6 /
4.8x2.0, all bridges >= 0.43 mm) whose alpha_leak is nearest the required one, accumulates the table's beta into the
reference phase and decays the power per radian by exp(-2 alpha rho). Slot axis a = (chi + phi_ray)/2,
chi = phi_ref - k_t.r (current direction `chi = 2a - phi_ray`: the rotating-axis sheet reflects the surface-wave E
about the slot axis); `J = sqrt(2 alpha W/rho) exp(-j phi_ref) (cos chi, sin chi)`. hand = -1 produces an LHCP beam
(RHCP fraction 0.025), hand = +1 the RHCP one. The far field is computed from the 1430 currents (grounded-slab dyad), so
purity is a result of the sum, not an input.

| eta_design | taper c | conversion | D_ap RHCP [dBi] | gain (launch 0.592) [dBi] | XPD at peak [dB] | AR at peak [dB] | RHCP SLL [dB] |
|---|---|---|---|---|---|---|---|
| 0.5 | 0.0 | 0.418 | 28.82 | 22.75 | 17.0 | 2.47 | -19.8 |
| 0.5 | 0.5 | 0.464 | 28.69 | 23.08 | 16.6 | 2.57 | -22.8 |
| 0.5 | 0.8 | 0.464 | 28.21 | 22.60 | 17.0 | 2.47 | -23.1 |
| 0.7 | 0.0 | 0.472 | 28.73 | 23.19 | 16.7 | 2.56 | -22.0 |
| 0.7 | 0.5 | 0.519 | 28.58 | 23.45 | 17.0 | 2.47 | -22.9 |
| 0.7 | 0.8 | 0.524 | 28.24 | 23.16 | 16.9 | 2.51 | -23.0 |
| 0.9 | 0.0 | 0.500 | 28.64 | 23.36 | 16.7 | 2.55 | -22.8 |
| **0.9** | **0.5** | **0.538** | **28.48** | **23.51** | **16.9** | **2.50** | **-22.8** |
| 0.9 | 0.8 | 0.548 | 28.28 | 23.39 | 16.7 | 2.54 | -21.9 |

(selected: highest gain with RHCP fraction > 0.5 and SLL <= -13 dB; conversion by ray 0.44-0.64.)

**Selected design vs the brief's targets** (aperture pattern; gain referred to the accepted power = D_ap x launch x conversion):

| Quantity | Predicted | Target | |
|---|---|---|---|
| Beam direction | 30 / 224 deg (grid 2 deg) | 30 / 225 | ok |
| Peak RHCP gain, launch 0.592 (bare slab) | **23.5 dBi** (loaded-slab launch 0.66-0.88: 24.0-25.2) | >= 24 | **short by 0.5 dB (bare)** |
| Aperture efficiency (gain / (4 pi A cos30 / lambda^2) = 23.5 / 29.0 dBi) | 0.28 (0.32 launch x conversion; D_ap alone 0.88) | >= 35 % | not met |
| Axial ratio at peak / max over the -3 dB beam | 2.50 / 3.07 dB | <= 3 dB | marginal: 0.07 dB over |
| Peak XPD | **16.9 dB** | >= 25 | **not met** |
| RHCP fraction of the aperture pattern; LHCP | 0.975 ; 0.025 | >= 0.9 in beam, LHCP <= 10 % | ok (in the aperture pattern) |
| RHCP fraction of power within 8 deg | 0.928 / 0.948 | >= 90 % | ok |
| Co-pol SLL | -22.8 dB | <= -15 | ok |
| Surface-wave power converted in the hexagon | **0.54** (0.44-0.64 by ray) | >= 0.85 | **not met** (table-limited, B2.4f) |
| Grazing + back power | <= 0.68 of accepted power, NOT MODELLED (pin direct 0.41 + unconverted 0.27) | <= 15 % | not met / unknown |
| Frequencies 10.33 / 10.75 GHz | NOT_AVAILABLE_FROM_CURRENT_DATA | | table at 10.54 only |

* **XPD is limited by the grounded slab's element factor, not by the hologram.** At theta = 30 deg the TM and TE
  fields of a circular current differ by ~2.5 dB (AR 2.5 dB = XPD 16.9 dB). Pre-compensating the current
  ellipticity at the beam direction (J' = g(k_t)^-1 J, an UPPER BOUND not realised by this cell family: it needs a second,
  independently controlled modulation) gives AR 1.25 dB, XPD 22.9 dB, RHCP fraction 0.989.
* The aperture model is idealised (phase-perfect, no feed-region scattering, no edge diffraction, no cell-to-cell
  phase error): D_ap 28.5 dBi is 85 % of the geometric limit. The earlier board's model predicted 25.8 and HFSS gave
  19.85 dBi; this gain is therefore an UPPER-BIAS number for the aperture part, while the leakage/launch parts are
  calibrated only against that one board (they under-predicted it).
* Budget of the accepted power (bare 0.592 launch): pin direct space wave 0.41 (outside the beam by assumption),
  unconverted at the edge 0.27, aperture-radiated 0.32, in the 8 deg beam 0.30 (RHCP 0.30), LHCP total 0.008.
  The pin's direct radiation is the largest single loss and is unaffected by the cells.

## S. Substrate study (`results/phase26/substrate/`: `prelim.json`, `sizing.json`, `table_*.json`, `verify_*_mesh32.json`, `reach.json`)

Labels: SOLVER_SUPERCELL, MODEL_ONLY, NOT_HFSS_VALIDATED, GAMMA_0_COLLINEAR_ONLY. Script `scripts/research/phase26_substrate_study.py`
(`size`, `prelim`, `table`, `reach`). No design was made on any new board.

**Method (every choice stated).** Same slotted-circle family. Pitch P = 2 pi / (3 (beta_target + k0/2)) so the shortest
modulation period of the hexagon is three cells; circle D = the LARGEST D whose uniform deep-slot cell still has a bound root
with beta <= beta_target + 0.06 k0 (beta_target = bare + 0.012 k0, as the 5880 1.575 control), gap >= 0.35 mm — the control
also sits at its largest feasible D, just below the resonance edge. For 5880 3.175 the first pick (D 5.25 mm, far from the
edge, alpha_leak <= 0.7) was superseded by D 5.65 (gap 0.35; `_superseded/`). Slots keep the 5880 family's L/D, W/D (deep,
mid). Table points N = 3, 4, 5, 6 (deep and mid) and 8 (deep) at mesh 24, complex root by analytic continuation (as B2.4e).
Quality gates for a table point: alpha_leak > 0, deg-6 vs deg-5 roots agree to 1e-3 k0, complex root within 0.03 k0 of the
real-axis root. Reachable conversion: best case, nearest-K grid of the QC-passing points, alpha constant along each ray,
`1 - exp(-2 alpha rho_edge)` averaged over the hexagon; gain estimate = 28.5 dBi (idealised aperture directivity of item C)
+ 10 log(launch x conversion). TE1 cutoff of the grounded slab = c / (4 h sqrt(eps - 1)): **no board is excluded** (all > 10.54 GHz).

| Board | TE1 cutoff | P / D [mm] | smallest feature [mm] | launch (pin, bare slab) | table points passing QC | alpha_leak [Np/m] | reachable conversion (best case) | gain estimate [dBi] |
|---|---|---|---|---|---|---|---|---|
| RT5880 1.575 (control, existing table) | 43.4 GHz | 6.4 / 6.0 | 0.43 bridge | 0.59 | 20/20 | 0.4-4.3 (max per K 2.5-4.3) | **0.60** (C design: 0.54) | 24.0 (C: 23.5) |
| **RT5880 3.175** | 21.5 GHz | 6.0 / 5.65 | 0.35 | **0.82** | 9/9 | 0.6-1.1 | 0.23 | 21.3 |
| RO3006 1.28 | 25.8 GHz | 6.1 / 4.275 | 0.32 | 0.71 | 3/9 | 0.08-6.9 (K 0.93: 6.9, 0.78: 4.1, 1.55: 0.08) | 0.49 — K gaps, **not reproduced at mesh 32** | 23.9 — unverified |
| RO3006 1.52 | 21.7 GHz | 6.0 / 4.05 | 0.28 | 0.79 | 3/8 (K 0.95 and 1.58 only) | 12.6-17.6 (mesh 32: 14.4-14.6) | 0.98 if that alpha held over all K — **not shown** | 27 — NOT CREDIBLE as an estimate |
| RO3010 1.27 | 19.5 GHz | 6.1 / 3.5 | 0.25 | 0.78 | 4/8 | 0.9-1.5 (N=8: 6.3, not reproduced at mesh 32) | 0.35 | 22.8 |
| RT6010 1.90 | 13.0 GHz (margin 1.23x; loaded-slab TE cutoff is lower: not checked) | 5.3 / 2.125 | 0.15 | 0.96 | 8/8 | 0.4-0.6 | 0.12 | 19.1 |

(smallest feature = min of circle gap, slot width, slot-end bridge. Mesh-stable features need >= 0.35 mm: only the two RT5880
boards qualify. The launch column is the thin-wire pin model at the slab surface; pin top flush vs 1.22 mm above differ by < 0.01.)

**Reading.**
* **RT5880 3.175 mm (first candidate): launch rises to 0.82, but the conversion falls to 0.23.** The thicker slab binds
  the wave more loosely and the cells, at the largest D the 0.35 mm gap allows, are further from resonance: alpha_leak
  0.6-1.1 Np/m against 2.5-4.3 on the 1.575 mm board. Net launch x conversion 0.19 against 0.35 — worse. The table is clean (mesh-
  insensitive roots, 9/9 QC); a larger pitch/diameter (>= 6.4 / 6.05) is not tried and could change this.
* **The high-eps boards only reach large alpha where the cell is almost resonant** (RO3006: beta jumps by 0.2-0.25 k0 from the
  uniform row, alpha 4-17 Np/m). Two of those points (RO3006 1.52, K 0.95 and 1.58) reproduce at mesh 32 to within 20 % in alpha;
  the RO3006 1.28 and RO3010 points do not (the root is lost). Features are 0.15-0.32 mm (mesh- and etch-fragile) and the
  uniform-cell beta itself moves 0.05 k0 between mesh 24 and 32 for RO3006 1.52. A credible estimate would need the
  column of K = 0.6-1.6 filled with reproducible points, which this sizing rule does not deliver.
* RT6010 1.9 mm: weak (the sizing rule drives D to 2.1 mm to hold beta; features 0.15 mm).
* **Nothing reaches the 0.85 conversion of the brief on a mesh-stable family**; the only large-alpha candidate (RO3006 1.52)
  is fragile and incompletely mapped. Honest recommendation: RT5880 1.575 remains the best *verified* board
  (0.60 / 24 dBi best case); a re-sized RO3006 1.52 family (larger bridges, K grid filled) is the one worth a second study.

## H. HFSS of the single hexagon (`hfss_results/hex150_tensor_cp_single/`, `hex150_tensor_cp_single_hfss_comparison.json/.md`)

Feed field of the user's note was left unfilled ("[scripted pin / my connector]"): recorded as `[user to state]`. The
antenna-parameter table is at **10.5 GHz** (design 10.54); the realized-gain CSV is the 10.5 GHz, 1-degree export. The
E-field line (surface-wave decay) was not delivered: `NOT_PROVIDED`.

| Quantity | HFSS (10.5 GHz) | Prediction (10.54) | State |
|---|---|---|---|
| Peak RHCP realized gain | **23.03 dBi** (table 23.05; directivity 23.16 dBi, total efficiency 0.975) | 23.5 (bare launch 0.592) .. 25.2 (loaded) | GAIN_WITHIN_PREDICTED_RANGE (-0.5 dB vs bare) |
| Peak direction | (30, 225) | (30, 224) | OK |
| XPD at the peak | **23.0 dB** | 16.9 | better by 6 dB; target 25 not met |
| AR at the peak / worst over the -3 dB beam | 1.23 / 1.58 dB (cut phi 224: <= 3 dB over 9.5-11.5 GHz at theta 30) | 2.50 / 3.07 | better; target met |
| RHCP share of the power within 8 deg | 0.994 | 0.979 | ok |
| RHCP SLL (outside 15 deg) | -18.5 dB | -22.8 | worse by 4 dB; target <= -15 met |
| S11 | -15.0 dB at 10.54; min -27.3 dB at 10.08; <= -10 dB over 8.81-11.59 GHz | NO_PREDICTION | ok |

Power budget (fractions of the radiated power; radiation efficiency 1.004, so also of the accepted power):

| Region | HFSS | Model (of accepted power, bare launch) |
|---|---|---|
| Main beam <= 8 deg (RHCP) | **0.271** | 0.302 |
| 8-20 deg | 0.024 | in the aperture pattern (0.03) |
| Far sidelobes theta <= 70 | 0.171 | in the aperture pattern |
| Grazing 70-90 deg | 0.264 | lumped with the back in the upper bound 0.68 (pin direct 0.41 + unconverted 0.27) |
| Back hemisphere | 0.269 | |
| **LHCP, all angles** | **0.310** | 0.008 (aperture only; the pin's direct radiation is not given a polarisation) |

* The beam, its purity and its AR are better than predicted; the gain is 0.5 dB under the bare-launch value, as the
  launch fraction (0.59-0.62) plus the 0.54 conversion say. The non-beam power is 0.70 (far sidelobes + grazing + back)
  against the model's 0.68 + aperture sidelobes: **the model's loss accounting is confirmed in total**; the split between
  the pin's direct radiation (0.38-0.41) and the unconverted edge power (0.27) is not separable from these exports (needs the E-field
  decay line).
* **0.31 of the power is LHCP**, essentially all outside the beam (RHCP is 0.994 of the power within 8 deg): the pin's direct
  space wave (a vertical monopole is not circularly polarised) and the edge diffraction. This, not the hologram, is why
  the integrated LHCP misses the brief's <= 10 %.

## I. Launcher redesign on the loaded surface (`launcher_redesign.json`, `scripts/research/phase26_launcher_redesign.py`) -- STOP, no export

Labels: MODEL_ONLY, HFSS_CALIBRATED_FEED_R_SCALE, MATCH_ESTIMATE_THIN_WIRE_TL, NOT_HFSS_VALIDATED.

**Which pin did the HFSS run have?** The S11 sweep (min -27 dB at 10.10 GHz, -15 dB at 10.54, band 8.8-11.6 GHz) is reproduced
by the thin-wire model only for a pin ~5.5 mm above the surface (7.1 mm from the ground, 0.635 mm radius: the OLD exporter default) and
an input resistance 1.64x the model's radiated-power resistance (rms 2.5 dB over 8-12 GHz, 6.6 dB without the factor; model
-13.8 dB at 10.54 vs HFSS -15.1). A 110 mil (2.794 mm) pin, as written in the README of the export, would be
capacitive with R of a few ohms and could not give -15 dB. **Please confirm the pin height used**; the redesign below
uses 7.1 mm / 0.635 mm as the baseline.

**Result.** Direct space-wave share = 1 - launch fraction (TM0 share of the pin's radiated power). Uniform shunt sheet at the
surface = the patterned surface around the feed; pin top tuned for X = 0 (series resonance) or fixed:

| surface loading at the feed [mS] | baseline pin (5.5 mm above the surface) launch / est. S11 | resonant pin (4.8 mm) launch / S11 | top-hat 1.5 mm (resonant 2.2 mm) launch / S11 | top-hat 3 mm (resonant 1.2 mm) launch / S11 |
|---|---|---|---|---|
| 0 (bare slab: the current keep-out, r < 6.5 mm) | 0.62 / -13.8 dB | 0.62 / -17 dB | 0.60 / -9.0 | 0.59 / -4.9 |
| 4 (about the hologram mean sheet) | 0.74 / -12.6 | 0.72 / -31 | 0.70 / -13.4 | 0.68 / -7.7 |
| 8 (the primary cell, D 6.0 4.8x2.0) | 0.89 / -11.2 | 0.87 / -22.5 | 0.83 / -23.6 | 0.79 / -13.3 |
| 12 (beyond the family, near-resonant cells) | 0.98 / -10.9 | 0.98 / -19.3 | 0.97 / -54 | 0.93 / -19.1 |

* **The lever is the surface around the pin, not the pin.** The pin height and a top-hat do not reduce the direct radiation at a
  matched condition (a hat lowers the launch fraction by 0.02-0.08 once the pin is retuned); loading the surface right up to the
  pin raises it from 0.62 to 0.74-0.89, i.e. the pin's direct radiation falls from 38 % to 11-26 % of the accepted power.
  With the HFSS gain as anchor (23.05 dBi at launch 0.62): +0.7 dB (4 mS) to +1.5 dB (8 mS) -> 23.8-24.6 dBi.
* S11 <= -10 dB is kept in every row with loading (the estimate is good to a few dB; HFSS decides). A resonant pin of 4.8 mm is the
  better-matched choice at every loading (-17...-31 dB estimated) with essentially the same launch fraction.
* **Caveats.** (1) The model's sheet is uniform and infinite; the real situation is a finite pad around the pin followed by the
  hologram, so these are bounds: bare slab (0.62, which the HFSS gain is consistent with) and infinite sheet; the radius needed
  for the effect is not modelled. (2) The loaded surface binds TM0 more tightly near the pad; a mismatch to the hologram's
  beta 1.06 k0 at the pad edge is not modelled. (3) The Phase 26 C design leaves a bare disk r < 6.5 mm around the feed: the
  cells inside it (about 8 cells) and a feed-cell keep-out of r < 3.2 mm are the candidate change.
* Recommended HFSS variants (not exported): V0 current design with the pin 4.8 mm above the surface; V1 as V0 with the keep-out
  reduced to the feed cell (r < 3.2 mm, cells filled with the hologram's states as computed); V2 as V1 with the cells of
  r < 30 mm replaced by an unslotted 6.0 mm circle (the high-loading "launch pad", ~12-16 mS).

### I2. Refit with the reported pin (6.0 mm above the surface = 7.575 mm from the ground, 0.635 mm)

The thin-wire model with this pin resonates at 9.4-9.6 GHz against HFSS 10.10 GHz (5 % low); free input-resistance factor 1.35, rms
5.4 dB over 8-12 GHz (6.5 dB unscaled): worse than the 5.5 mm / 0.635 mm fit (1.66, 2.4 dB). With a larger EFFECTIVE radius (0.305 / 0.635 / 1.0 / 1.5 mm:
rms 6.7 / 5.4 / 4.4 / 3.4 dB, resonance 9.2 / 9.4 / 9.6 / 9.8 GHz) the fit improves: the feed (connector body, coax annulus) behaves like a
thicker pin. **The resistance factor is still needed: 1.30 (1.5 mm effective radius; 1.35 for 0.635 mm).** Calibration used from here
on (HFSS_CALIBRATED, one board): top 6.0 mm, effective radius 1.5 mm, R x 1.30: model -13.5 dB at 10.54 (HFSS -15.1), baseline launch 0.627.
Table with that calibration (`launcher_redesign.json`): launch 0.63 / 0.74 / 0.90 / 0.98 at 0 / 4 / 8 / 12 mS of uniform loading; resonant pin
4.77 mm: 0.62 / 0.72 / 0.87 / 0.98 with est. S11 -12 / -17 / -28 / -42 dB. The 110 mil (2.794 mm) reading is dropped from the export READMEs.

## J. V2: graded launch pad (`v2_pad_design.json`, `v2_pad_loading.json`, `hfss_exports/hex150_tensor_cp_single_v2_pad/`)

Mean loading (Bloch tensor, kt 1.06-1.10 k0, mesh 24): unslotted circle D 6.0 **4.4 mS**, D 6.1 4.7 mS, D 6.2 no bound root; slotted cells
(angle-averaged) 4.4x1.2 5.4, 4.6x1.6 6.5, 4.8x1.6 6.7, 4.8x2.0 7.8 mS; cross-slotted D 6.0 L 4.8: W 0.8 / 1.2 / 1.6 / 2.0 -> 6.4 / 6.8 / 7.3 / 7.6 mS.
The highest isotropic loading at stable bridges is 7.6 mS, only 0.7 mS above the hologram's local mean (6.9 mS): V2 mainly fills the bare keep-out.
A new additive geometry `CrossSlottedCircularPatch` (tested: 4-fold symmetric, isotropic tensor). Export: pad r < 30 mm, ~68 cross cells, rest of
the hologram unchanged. Predictions (gain 22.7 / 23.5 / 24.1 dBi for launch unchanged / central / infinite-sheet; grazing + back 0.56 / 0.53 / 0.50; LHCP 0.33 / 0.31 / 0.29;
S11 -13.5 ... -11.5 dB est.) are in the export README. A plain 12-16 mS "unslotted" pad proposed earlier is withdrawn (4.4 mS only).

## K. Aperture-size study, model only (`scaling_study.json`, `scripts/research/phase26_scaling_study.py`)

Single hexagon, side (= circumradius) 150/200/250/300 mm, same cell family (nearest-grid SOLVER_SUPERCELL table, mesh-stable slots),
pin, board and centre feed, theta 30 / phi 225. Design march re-run per size (eta_design 0.7/0.9/0.99 x taper 0.5/0.8; best gain with SLL <= -13 dB);
the far field of the cell currents is summed (5692 cells at 300 mm). HFSS anchor: 150 mm = 23.03 dBi at launch 0.62 -> offset -0.68 dB to the
model, main-beam factor 0.857; grazing+back and LHCP shares = the HFSS split of the non-beam power (0.731, 0.425), an ASSUMPTION until the E-field
line separates pin-direct from edge power. (A first pass truncated the lattice at +/-192 mm for the 250/300 mm sizes; fixed and re-run.)

| Side [mm] | cells | conversion (rays) | D_ap RHCP [dBi] | **gain [dBi] (anchored)** | main <= 8 deg | grazing+back | LHCP | XPD / AR peak, worst [dB] | SLL [dB] | 4-hex board [mm] |
|---|---|---|---|---|---|---|---|---|---|---|
| 150 (HFSS) | 1430 | 0.54 (0.44-0.64) | 28.5 | **23.0** | 0.27 | 0.53 | 0.31 | 16.9 / 2.5, 3.1 | -22.8 | 650 x 525 |
| 200 | 2560 | 0.63 (0.53-0.73) | 30.9 | **26.1** | 0.32 | 0.50 | 0.29 | 16.6 / 2.6, 3.0 | -23.7 | 866 x 700 |
| 250 | 3936 | 0.69 (0.59-0.78) | 32.7 | **28.4** | 0.36 | 0.47 | 0.27 | 17.0 / 2.5, 2.8 | -25.9 | 1083 x 875 |
| **300** | 5692 | 0.74 (0.65-0.82) | 34.2 | **30.2** | 0.38 | 0.45 | 0.26 | 16.9 / 2.5, 2.7 | -27.3 | **1299 x 1050** |

* ~30 dBi realized needs the **300 mm side** (30.2 dBi; 250 mm gives 28.4). The 4-hex board (2 x 2 honeycomb, extent 2.5 sqrt(3) s x 3.5 s, margins
  excluded; the 150 mm board of Phase 25 is 650 x 525 mm) is then about **1.30 x 1.05 m**.
* Gain grows ~ +3.1 dB per step of 50 mm at the start and +1.8 dB at the end: the aperture efficiency vs 4 pi A cos(30)/lambda^2 stays 0.26-0.33 because
  the launch fraction (0.62) and the table-limited conversion (alpha <= 3-4 Np/m) are fixed; a larger aperture converts more of the surface wave
  (0.54 -> 0.74) but the pin's direct radiation (0.38 of the accepted power) does not shrink, so grazing+back stays at 0.45-0.53 and LHCP at 0.26-0.31.
* Model caveats carried to the larger apertures: idealised phase (no cell-to-cell phase error over 300 mm = 8.5 lambda; the nearest-K table error
  and the oblique-modulation approximation accumulate with path length), constant -0.68 dB calibration offset, pin direct and edge power not yet separated.

## L. 4-hexagon tensor-CP board: scope (`fourhex_scope.json`, `scripts/research/phase26_fourhex_scope.py`) -- NOT exported

Layout of `four_port_hex150_gap0` (feeds (64.952, 112.5), (-194.856, 112.5), (-64.952, -112.5), (194.856, -112.5) mm; hexagon side 150 mm, gap 0),
beams theta 30 at phi 45 / 135 / 225 / 315, all RHCP, pin 6.0 mm above the surface, same cell family (nearest-grid table). Model only.

**Lattices (this decides the layout).** With every port's lattice centred on its own feed (as in `four_port_hex150_gap0`, where cells
were 5.85 mm squares with a 0.55 mm gap) the D 6.0 mm circles of adjacent hexagons overlap: **64 cells of one hexagon lie closer than
6.0 mm to a cell of its neighbour (minimum centre distance 2.2 mm)** because the feed spacing (129.9, 225) mm is not a multiple of 6.4 mm. A
shared 6.4 mm lattice with cells assigned by hexagon (centre inside) has no overlap (minimum distance 6.4 mm) and 5684 cells (4 x 1408-1434).
=> **use ONE lattice for the board**; the feed pins then fall off cell centres: remove cells within 3.2 mm of a feed (done).

**Per-port prediction, driven hexagon alone** (anchored to the 150 mm HFSS run; same table, march and pattern as the single hexagon):

| Port | beam phi | cells | conversion | gain [dBi] | AR peak / worst in -3 dB beam [dB] | XPD [dB] | SLL [dB] |
|---|---|---|---|---|---|---|---|
| 1 | 45 | 1434 | 0.54 | 23.0 | 2.6 / 3.2 | 16.7 | -22.8 |
| 2 | 135 | 1408 | 0.54 | 22.95 | 2.4 / 3.0 | 17.1 | -22.5 |
| 3 | 225 | 1434 | 0.54 | 23.0 | 2.6 / 3.2 | 16.7 | -22.8 |
| 4 | 315 | 1408 | 0.54 | 22.95 | 2.4 / 3.0 | 17.1 | -22.5 |

(Model XPD/AR are pessimistic vs the single-hex HFSS: 23 dB and 1.3 dB. The hexagon is symmetric under 180 deg rotation: ports 1 = 3, 2 = 4;
ports 1/2 differ only by the angle of the beam to the hexagon's vertices.)

**Unconverted surface wave crossing into the neighbours.** Of the accepted power, 0.62 launched x (1 - 0.54) = 0.29 reaches the hexagon edge. The
neighbours' cells are rotating-slot cells of THEIR design: the wave arrives travelling toward their feed (counter to their design), so
it is NOT absorbed by loss alone: it excites the opposite circular harmonic with the cells' own alpha (ASSUMPTION: the same alpha as designed)
and radiates as an LHCP-dominated secondary lobe, not as part of the beam. Model, driven port alone, all three neighbours' cells excited
(`crossing` in the JSON):

| Driven port | pickup power / main aperture power | share of accepted power | strongest pickup lobe (theta, phi), level vs main peak | main beam with pickup: gain change, AR worst, XPD |
|---|---|---|---|---|
| 1 and 3 | 0.126 | 0.042 | LHCP (29, 76) / (29, 256): -20.0 dB (RHCP lobes <= -29.7 dB) | -0.5 dB (directivity), 4.3 dB (from 3.2), 16.5 dB |
| 2 and 4 | 0.093 | 0.031 | LHCP (21, 206) / (21, 26): -21.6 dB (RHCP <= -29.4 dB) | -0.4 dB, 3.4 dB (from 3.0), 16.0 dB |

The pickup is 3-4 % of the accepted power, 77-79 % of it LHCP, in lobes 20-22 dB below the main peak and 30-50 deg from the driven beam; it costs about 0.5 dB
of aperture directivity (the main beam's 8 deg share 0.95 -> 0.85-0.87), leaves the SLL (-22 dB) and XPD alone and worsens the worst AR over the beam by 0.4-1.1 dB.
**A boundary row of cells or a gap strip is not indicated** by these numbers: a bare strip carries the wave unchanged (alpha 0, so more arrives),
and no cell in this family absorbs it. Stopping the crossing needs lower edge power (higher alpha near the boundary, table-limited at ~3-4 Np/m) or an
absorber; neither is justified for a -20 dB LHCP lobe. If all four ports are fed at once, port k's LHCP pickup lands near the other
beams and adds to their cross-polar floor at the -20 dB level.

**Port-to-port coupling (estimate).** Capture width of a launcher with TM0 share 0.62 = 2.65 mm; the driven wave's power per radian decays
through its own hexagon (unconverted 0.40-0.60 on the rays to the neighbours' feeds) and through the neighbour's cells (alpha of the cells) over
the 259.8 mm to the neighbour's feed (cylindrical spreading included): **raw estimate -36.3 ... -37.9 dB for all five adjacent pairs**
(1-2, 1-3, 1-4, 2-3, 3-4; 2-4 is 450 mm). The same estimate applied to the HFSS 1-bit board gives -32.6 dB where the HFSS `.s4p`
has -51.5 ... -73.5 dB at 10.54 GHz (<= -41.7 dB over 8-12 GHz): the raw method is ~26 dB pessimistic there, so applying that bias gives
**about -62 ... -64 dB**. Expected range for the tensor board: -37 dB (conservative, raw) to -60 dB; both are below the -10 dB-class concerns, but
the bias is one board's and unexplained, so it is an estimate, not a result. Reciprocal pair values are equal; the unconverted power
is 0.4 (diagonal pairs) and 0.6 (the pair 1-4 across the vertical edge: the apothem ray is shorter).

**HFSS size (order of magnitude).** 5684 slotted circles (each ~28 arc facets + 4 slot edges, features 0.4 mm bridges) on a 650 x 525 x 1.575 mm board with a
pad of 15 mm and a 30 mm top air box: about 4.3 x the single-hexagon volume and 4 x its cells. Scale your single-hexagon run: tetrahedra and RAM x 4-5
(a few 10^7 tetrahedra and >100 GB at the 0.2 mm feature mesh the single run needed if it used lambda/3 refinement), solve time x 5-8 for the adaptive
passes. Options to keep it tractable: one port at a time with the others terminated (needed anyway for the gain), initial mesh seeded on the patch
sheets (0.5 mm) with the substrate at 1 mm, or domain decomposition / HPC. Not measured: take the real figure from your single-hex solve.

## M. Edge finding from the single-hexagon E-line, and board A (static beams) -- EXPORTED (`fourhex/`, `hfss_exports/hex4_tensor_cp_static/`)

**E-line** (`hfss_results/hex150_tensor_cp_single/efield_line_magE_phi225_1mm.csv`, |E| at 1 mm above the surface, phi 225, 1 mm steps, 10.5 GHz,
board edge ~142 mm). Binned `rho |E|^2` (proportional to the power through the arc; bound wave, leaky field and pin space wave all mix at 1 mm, so it
oscillates by factors of ~2 bin to bin): 20-40 mm 1.1e4; 40-60 mm 0.8e4; 60-80 mm 1.7e4; 80-100 mm 0.6e4; 100-120 mm 1.2e4; 120-135 mm 0.77e4;
135-142 mm 0.22e4; 150-170 mm 0.03e4. **Reading:** substantial power reaches the rim and then drops ~10x outside the board (free-edge diffraction);
rim/inner ratio 0.7 (120-135 vs 20-40 mm) against the model's W(127)/W(30) = 0.54 along this ray: consistent within the bin scatter, NOT
discriminating (alpha cannot be fitted). The model's edge power (0.46 of the launched wave, 0.29 of the accepted power) is therefore
used with a sensitivity scale s = 1 (design march) and s = 1.7 (upper reading); neither is a fit.

**Where the unconverted power goes (design march, fraction of the launched surface-wave power; x 0.62 for the accepted power):**

| hexagon | free (outer) edge | to hex 1 | to hex 2 | to hex 3 | to hex 4 | total |
|---|---|---|---|---|---|---|
| 1 | 0.252 | - | 0.065 | 0.065 | 0.079 | 0.462 |
| 2 | 0.331 | 0.065 | - | 0.065 | 0 | 0.462 |
| 3 | 0.252 | 0.065 | 0.079 | - | 0.065 | 0.462 |
| 4 | 0.331 | 0.065 | 0 | 0.065 | - | 0.462 |

Hexagon 1 and 3 hand 0.21 of the launched wave to neighbours (0.13 of accepted); 2 and 4 hand 0.13 (0.08 of accepted). The rest leaves the board edge.

**Board A prediction** (shared 6.4 mm lattice, 5685 cells, port driven alone; HFSS-anchored, s = 1 / 1.7; `fourhex_design.json`; 1 deg grid):

| port | beam | model peak (theta, phi) | gain dBi | XPD dB | AR peak / worst in beam | pickup / main power | neighbour LHCP lobe (vs main RHCP peak) |
|---|---|---|---|---|---|---|---|
| 1 | (30, 45) | (29.7, 45.0) | 23.1 / 23.1 | 16.4 / 16.3 | 2.65 / 3.9 | 0.13 / 0.22 | -19.9 / -17.6 dB at (30, 75), 15 deg off |
| 2 | (30, 135) | (29.7, 135.2) | 23.0 / 23.0 | 16.4 / 16.0 | 2.65 / 3.4 | 0.09 / 0.16 | -21.7 / -19.4 dB at (21, 206), 30 deg off |
| 3 | (30, 225) | (29.7, 225.0) | 23.1 / 23.1 | 16.4 / 16.3 | 2.65 / 3.9 | 0.13 / 0.22 | -19.9 / -17.6 dB at (30, 255), 15 deg off |
| 4 | (30, 315) | (29.7, 315.2) | 23.0 / 23.0 | 16.4 / 16.0 | 2.65 / 3.4 | 0.09 / 0.16 | -21.7 / -19.4 dB at (21, 26), 30 deg off |

**Correction to section L:** the strongest neighbour LHCP lobes are 15 deg (ports 1, 3) and 30 deg (ports 2, 4) from the driven beam, not "30-50 deg"
(that range is the pickup's RHCP lobes at 39-40 deg, -27...-30 dB). With s = 1.7 the lobes rise to -17.6 / -19.4 dB. XPD/AR values are the pessimistic
model values (HFSS single hexagon: 23 dB / 1.3 dB). Coupling and size notes as in L. Export README and `predictions.json` carry the check items
(gain +-1 dB, neighbour LHCP level and direction, adjacent-pair coupling, AR, S11).

## N. Board B (Butler-ready): same lattice, cells and pin, all four hexagons at (30, 225) RHCP, phase-coherent -- SCOPE ONLY, NOT exported

**Design.** The same march per hexagon (eta_design 0.9, taper 0.5, hand +1). The only change from the single-hexagon design is the phase
reference of the polarisation rotation `chi = phi_ref - kt u.(r - r_ref)`: `r_ref` is the board centre instead of the hexagon's feed. Equal-phase feed
then gives ONE continuous aperture phase `exp(-j kt u.r)` over the board. Consequence: every slot of hexagon j is rotated by the constant
`-1/2 kt u.r_j` (mod 180 deg) relative to the single-hexagon design: port 1 +37.0 deg, port 2 -4.3, port 3 -37.0, port 4 +4.3 (checked on port 3
against board A: -37.03 deg uniform, `board_B_minus_A_slot_angle_offset_deg`). Cell sizes, alphas, positions and the pin are unchanged.
Export files (when asked): `results/phase26/fourhex/fourhex_B_cells.json` (5685 cells).

**Port order for the Butler outputs.** The four feed x positions are -194.9 (port 2), -64.95 (port 3), +64.95 (port 1), +194.9 (port 4) mm: equally spaced by
129.9 mm. A Butler output n has phase n*step, so the outputs must be ordered along the x axis: **Butler outputs 1, 2, 3, 4 = ports 2, 3, 1, 4.**
That is the only ordering (up to reversal) for which the phase is linear in position (the projection of the rhombic 2x2 arrangement onto x is
the only equally spaced one), so the beam moves in `u = sin(theta) cos(phi)` only (v fixed). The staggered rows (y = +/-112.5 mm) need no
correction: the board-centre phase reference already contains them. Reversal of the order = reversal of all steps.

**Per port alone (a port driven alone, others terminated; s = 1 / 1.7):**

| port | peak (theta, phi) | gain dBi | XPD | AR peak / worst |
|---|---|---|---|---|
| 1 | (28.7, 224.2) | 22.9 / 22.9 | 17.0 / 16.4 | 2.5 / 3.4 |
| 2 | (29.9, 223.2) | 22.8 / 22.8 | 17.0 / 17.1 | 2.5 / 3.9 |
| 3 | (30.1, 226.2) | 22.6 / 22.5 | 17.3 / 16.9 | 2.4 / 3.1 |
| 4 | (29.4, 226.7) | 22.9 / 22.9 | 17.8 / 17.9 | 2.3 / 3.3 |

(each hexagon alone has a beam 12 deg wide and sits within 1.6 deg of the design direction; pickup/main power 0.11-0.15, s = 1.7: 0.19-0.26.)

**Ideal 4 x 4 Butler, equal amplitude (0.5 each, 1 W accepted in total), model, s = 1 (s = 1.7 in brackets where different):**

| input (step) | port phases (deg) p1 / p2 / p3 / p4 | array theory (theta, phi) | model peak (theta, phi) | RHCP gain dBi | AR peak / worst | XPD | strongest other RHCP lobes (dB re peak; offset in u, v) |
|---|---|---|---|---|---|---|---|
| equal phase | 0 / 0 / 0 / 0 | (30.0, 225.0) | (29.9, 225.0) | 28.6 (28.4) | 2.5 / 3.0 | 16.9 | -11 (-0.01, +0.09), -14 (+0.09, -0.05), -14 (-0.09, +0.06) |
| 1L (-45) | -90 / 0 / -45 / -135 | (28.75, 227.3) | (28.9, 227.2) | 28.6 | 2.3 / 2.7 | 17.5 | -8 (-0.08, -0.04), -9 (-0.09, +0.06), -11 |
| 1R (+45) | 90 / 0 / 45 / 135 | (31.3, 222.9) | (30.9, 223.0) | 28.3 | 2.6 / 3.1 | 16.7 | -6 (+0.09, +0.04), -8 (+0.09, -0.06), -11 |
| 2R (+135) | -90 / 0 / 135 / 45 | (34.1, 219.1) | **(27.1, 222.6)** | 26.9 (27.0) | 2.1 / 3.0 | 18.3 | **-0.9 at (31, 231)**, -3 at (33, 220) |
| 2L (-135) | 90 / 0 / -135 / -45 | (26.5, 232.5) | **(29.2, 218.7)** | 26.2 | 2.3 / 3.4 | 17.7 | **-0.1 at (33, 227)**, -1.4 at (27, 231) |

**Beam direction.** The 1D steering law is `du = -step / (k d)` with the outputs 129.9 mm apart: ideal beam shifts are **+/-1.7 deg for +/-45 deg and
+/-5.1 deg for +/-135 deg** (in angle; u shifts +/-0.027 and +/-0.082) -- not +/-0.9 / +/-2.7 deg. Your numbers are exactly what the 259.8 mm spacing of
the same-row pair (ports 2-1 or 3-4) gives (lambda/259.8 = 0.110 in u -> 6.3 deg grating lobes): **the 7 deg grating-lobe estimate is confirmed,
but it belongs to the tile period, not to the 129.9 mm output spacing.** The model confirms the +/-45 deg inputs to 0.4 deg (28.9 vs 28.75; 30.9 vs 31.3).

**Grating lobes / tile lattice.** The four hexagons are tiles of period 259.8 mm (a1 = (259.8, 0), a2 = (129.9, 225)): lobes of the tile-amplitude
modulation sit at du = lambda/259.8 = 0.110 (6.3 deg), seen at -11...-14 dB with equal phase and at -6...-9 dB for the +/-45 deg inputs (the steered beam moves
toward one of them, where the tile element pattern is higher). For +/-135 deg the intended beam (shifted 5 deg) and the neighbouring tile lobe have comparable level:
**the model gives no clean single beam at the nominal 5 deg shift**: the strongest lobe is at (27, 223)/(29, 219) (gain 26.9/26.2 dBi), with a second lobe -0.9/-0.1 dB
at 8-9 deg away. Honest summary: 1L and 1R work as a +/-1.7 deg scan at 28.3-28.6 dBi with first lobe 6-8 dB down; 2R and 2L do not produce a single beam
on this board (peak 26-27 dBi, split in two lobes). The 4 x 4 Butler's +/-135 deg beams therefore need either a tile-period amplitude equaliser (a taper that
continues across the shared edges) or acceptance of the split beam. Mutual coupling at -37 ... -60 dB moves the port phases by < 1 deg: negligible.

**Ports-to-Butler Edit Sources table** (HFSS e^{j omega t}, which is the convention of the model: the wave is exp(-j beta rho); equal magnitude, e.g. 1 W each --
HFSS normalises gain to the total incident power, check the "Incident power" row of the antenna-parameter table):

| Butler input | Port 1 | Port 2 | Port 3 | Port 4 | expected beam (theta, phi) |
|---|---|---|---|---|---|
| 1L (-45) | -90 deg | 0 deg | -45 deg | -135 deg | (28.9, 227.2) |
| 2R (+135) | -90 deg | 0 deg | 135 deg | 45 deg | (27.1, 222.6) + (31, 231) lobe |
| 2L (-135) | 90 deg | 0 deg | -135 deg | -45 deg | (29.2, 218.7) + (33, 227) lobe |
| 1R (+45) | 90 deg | 0 deg | 45 deg | 135 deg | (30.9, 223.0) |
| equal phase (reference) | 0 | 0 | 0 | 0 | (29.9, 225.0), 28.6 dBi |

First check in HFSS: 1L must move the beam to smaller theta/larger phi (28.9, 227); if it moves the other way, negate all phases (convention check).
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, SOLVER_SUPERCELL; the pickup uses the cells' designed alpha for the crossing wave (ASSUMPTION).

## O. Board B: tile-amplitude equalisation (`fourhex/tile_equalise.json`, `scripts/research/phase26_tile_equalise.py`) -- TARGETS NOT REACHED

Question: can flattening the aperture amplitude across the shared edges bring the tile/grating lobes to <= -15 dB (equal phase) and <= -10 dB
(+/-45 deg Butler inputs)? Metric: strongest RHCP local maximum 4-20 deg from the peak, dB re peak, main patterns without the pickup
(baseline -10.3 here vs -11 in section N, which also had the pickup and a 9 deg cone). Additive solver hook: `march_rays(..., density_fn=)`
(default `None` is bit-identical; tested).

| design (target density change only; cells stay table-limited) | radiated eta | equal-phase gain / lobe | 1L gain / lobe | 1R gain / lobe | 2R / 2L lobe |
|---|---|---|---|---|---|
| baseline (eta_design 0.9, taper 0.5) | 0.54 | 29.1 / -10.3 | 28.7 / -7.4 | 28.4 / -5.6 | -0.9 / -0.1 |
| uniform target density | 0.50 | 28.9 / -10.8 | 28.6 / -7.8 | 28.2 / -6.0 | -0.9 / -0.1 |
| outward boost b = 4 (flatten the depletion) | 0.39 | 28.0 / -12.6 | 27.5 / -8.5 | 27.2 / -6.8 | -1.0 / -0.1 |
| taper only toward free (outer) edges, 60 mm | 0.52 | 28.9 / -10.4 | 28.6 / -7.5 | 28.3 / -5.7 | -0.9 / -0.1 |
| low target density eta_design 0.4 / 0.3 | 0.37 / 0.29 | 27.8 / -12.4; 26.7 / -12.5 | 27.4 / -8.4; 26.3 / -8.5 | 27.1 / -6.8; 26.0 / -6.8 | -1.1 / -0.0 |
| board-wide cos^2 taper (330 x 270 mm), eta_design 0.9 / 0.3 | 0.48 / 0.30 | 28.1 / -11.0; 26.2 / -11.3 | 27.8 / -7.8; 25.8 / -8.0 | 27.4 / -5.9; 25.5 / -6.2 | -1.0 / -0.2 |
| quiet zone: no radiation within 50 mm of each pin (ideal: cells kept, currents removed) | 0.38 | 27.4 / **-13.2** | 26.9 / -8.6 | 26.6 / -7.0 | -1.1 / -0.0 |

**Result: no variant reaches -15 dB (equal phase) or -10 dB (+/-45 deg); the best is -13.2 / -8.6 / -7.0 dB, at 1.5-2 dB gain cost.**

**Why (two floors that amplitude shaping cannot move):**
1. *Outline floor.* An IDEAL uniform-amplitude, ideal-phase circular aperture on the same 5685 cells has its first lobes at **-12.3 dB** (the board outline
   is 650 x 525 mm, ~23 lambda: uniform-aperture sidelobes). Only an ideal smooth board-wide taper would reach -18.9 dB (reference in the JSON), and the
   feed/cell family cannot produce it: the table's alpha window is 0.42-4.3 Np/m, so the radiated density near each pin cannot go below ~90 (units of the design amplitude^2)
   while it saturates at 26-40 over the outer 60 % of each tile (alpha 3.4 Np/m x depleted W / rho). The envelope is peaked at the four pins; lowering the
   target only trades gain for a -12.5 dB lobe.
2. *Staircase floor for the Butler inputs.* With uniform amplitude and the SAME phases the steered inputs still give 1L -8.4 / 1R -6.7 dB and the +/-135 deg inputs
   no clean beam (-1.1 / -0.0 dB): one feed phase per 260 mm tile is a 4-step staircase of a linear phase ramp (45 deg or 135 deg per 130 mm of tile centre),
   which is a quantisation lobe, not an amplitude effect. No amplitude equalisation changes this.

**The equal-phase case on its own, "combined-beam board" (1:4 in-phase divider, 0.25 W per port; model with pickup, section N):** one beam at (29.9, 225.0),
**28.6 dBi RHCP (28.4 with s = 1.7), AR 2.5 dB at peak / 3.0 worst in beam, XPD 16.9 dB** (model; HFSS single hexagon was 23 dB / 1.3 dB), +5.5 dB over one
port (ideal +6.0), first lobes -11 dB at 5-6 deg (outline sidelobes of the 23-lambda board; -13.2 with the quiet-zone variant at 1.2 dB less). The RHCP lobes beyond 20 deg are
<= -26 dB. This board needs only the in-phase 1:4 divider (no Butler) and is the cleaner use of board B.

## P. Decision: board B = combined-beam board; Butler scanning on this layout is shelved

* **Board B** is the "combined-beam board": all four ports equal amplitude and equal phase (an ideal 1:4 in-phase divider, 0.25 W per port), one beam at
  (30, 225) RHCP. Design = the baseline of N/O (eta_design 0.9, taper 0.5, global phase reference; slot rotations vs the single hexagon +37.0 / -4.3 / -37.0 / +4.3 deg
  for ports 1..4): **28.6 dBi (model), AR 2.5 / 3.0 dB, XPD 16.9 dB (model), first lobes -11 dB.** The quiet-zone variant (-13.2 dB lobes, -1.2 dB gain) is recorded
  as an OPTION only, not exported.
* **Butler scanning on the 4-hex layout is shelved.** Reason: quantisation lobes from one feed phase per 260 mm tile (a 4-step staircase of the linear ramp;
  section O: 1L -8.4 / 1R -6.7 dB even with ideal uniform amplitude, +/-135 deg inputs no single beam). A clean Butler scan needs feed spacing of about lambda/2
  (a different architecture, not this board). Sections N's Butler tables are kept for reference only; no Butler Edit Sources table is to be used.
* **Export of board B is held** until the board A HFSS results are in (board A being run by the user).
