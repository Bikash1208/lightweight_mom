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
