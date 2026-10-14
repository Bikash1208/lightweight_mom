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
