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
