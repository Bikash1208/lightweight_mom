# Phase 26 — leakage-controlled, tensor-impedance, pure-CP hexagon (one port first)

**Status:** active brief, written 2026-10-14. Read `docs/START_HERE.md`, `docs/STATUS.md`, `CLAUDE.md` first. Branch
`phase26-leaky-tensor-cp` from **`phase25-butler-sector`** (it holds the Phase 25 forward model, the locality checks and the
STATUS caveat). Phase 25 B/C (Butler sector scanner) is **paused**: move its brief to `docs/briefs/_parked/` with a
one-line note pointing here. All repo rules apply (one commit per item, fast suite + snapshot gates green, physics
tests, honest labels, additive-only solver changes with bit-identical old results).

## 1. What HFSS showed (first real validation of the app's boards)
User's HFSS of `hfss_exports/four_port_hex150_gap0` (files in `hfss_results/hex150_gap0/`: `.s4p` stays local/gitignored;
commit the antenna-parameter table and the realized-gain CSV and index them in `docs/validation/README.md`):
- **S-parameters (8–12 GHz):** S11…S44 −13.9/−17.6/−16.4/−14.6 dB at 10.54 GHz, −10 dB band ≈ 8.5–11.4 GHz (minimum
  near 9.6–10 GHz → monopole slightly long); coupling ≤ −50 dB at 10.54 GHz, ≤ −41.7 dB anywhere in 8–12 GHz.
  **No gap is needed for isolation.**
- **Port 3 alone, 10.5 GHz:** beam at **(θ 30°, φ 224°)** (design 30°, 225°); peak realized gain **19.85 dBi**
  (directivity 20.0 dBi) vs model **25.8–26.0 dBi**; total efficiency 97.7 %; aperture efficiency ≈ 11 % (model ≈ 40 %);
  XPD at the peak 22 dB; LHCP shoulders −5.8 / −7.7 dB beside the beam; co-pol lobes beyond the first −16 to −18 dB.
- **Power budget (integrated from the CSV, 2° grid, ±5 %):** main beam (≤ 8°) ≈ 20 %; within 20° ≈ 36 %; θ > 70° ≈ 40 %;
  back hemisphere ≈ 15 %; ≈ half the radiated power is LHCP.
- **Diagnosis:** the hologram pattern is right (direction, beam XPD, far lobes), but the surface wave is **not converted
  within the hexagon** — it reaches the edges and diffracts (grazing + back radiation). Consistent with Phase 25's finding
  that the 0.55 mm-gap 1-bit pair modulates far more weakly than the cell-wise model assumes. The half-power LHCP is the
  scalar-cell CP ceiling (STATUS) plus the unpolarised edge radiation.

## 2. User requirement (hard)
**Pure circular polarisation (RHCP) — no LP/LHCP mixing.** Measurable targets (model, then HFSS):
- axial ratio ≤ 3 dB at the beam peak and over the −3 dB beamwidth;
- ≥ 90 % of the power in the main-beam region is RHCP; integrated LHCP over the full sphere ≤ 10 % of radiated power;
- peak XPD ≥ 25 dB.
Supporting targets: ≥ 85–90 % of the launched surface-wave power radiated inside the hexagon (model); HFSS power at
θ > 70° plus the back hemisphere ≤ 15 %; single-port realized gain ≥ 24 dBi (aperture efficiency ≥ 35 %); co-pol SLL ≤ −15 dB.
Same board: one 150 mm-side hexagon, monopole at the centre, RT/duroid 5880 1.575 mm, 10.54 GHz, beam θ 30°, φ 225°.
Scalar cells cannot meet the CP target → this phase needs **tensor (anisotropic) cells**.

## 3. Method (literature)
- Tensor holographic impedance surfaces and CP from a TM surface wave: Fong et al., IEEE TAP 58(10) 2010; Minatti et al.,
  "Synthesis of modulated-metasurface antennas with amplitude, phase, and polarization control," IEEE TAP 64(9) 2016
  (review [10]); Teniou et al., IEEE TAP 65(11) 2017 and 67(7) 2019 (review [23], [24]); Faenzi et al., Sci. Rep. 2019
  (review [12]). Adiabatic Floquet-wave analysis for local leakage: Minatti et al., IEEE TAP 64(9) 2016 (review [11]).
- Penetrable-sheet vs impenetrable-impedance definition matters (Martini et al., EPJ AM 2020, review [27]) — state which.

## Work items (stop points marked ■)

**A. Calibrate the model against today's HFSS (must come first).**
- Add a **leakage model** to the forward prediction: local leakage constant α(ρ) from the local modulation (adiabatic
  Floquet / Oliner–Hessel closed forms for the radiating n = −1 harmonic), surface-wave power decaying along each ray,
  the unconverted remainder reaching the edge reported (and, optionally, an edge-diffraction estimate).
- Re-predict the **existing** hex150 port-3 board with it (1-bit 0.55 mm-gap pair, its real per-cell contrast taken from
  the Phase 25 supercell MoM rather than the isolated-cell Γ). It should move toward HFSS: ~20 dBi, ~40 % grazing.
  Report predicted vs HFSS side by side (named states). **■ Stop if the calibrated model still misses HFSS by > 2 dB in
  directivity — report why before designing anything new.**

**B. Tensor cell family and extraction.**
- Cell: rotated rectangular patch (the `Geometry.rotation` field already exists) — sweep L, W, rotation on the 6.4 mm cell
  (or justify another period), staying in the **non-resonant** range where Phase 25's locality checks hold.
- Extraction: from the MoM's full 2×2 normal-incidence reflection (xx, yy, xy), derive the 2×2 **tensor sheet admittance**,
  then the surface-wave dispersion for any in-plane direction with a tensor transverse-resonance (sheet) model — the
  tensor generalisation of `dispersion/sheet_model.py` (additive; scalar case bit-identical). Do **not** rely on
  `dispersion/mom_scaffold.py` (its MoM eigen-search only finds the bare-slab pole, documented there).
- Tests: rotation equivariance (cell rotated by α ≡ tensor rotated by α); square patch → isotropic tensor; reciprocity
  (symmetric tensor for a lossless reciprocal cell); scalar limit reproduces the existing `sheet_model` numbers exactly.
- Database: (L, W, rotation) → tensor; report coverage of the tensor space needed by item C. **■ Stop with the coverage
  map and the extraction validation.**

**C. CP tensor synthesis with leakage/amplitude control (one hexagon).**
- Target aperture field: pure RHCP pencil at (30°, 225°) with an amplitude taper (choose for SLL ≤ −15 dB) and leakage
  profile α(ρ) so ≥ 85–90 % of the surface wave is radiated within the hexagon.
- Tensor modulation per Minatti 2016 / Fong 2010 (state the formulas used and the impedance definition); map to
  (L, W, rotation) per cell via the B database; continuous geometry (no 1-bit).
- Predict with the calibrated forward model (A): gain, AR over the beam, RHCP fraction, integrated LHCP, XPD, SLL,
  grazing/back power, at 10.33/10.54/10.75 GHz. Tests: the synthesis recovers a known tensor surface in a toy case; the
  predicted polarisation of a uniform-rotation test surface matches analytic expectation.
- **■ Stop with the predicted table vs the §2 targets.**

**D. Export the single hexagon** (exact STEP + AEDT script + layout, rotated rectangles drawn exactly) to
`hfss_exports/hex150_tensor_cp_single/` with a README: monopole details (note the slightly-long pin from §1 — suggest a
trimmed pin length or a parametric pin), HFSS setup (one port, radiation boundary/PML, far-field sphere θ 0–180° 1°,
φ 0–360° 1°, antenna parameters + RHCP/LHCP/axial-ratio exports), and the predicted numbers to compare.

**E. Close-out:** STATUS (capability rows: leakage model, tensor cells, CP synthesis — model-only until HFSS; the
§1 HFSS validation row for hex150_gap0), CHANGELOG, docs/validation index, move this brief to history.

## Do not
- Use 1-bit/coded cells for this design; claim CP purity or gain without the calibrated model (A) and, finally, HFSS.
- Change existing results (snapshot gates), rewrite history, or start the 4-hexagon board before the single hexagon is
  validated in HFSS.
