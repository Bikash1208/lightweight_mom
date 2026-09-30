# Phase 15 brief — per-beam polarization, pencil vs vortex beams, vector pattern model, multi-port, coded method

Read first: `CLAUDE.md`, `docs/theory/CONVENTIONS.md` (exp(+jωt)!), `src/hologram/hmsa_synthesis.py`,
`src/pattern/aperture_pattern.py`, `src/app/routers/hmsa.py`, `src/app/static/index.html`, and the reference paper
Ghosh, Ghosh, Saha, IEEE OJAP 2026, DOI 10.1109/OJAP.2026.3694994 (user has it; copy to `docs/references/` if given):
eqs. (1)–(17) and Table 1 are the methodology source for this phase. Commit per part; fast suite green each time.
Existing single-polarization designs must reproduce bit-identically when the new options are left at their defaults.

## Why (measured by the user in HFSS)
Every app design fed by the monopole radiates a **vortex (doughnut) beam** in HFSS although the app's pattern shows a
pencil beam. Cause: the monopole's TM₀ surface wave has a **radial** tangential field, so the surface current is
radial, J = ρ̂ · J_s(x,y), with ρ̂ = x̂ cos φ′ + ŷ sin φ′ = ½[e^{+jφ′}(x̂ − jŷ) + e^{−jφ′}(x̂ + jŷ)] (φ′ = azimuth about
the feed). Each circular component carries e^{±jφ′} → an l = ±1 OAM vortex at the beam direction. The app's pattern
model is **scalar**, so it cannot see this. Fix = polarization-compensation term φ_ext in the hologram (paper Table 1;
Minatti et al. TAP 2016), plus a **vector** pattern model so the app predicts what HFSS shows.

## Part A — per-beam polarization and beam type in the synthesis
For each object beam b: direction (θ_b, φ_b), relative power, **polarization p̂_b** ∈ {LHCP, RHCP, HP (x̂), VP (ŷ),
slant ±45°, "Natural (radial — no compensation)"} and **beam type** ∈ {Pencil, Vortex with OAM charge l_b (±1, ±2, …)}.
- Compensation term, derived (not copied) from projecting the wanted polarization onto the local current direction:
  the field radiated into p̂_b is ∝ (p̂_b* · ρ̂) J_s, so choose
  **φ_ext,b(x,y) = exp(−j·arg(p̂_b* · ρ̂(x,y))) · exp(+j·l_b·φ′(x,y))**, with l_b = 0 for Pencil.
  - CP: p̂*·ρ̂ = e^{±jφ′}/√2 → φ_ext is a spiral e^{∓jφ′} (paper Table 1 writes e^{∓j(tan⁻¹((y−y₀)/(x−x₀)) + π/4)};
    the π/4 is a constant — keep it configurable).
  - HP: p̂·ρ̂ = cos φ′ → φ_ext = 0/π flip across the line through the feed perpendicular to x (paper: e^{±jnπ} by sign
    of y_j — **check which half-plane convention is correct in our coordinates by the vector test in Part B**).
  - VP: sign of sin φ′. Slant: sign of cos(φ′ − 45°).
  - "Natural": φ_ext = 1 (today's behaviour → vortex). Vortex with charge l: add e^{jlφ′} on top of the compensation.
  Derive the handedness signs from CONVENTIONS.md (exp(+jωt), IEEE RHCP definition) and pin them with tests; the
  paper's sign convention may differ.
- Hologram (both methods): Z(x,y) = jX̄ [1 + M · Re{ Σ_b w_b · φ_ext,b · φ_ref* · φ_obj,b }] (full-aperture
  multiplexing, paper eq. 16) with power weights w_b normalised; **option: sector multiplexing** (paper eq. 17 —
  each beam owns a sector/quadrant of the aperture; equal sectors by default, user-assignable) — report the gain
  difference (paper: sector ≈ +3 dB per LP beam).
- Feed position (x₀, y₀) is taken from the reference field's NaN/singularity location (default) or entered.

## Part B — vector aperture/pattern model (replace the scalar one for HMSA)
- Aperture current: J(x,y) = ρ̂(x,y) · J_s(x,y), J_s = E_ref · (X_achieved/X̄ − 1) (method A) — keep the documented
  first-order model; allow using the magnitude file for amplitude, and optionally an exponential leakage decay
  e^{−α ρ} (α from the paper's eq. 11–12 form, default 0 = off, documented).
- Far field of a horizontal current sheet on the **grounded substrate**: use this project's own layered Green's
  function / TL model to get the exact TE/TM spectral response at (k_x, k_y) = k₀ sin θ (cos φ, sin φ) (reciprocity:
  plane wave from (θ, φ) → tangential field at z = 0), rather than a free-space dipole with an image. Output E_θ, E_φ;
  then co/cross by Ludwig-3 for LP and E_R/E_L for CP (exp(+jωt) convention), axial ratio, XPD per beam.
- Show per beam: achieved direction, polarization, XPD / axial ratio, and **on-axis null depth + phase winding
  number** around the beam (so "vortex" vs "pencil" is measured, not assumed). Plots: co-pol and cross-pol u–v maps,
  cuts through each beam, and for vortex beams a phase plot around the beam axis.
- **Validation tests (must pass):**
  1. "Natural" polarization, single beam θ=20°, φ=0 on the user's walk-through board → vector model shows a null at
     the beam direction with |winding| = 1 in the circular components (the vortex HFSS showed).
  2. Same with p̂ = RHCP (and LHCP, HP, VP) → pencil beam, peak at the beam direction, XPD ≥ 15 dB (report actual).
  3. Vortex l = +2 requested → null with winding +2 in the requested polarization.
  4. Scalar-model regression: with the vector model disabled, previous results unchanged.
  5. Reproduce the paper's Design-1 **qualitatively** (25 GHz, RT/duroid 5880, L_g = 3.6 mm, 88×88 cells, X̄ =
     165.43 Ω, M = 0.45, beams HP (−35°,0), LHCP (−18°,0), RHCP (25°,0), VP (45°,0)) with a synthetic TM₀ reference
     (β from the analytic slab model at 25 GHz): four beams at the right angles with the right polarizations. Don't
     claim gain agreement — the paper uses its own unit cell (gap-parameter square patch, HFSS eigenmode Z_s).

## Part C — multi-port (multiple monopole feeds)
Primary reference: Mishra, Barad, Sarkar, "A Multi-Port Multi-Beam Binary Coded mm-Wave Metasurface Antenna for
Full-Duplex Vehicular Communications," IEEE TVT 2026, DOI 10.1109/TVT.2026.3733618 (`docs/references/`): five coaxial
monopoles on one grounded RT/duroid 5880 board (26 GHz), **each port feeding its own sub-aperture** (ports 1–4:
5λ×5λ single beams; port 5: 10λ×10λ quad beam), inter-port isolation > 32 dB measured, AR < 3 dB, 13.8–16.7 dBic.
- Model a board as a set of **sub-apertures** (rectangles, user-drawn on a grid or from presets: 2×2 + 1 large,
  1×N strip, custom), each with: its own feed position (monopole at the sub-aperture centre by default, or imported
  `.fld` per port), its own beam list (direction, polarization, pencil/vortex), and its own method. The reference wave
  of each sub-aperture is radial about **its own** feed (φ′, ρ̂, φ_ext all per port).
- Also support the single-aperture **superposition** alternative (Ghosh OJAP eq. 5: several feeds illuminate one
  shared hologram) as an option, with the documented efficiency penalty.
- Outputs per port: pattern (vector model), gain estimate, polarization/AR, and an **isolation estimate** between
  ports (first-order: surface-wave coupling across sub-aperture boundaries — at least flag adjacent-edge distance in
  λ_sw and recommend a gap/guard row; say plainly that true isolation needs HFSS). Aperture efficiency per port.
- Export: one STEP with all sub-apertures + one marker per feed; the HFSS README lists every port position and
  suggests one coax/lumped port per monopole and the S-parameter matrix to check isolation.
- Validation: reproduce Mishra's layout qualitatively — 5 ports, beams (15°,0°), (30°,90°), (45°,180°), (60°,270°)
  from ports 1–4 and all four from port 5 — correct directions and LHCP from the vector model.

## Part E — coded (1-bit / N-bit) design method (Mishra TVT 2026 methodology)
Add a third design method next to "Surface-wave (eigenmode)" and "Reflection phase": **Coded (phase-encoding)**.
- Continuous aperture phase per cell (Mishra eq. 1):
  ψ_pq = −k₀(x_p sin θ₀ cos φ₀ + y_q sin θ₀ sin φ₀) + β_sw·√(x_p² + y_q²) + χ·tan⁻¹(y_q/x_p),
  χ = ±1 for LH/RH CP (the same spiral compensation as Part A — unify the code so both methods use one φ_ext
  implementation; χ = 0 with the half-plane flip for LP). Use β̄ (loaded) rather than k₀√ε_eff when the dispersion
  table exists, and report which was used.
- Multi-beam: f = Σ_m w_m cos ψ_m (Mishra eq. 3). 1-bit: B = 1 if f > 0 else 0 (eq. 2). N-bit: quantize arg of the
  complex sum Σ_m w_m e^{jψ_m} into 2^N levels.
- Unit cells: 1-bit needs two cells with ~180° reflection-phase difference and equal magnitude (Mishra: full patch vs
  patch with annular slot ring). The app's library currently has only square patches — add a **two-state cell
  library** option: (a) pick the best-separated pair from the existing square-patch library and report the achieved
  phase difference honestly, and (b) support a second geometry type "patch with annular slot" (parameters a, b, R, r)
  in the unit-cell solver/geometry/STEP exporter, validated like the square patch; (c) allow uploading an HFSS
  two-state table. Show the achieved difference and warn if it is far from 180°.
- Pattern: same vector model (the aperture field is the binary/quantized pattern times the radial current). Report
  quantization loss vs continuous.
- Export: STEP/DXF with the two cell types; patch-layout PNG in copper/grey shows both geometries.

## Part C′ — keep for superposition users
If a single shared aperture with several feeds is chosen, each feed's contribution uses its own ρ̂ in the vector model.

## Part D — app (Step 2 and Step 3), plain language
- Step 2 beam table gets two columns: **Polarization** (LHCP / RHCP / Horizontal / Vertical / Slant +45 / Slant −45 /
  Natural) and **Beam type** (Pencil / Vortex, with "charge" when Vortex). Defaults: Pencil, Horizontal (a linear
  default avoids surprising handedness sign errors; the user picks CP explicitly).
- Design method gains **Coded (1-bit / N-bit)**. A **Ports** section (Single / Multi-port: pick a layout preset or
  draw sub-apertures; each port gets its own beam table). Under Advanced: Multiplexing (Full aperture / Sectors),
  Leakage decay (off/on + value), Pattern model (Vector — default / Scalar legacy).
- Step 3: per-beam result rows (direction, polarization achieved, XPD or axial ratio, "Pencil ✓" or "Vortex l=…").
  Plain notes, e.g. "Vortex expected — Natural polarization keeps the monopole's radial currents."
- Patch-layout PNG, STEP (rounded mm, sheets, one Patches body), CSV/DXF work unchanged for all options.

## Deliverables
Code + tests as above; `docs/phases/PHASE15_POLARIZATION_VORTEX.md` with derivations, sign conventions, validation
numbers, and limits (first-order model; scalar isotropic cells limit XPD — tensor/anisotropic cells are the next step
per Fong 2010 / Minatti 2016); updated USER_GUIDE; fresh screenshots; HFSS notes updated with how to check
polarization (Ludwig-3 / RHCP-LHCP far-field outputs) and vortex (phase of E_R around the beam).

## Methodology sources (journal)
- Ghosh, Ghosh, Saha, IEEE OJAP 2026, DOI 10.1109/OJAP.2026.3694994 — φ_ext per polarization, full-aperture vs
  quadrant multiplexing, monopole-fed square-patch HMS (Table 1, eqs. 13–17).
- Minatti, Caminita, Martini, Sabbadini, Maci, "Synthesis of Modulated-Metasurface Antennas With Amplitude, Phase,
  and Polarization Control," IEEE TAP 2016 (ieeexplore 7508983).
- Fong, Colburn, Ottusch, Visher, Sievenpiper, "Scalar and Tensor Holographic Artificial Impedance Surfaces,"
  IEEE TAP 2010.
- "Characterization of Single-Feed Multiple Dual-Polarized Dual-Beam Scalar Impedance-Modulated Metasurface
  Antennas," IEEE TAP 2025 (ieeexplore 11008880).
- Weng et al., "Tensor-Free Holographic Metasurface Leaky-Wave Multi-Beam Antennas with Tailorable Gain and
  Polarization," Sensors 2024.
- Mishra, Barad, Sarkar, "A Multi-Port Multi-Beam Binary Coded mm-Wave Metasurface Antenna for Full-Duplex Vehicular
  Communications," IEEE TVT 2026, DOI 10.1109/TVT.2026.3733618 — multi-port sub-apertures, 1-bit coding, CP spiral term.
- S. Ghosh, H. Priyadarshnamm, C. Saha, "Characterization of Single-Feed Multiple Dual-Polarized Dual-Beam Scalar
  Impedance-Modulated Metasurface Antennas," IEEE TAP, vol. 73, no. 9, pp. 6334–6343, Sep. 2025.
- "Array Surface-Wave Launcher for the Efficient Generation of Shaped Beam and Multibeam With Metasurface," IEEE
  (ieeexplore 9465675) — multi-feed.
- OAM from holographic leaky-wave metasurfaces: Sci. Rep. 2020 (s41598-020-64278-9); Optics Letters 2019
  (ol-44-6-1452, anisotropic holographic OAM); Appl. Sci. 2021 (monopole-fed CP vortex).
