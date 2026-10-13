# Phase 25 — Butler-fed sector scanner: one shared multiplexed hologram, 4 monopoles (model study)

**Status:** active brief, written 2026-10-13. Read `docs/START_HERE.md`, `docs/STATUS.md`, `CLAUDE.md` first. Branch
`phase25-butler-sector` from `main` (`v1.1-handoff`). All repo rules apply (one commit per item, fast suite green before
each, physics tests, honest labels, no solver-core changes, model-only until the user's HFSS run).

## Idea and constraint (from the user)
Hardware is fixed to **4 monopoles + a 4×4 Butler matrix (optionally with tunable phase) + one passive metasurface**,
nothing else. Physics limit: 4 ports = 4 independent beam "slots"; every excitation is a blend of the 4 basic patterns.
Chosen use of the slots: **sector scanner** — design ONE shared hologram over the whole board so that Butler input k
gives beam k, with the 4 beams adjacent (about one beamwidth apart). Fixed Butler → 4 switchable neighbouring beams;
tunable Butler (blending adjacent inputs with a relative phase ψ) → continuous steering across the sector.
This replaces the four separate hexagon holograms (where a Butler only mixes fixed, non-overlapping beams).

Reference papers (cite in docs):
- M. Bodehou, E. Martini, S. Maci, I. Huynen, C. Craeye, "Multibeam and Beam Scanning With Modulated Metasurfaces,"
  IEEE TAP 68(3):1273–1281, 2020 — joint multibeam synthesis; steering by relative feed phase gave ~20° scan limited to
  the 3 dB beamwidth with high sidelobes; more beams → lower performance. **This is the benchmark to compare against.**
- Y. B. Li, X. Wan, B. G. Cai, Q. Cheng, T. J. Cui, Sci. Rep. 4:6921, 2014 — interferogram superposition with a
  monopole reference.
- B. H. Fong et al., IEEE TAP 58(10):3212–3221, 2010 — holographic impedance surface design.
- Crosstalk benchmark: L. Yang et al., FITEE 24(12):1776–1790, 2023 (−9 dB max between 4 channels).

## Fixed starting values (parameters, not constants)
10.54 GHz, RT/duroid 5880 1.575 mm, 6.4 mm cells, Phase 22 automatic cell states; board = the Phase 23 `hex4` cluster
outline (side 150 mm) as ONE aperture; RHCP beams. Monopole layouts to compare: (L1) the four current hexagon centres;
(L2) a compact line of 4 at λ_sw/2 spacing at the board centre. Default sector: θ = 30°, centred on φ = 45°, with the
4 beam azimuths spaced so adjacent beams cross at ≈ −3 to −4 dB for the full-board beamwidth (compute it; report the
resulting sector width). Ideal 4×4 Butler (standard hybrid/45° form; document its matrix and port↔phase-step mapping);
tunable mode = blends of adjacent inputs a(ψ) = (B_k + e^{jψ} B_{k+1})/√2, ψ ∈ [0, 360°).

## Work items

**A. Forward model for an arbitrary excitation (prerequisite — the current pattern path uses the design-time
`aperture_field`, which cannot predict a different excitation).**
- Given a FIXED cell map and an arbitrary monopole excitation vector a = (a₁…a₄): incident surface wave
  E_inc = Σ a_m E_m (reuse `hologram.multi_feed` monopole fields, complex), radiated aperture current from the
  first-order holographic model (state/modulation × E_inc, same convention the synthesis assumes), then the existing
  vector pattern. Document the approximation (first-order, no feed–feed coupling, local periodicity).
- Validate: (i) for the designed excitation it reproduces the existing synthesis prediction to numerical precision;
  (ii) on a small aperture, compare against the existing full-array MoM verify path for the designed state and one
  non-designed state; report the difference honestly (named state, not pass/fail).

**B. Multiplexed synthesis for Butler states.**
- Reference wave per Butler input: E_ref,k = Σ_m B_km E_m. Object wave per input: beam k.
- Method 1 (superposition, Li 2014 / Fong 2010): modulation ∝ Σ_k Re(E_obj,k E_ref,k*), then quantise (1-bit / 2-bit
  Phase 22 states / continuous upper bound).
- Method 2 (joint, Bodehou-style): choose the cell map that minimises, in least squares over a direction grid, the
  error between the forward-model pattern for each Butler input and its target beam (with the other beams' directions
  as penalised "crosstalk" points). Reuse the Phase 23 cell-wise linear field update and optimiser where possible;
  deterministic, seeded, no ML.
- Tests: forward model linearity (pattern(a+b) = pattern(a)+pattern(b)); single-feed degenerate case reproduces the
  existing single-port design; Butler matrix unitary; seeds reproducible.

**C. Sector-scanner evaluation (held-out 0.5° × 1° grid, at 10.33/10.54/10.75 GHz with states re-solved).**
For layouts L1 and L2 × methods 1/2 × {1-bit, 2-bit, continuous}:
- Per Butler input: beam direction vs target, directivity/gain, co-pol SLL, **crosstalk** (level at the other three
  target directions, relative to own peak), XPD, efficiency.
- Tunable sweep ψ for each adjacent pair: peak-direction trajectory, gain vs ψ (dip at crossover), SLL vs ψ;
  the continuous sector where gain stays within 3 dB of the maximum. Compare with Bodehou 2020 (~20°, high SLL).
- Direct single-monopole drive (no Butler) for completeness: show it is the predicted blend of the 4 beams (linearity check).
- Recommend the best combination and state plainly whether the sector scanner is worth building.

**D. Export of the recommended design (only if C recommends it).** Exact STEP + AEDT script + layout
(`hfss_exports/butler_sector_hex150_*`), README with: monopole positions, the ideal Butler matrix, the HFSS
**Edit Sources** amplitude/phase table for each Butler input and for ψ = 0/90/180/270° blends of one adjacent pair,
predicted beams/crosstalk to compare. Include an S-parameter note: the 4 monopoles now share one surface, so coupling
(S_mn) will be higher than on the separate-hexagon board — HFSS must report it.

**E. Close-out.** STATUS row (model-only), CHANGELOG, USER_GUIDE section if any UI was added (UI is optional this
phase — scripts + results are enough), results README, move this brief to `docs/history/briefs/`.

## Stop points
- After A: report the forward-model validation numbers before building on it.
- After C: give me the comparison table and your recommendation before doing D.

## Do not
- Change the solver core, existing synthesis results (rectangle/hex snapshot gates must stay identical), or claim any
  number not computed on the held-out grid. Don't model the Butler's own losses/imperfections beyond an ideal matrix
  unless asked (note it as a limitation).
