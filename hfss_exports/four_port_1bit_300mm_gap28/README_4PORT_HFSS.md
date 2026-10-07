# 4-port 1-bit coded HMSA, 300 mm sub-apertures — HFSS check (gap 28 mm)

Minimum-gap study, board B of 2 (the other is `four_port_1bit_300mm_gap0/`). Same design as
`four_port_1bit_gap28/` but with each sub-aperture enlarged from 142.2 mm (5 λ0) to **300.8 mm (10.6 λ0, 47 × 47 cells)**.

10.54 GHz · Rogers RT/duroid 5880 1.575 mm · PEC ground · 6.4 mm cells · board **629.6 × 629.6 mm** (98 × 98 cells), centred at (0, 0).
Four 300.8 mm sub-apertures in a 2×2, **28 mm apart**; one monopole at each sub-aperture centre.

| Port | Monopole (x, y) mm | Beam (θ, φ) | Polarisation |
|---|---|---|---|
| 1 | (164.4, 164.4) | (30°, 45°) | RHCP |
| 2 | (-164.4, 164.4) | (30°, 135°) | RHCP |
| 3 | (-164.4, -164.4) | (30°, 225°) | RHCP |
| 4 | (164.4, -164.4) | (30°, 315°) | RHCP |

Cells: "0" = empty (bare substrate), "1" = 5.56 mm square copper patch with an annular slot (inner radius 1.179 mm,
outer radius 1.922 mm). Patches are zero-thickness sheets in one body "Patches" → Perfect E; "Ground" → Perfect E;
"Substrate" → RT/duroid 5880. `array.step` is gitignored (regenerate with the command below); everything else is committed.

## Model prediction (this project's vector aperture model, not HFSS)

Per port, fine 0.5° × 1° far-field grid (`design_report.json` → `ports[*].fine_grid`). SLL is the highest co-pol
(RHCP) sidelobe on the cuts through the beam peak along the board's x and y directions (`port_cuts.png`) — a square
aperture's first sidelobes sit there; the diagonal (φ = beam φ) plane is much lower (≈ −23 dB).

| Port | Directivity (dBi) | Peak (θ, φ) | Pointing error | SLL along x / y (dB) | XPD at peak (dB) | Aperture eff. (upper bound) |
|---|---|---|---|---|---|---|
| 1 | 27.6 | (29.5°, 45.0°) | 0.5° | -12.0 / -11.5 | 17.0 | 0.40 |
| 2 | 27.6 | (30.0°, 135.0°) | 0.0° | -11.9 / -11.8 | 16.9 | 0.37 |
| 3 | 27.6 | (29.5°, 225.0°) | 0.5° | -12.1 / -11.4 | 17.2 | 0.36 |
| 4 | 27.6 | (30.0°, 315.0°) | 0.0° | -11.8 / -11.8 | 17.2 | 0.40 |

Isolation (model): Not flagged: gap = 1.00 λsw (λsw = 27.5 mm). The model has **no coupling prediction** — S21/S31/S41 need HFSS.

Reading these: the ≈ −11 to −12 dB SLL is set by the near-uniform illumination of a square sub-aperture plus 1-bit
quantisation, not by its size — enlarging 142 → 301 mm narrows the beam (21.4 → 27.6 dBi) but lowers the SLL by
only ≈ 0.5–1 dB (same script, 142.2 mm boards on the same fine grid: −10.0 to −11.5 dB). Lower SLL would need amplitude tapering, which a phase-only 1-bit hologram does not provide directly. The
cross-pol (LHCP) shoulders at ≈ −6 dB beside the beam are the documented scalar-cell polarization ceiling
(`docs/STATUS.md` capability matrix), not a size or gap effect.

## HFSS

**Before the full board:** one cell of each state (Floquet port, PEC ground, 6.4 mm period) should be ~180° apart at
10.54 GHz (MoM: +138° and −57° at the patch plane, 165° apart). The slotted cell sits on a steep resonance; if HFSS
disagrees, adjust the patch size b first.

**Full board:** same monopole as your reference-wave model at each position above; radiation boundary ≥ 15 mm around
and 30 mm above; 10.54 GHz. Solve the 4×4 S-matrix (S21 = adjacent along x, S41 = adjacent along y, S31 = diagonal)
and the far field one port at a time (Edit Sources). The board is ≈ 22 λ0 square — expect a large mesh;
FE-BI or the IE solver (zero-thickness sheets on an infinite-ground approximation) may be the practical choice.
Compare with the other gap version: the gap-study question is how much S21/S41 improve from 0 → 28 mm, and whether
the far-field SLL/pointing above survive when the neighbours are present.

Regenerate: `python scripts/adhoc/coded4.py hfss_exports/four_port_1bit_300mm_gap28 28 rhcp 30 model 300.8`
