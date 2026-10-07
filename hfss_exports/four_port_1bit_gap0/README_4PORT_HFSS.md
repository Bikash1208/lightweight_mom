# 4-port 1-bit coded HMSA — HFSS check (gap 0 mm)

10.54 GHz · Rogers RT/duroid 5880 1.575 mm · PEC ground · 6.4 mm cells · board 284.4 × 284.4 mm, centred at (0, 0).
Four 142.2 mm (5 λ0) sub-apertures in a 2×2, 0 mm apart; one monopole at each sub-aperture centre.

| Port | Monopole (x, y) mm | Beam (θ, φ) | Polarisation |
|---|---|---|---|
| 1 | (71.1, 71.1) | (30°, 45°) | RHCP |
| 2 | (−71.1, 71.1) | (30°, 135°) | RHCP |
| 3 | (−71.1, −71.1) | (30°, 225°) | RHCP |
| 4 | (71.1, −71.1) | (30°, 315°) | RHCP |

Cells: "0" = empty (bare substrate), "1" = 5.56 mm square copper patch with an annular slot (inner radius 1.179 mm,
outer radius 1.922 mm). Patches are zero-thickness sheets in one body "Patches" → Perfect E; "Ground" → Perfect E;
"Substrate" → RT/duroid 5880.

**Before the full board:** simulate one cell of each state (Floquet port, PEC ground, 6.4 mm period) and confirm the
two reflection phases are ~180° apart at 10.54 GHz (MoM: +138° and −57° at the patch plane, 165° apart). The slotted
cell sits on a steep resonance; if HFSS disagrees, adjust the patch size b first.

**Full board:** same monopole as your reference-wave model at each position above; radiation boundary ≥ 15 mm
around and 30 mm above; 10.54 GHz. Solve the 4×4 S-matrix (S21, S31, S41 = isolation) and the far field one port at a
time (Edit Sources). Compare this board with the other gap version: the prediction is that S21/S41 differ by only
~1–2 dB, i.e. a bare gap is not what provides isolation.
