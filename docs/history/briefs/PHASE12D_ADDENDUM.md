# Phase 12D addendum — MoM mesh convergence at resonance (read after PHASE12D_PROMPT.md)

## New evidence (measured 2026-09-28, reproduce it)
HFSS vs MoM (40 cells, 0.2 mm, floquet 40, ref plane 5 mm), grounded 8 mm cell:
- G1 (L = 2.0): phase diff −0.13° → −0.52° over 10–15 GHz (agrees).
- G4 (L = 6.0, resonance): −4.6° to **−13.8°**, peaking at resonance (~12.25 GHz); |S11| within ~4e-4.
- G6 (L = 7.6): −9.6° at 10 GHz (near its resonance) → −0.8° at 15 GHz.
Pattern: HFSS resonates ~1% lower in frequency than the MoM.

MoM mesh refinement for G4 at 12.5 GHz (aligned meshes, floquet order = cells):
| cells | dx (mm) | S11_xx phase |
|---|---|---|
| 40 | 0.200 | −168.74° |
| 48 | 0.167 | −169.69° |
| 56 | 0.143 | −170.37° |
| 64 | 0.125 | −170.89° |
| 72 | 0.111 | −171.29° |
| 80 | 0.100 | −171.62° |
| 96 | 0.083 | −172.11° |
Fit φ = φ∞ + a·dx^p gives **p ≈ 0.96 (first-order convergence) and φ∞ ≈ −174.65°** — the 0.2 mm mesh is ~6° off at
resonance. HFSS gives −182.25° here; the extrapolated MoM closes the gap from 13.5° to ~7.6°. Floquet order 40→80 at
40 cells moved it only +0.3°.

## Tasks
1. **Convergence at resonance.** First-order convergence is the signature of unresolved edge-current singularities
   (rooftops on a uniform mesh). Implement and validate at least one of: (a) Richardson extrapolation in dx using
   the measured order (report φ∞ with an error bar from the fit); (b) a non-uniform mesh graded toward patch edges
   (keeping edges on mesh lines); (c) edge-singular basis functions. Show G4, G6 and one off-resonance case converge
   to a stated tolerance (e.g. 0.5°) and report cost.
2. **Re-state earlier numbers.** The Phase 12B/12C 8 mm results used 0.2 mm cells: quantify their resonance-region
   error (phase coverage, resonant L, and the Phase 12C/12D step-size errors, which compare embedded vs isolated
   on the same mesh and may partly cancel — check). Don't edit earlier outputs; write the corrections in 12D.
3. **HFSS comparison** (Task C of 12D): include G4 (`g4.s36p`) and G6 (`g6.s36p`) alongside G1, comparing HFSS to
   both the 0.2 mm MoM and the converged/extrapolated MoM. Treat the remaining ~7–8° at G4 as unexplained until
   tested; candidate causes: HFSS patch is 35 µm thick (MoM zero-thickness), HFSS mesh convergence at resonance.
   A pending user HFSS run of G4 with a zero-thickness Perfect E sheet will discriminate — ingest it if present
   (`g4_sheet.s36p`). Single-cell HFSS gate: validated off-resonance only unless the resonance gap is explained
   and within tolerance.

## Update — zero-thickness HFSS run of G4 (2026-09-28 17:13)
`hfss_results/g4.s36p` was **overwritten** by the zero-thickness Perfect E sheet run (patch at z = 0, no ThickenSheet;
same port/mode/mesh settings). The 35 µm run survives only as `G4_thick35um_comparison.csv` (per-frequency HFSS vs
MoM values); the sheet run's comparison is `G4_sheet_comparison.csv`. Treat `g4.s36p` as the sheet case.

Phase of S11_xx at the 5 mm plane:
| f (GHz) | MoM 0.2 mm | MoM extrapolated (dx→0) | HFSS sheet | HFSS 35 µm |
|---|---|---|---|---|
| 11.5 | −101.77 | −106.79 | −110.33 | −114.09 |
| 12.5 | −168.74 | −174.65 | −178.59 | −182.25 |
| 13.5 | +131.38 | +126.73 | +123.77 | +121.63 |
(Extrapolation: aligned meshes 48/64/80 cells, fitted φ∞ + a·dx^0.96; at 12.5 GHz also 40–96 cells.)

Budget at resonance (12.5 GHz, total 13.5°): MoM mesh ≈ 5.9°, patch thickness ≈ 3.7°, **residual ≈ 3–4°**
(consistent across 11.5–13.5 GHz, ≈ 0.5 % resonance-frequency shift). Candidates for the residual: HFSS mesh
convergence at resonance, uncertainty in the MoM extrapolation. Report it as such; the single-cell HFSS gate should
compare HFSS-sheet against the converged MoM with a stated tolerance and record thickness as a known model difference.
