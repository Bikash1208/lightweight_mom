# 5.8 GHz meander cell screen (MODEL_ONLY)

RT5880 εr 2.2, tanδ 0.0009, 3.175 mm over PEC; period 11.6 mm, disc radius 5.6 mm. The straight slot has a 3.0 × 0.4 mm cut. The proposed baseline has 3 turns, a 1.2 mm copper finger, 0.4 mm slot, and 3.0 mm full horizontal run along x. This full run is how `tooth_depth` was interpreted. MoM uses a uniform rooftop mesh with at least three cells across the narrowest gap (88 cells/axis for 0.4 mm; 0.2 mm would require at least 174 and was not solved). The β values use the normal-incidence tensor-sheet approximation, not a Bloch eigenmode. No HFSS cell solve has checked these rankings.

| Cell | β/k₀ x | β/k₀ y | Δβ vs straight x/y | M_t | α rad potential (Np/m) | α diel mean (Np/m) | rad/diel | Mesh |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline_n3_t1p2_s0p4 | 1.0767 | 1.0655 | +0.0139/+0.0021 | 0.0641 | 0.0558 | 0.0226 | 2.47 | 88² |
| straight_w0p4 | 1.0627 | 1.0634 | +0.0000/+0.0000 | 0.0042 | 0.0002 | 0.0195 | 0.01 | 88² |
| meander_n2_t1p2_s0p4 | 1.0671 | 1.0650 | +0.0043/+0.0015 | 0.0124 | 0.0021 | 0.0206 | 0.10 | 88² |
| meander_n4_t1p2_s0p4 | 1.1142 | 1.0659 | +0.0515/+0.0025 | 0.2647 | 0.7046 | 0.0291 | 24.21 | 88² |
| meander_n3_t0p8_s0p4 | 1.0685 | 1.0652 | +0.0058/+0.0017 | 0.0196 | 0.0053 | 0.0209 | 0.25 | 88² |
| meander_n3_t1p2_s0p6 | 1.0861 | 1.0658 | +0.0233/+0.0024 | 0.1145 | 0.1693 | 0.0244 | 6.94 | 88² |
| meander_n3_t1p2_s0p4_rot90 | 1.0655 | 1.0767 | +0.0028/+0.0132 | 0.0641 | 0.0558 | 0.0226 | 2.47 | 88² |
| island_n3_t1p2_s0p4 | 1.0736 | 1.0656 | +0.0109/+0.0021 | 0.0467 | 0.0298 | 0.0220 | 1.35 | 88² |

M_t = |X_x - X_y|/|X_x + X_y|, where X is the effective surface reactance. The α rad value is a first-order **30° rotating-axis modulation potential** inferred from the local tensor. For the infinite uniform bound-mode cell, the modeled radiating leakage is zero; this nonzero potential requires spatial hologram modulation. The dielectric α comes from the lossy-slab transverse-resonance root; copper is PEC, so conductor loss is omitted. The ratio is a potential, not a measured radiation efficiency. The large four-tooth α is particularly sensitive to the first-order approximation. At p=11.6 mm, spatial dispersion and three-cell-wide gaps need HFSS and mesh-convergence checks.

The removed-tooth island row is geometry only: no diode, bias or switching claim. Relative to the static baseline, its β/k₀ shifts by -0.0031 along x and +0.0001 across y in this model; this does not establish a radiating diode state. [Finalists](finalists.json) list three unique static parameter sets, selected by M_t then model loss ratio. The [ideal aperture estimate](../aperture/REPORT.md) is separate and does not use these cell responses.

## HANDOFF

- STATUS: Cell screen complete, MODEL_ONLY.
- Cell JSON: `finalists.json` (2–3 static cells); removed-tooth row is not a finalist.
- β/k₀ along/across, α, M_t, loss ratio: see table and per-cell JSON.
- HFSS or MODEL_ONLY: MODEL_ONLY; HFSS cell and array checks pending.
- Next request: HFSS eigenmode/periodic-cell check of the finalist β and loss, then a physical CP/leakage board study if the cell passes.
