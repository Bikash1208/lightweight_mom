# SIDEQUEST_CP_HYBRID - CP option (i): RHCP beam 1 (tensor cells) + LP beam 2 (GB cells), model + HFSS package (2026-11-06). MODEL_ONLY
Package: `hfss_exports/sidequest_board_CP_LP_hybrid/` (STEP 38 MB gitignored, synced to `C:\hmsa_sync\data\sidequest\board_CP_LP_hybrid`; `HFSS_SETUP.md`, `predictions.json`, `cells.csv`, `sheets.csv`, `states.csv`, `geometry_checks.json` all PASS, EXACT_DIMS PASS). Scripts `sidequest_cp_hybrid.py` (model, placement, joint optimisation, sweep), `sidequest_cp_hybrid_package.py`, `sidequest_cp_hybrid_docs.py`.
## Model
Phase 26 tensor engine for the CP part (ray-march over the SOLVER_SUPERCELL table, vector far field, gain = D_ap x launch 0.592 x eta_conv x anchor -0.48 dB from the HFSS 23.03 dBi board), plus the 180 GB crossed-slit cells as sheet-current elements with the SWHB Born amplitude (same radiated-power-per-area units as the tensor `amp`), local mean sheet from the tensor beta. All cells except the GB sites are tensor cells (the former GA region and the static open/strap cells are replaced). Tensor beta/k0 at the cells 1.051 / 1.062 / 1.073 (min / mean / max) against 1.056-1.063 for the GB/static cells of the calibrated A/B model: **no wave-speed jump**. Relative sign of tensor and GB currents: physical argument (a slit across E raises the sheet susceptance: K2 library open cross slit 8.69 mS vs solid disc 5.18 mS) gives -1; the opposite sign is carried as a falsifier.
## Results (rel -1; bands: beam 1 +-1 dB, beam 2 +-3 dB, SLL +-3 dB incl. the +2.5 dB offset)
| | literal replacement (GB cells as in the A/B board) | **jointly optimised (package)** |
|---|---|---|
| beam 1 RHCP, GB OFF / ON [dBi] | 22.5 / 22.5 | **22.4 / 22.2** (AR 2.1 / 2.1 dB, worst in the -3 dB beam 2.8 / 2.9, XPD 18 dB, HPBW 7.0 deg; Phase 26 HFSS AR was 1.2-1.3 where the model said 2.5) |
| beam 2 at (25, 180), GB OFF / ON [dBi] | 1.7 / 4.5 | **-3.5 / 7.7** (LP-V: tilt 6.5 deg from phi-hat, phi-hat / theta-hat 18 dB, HPBW 9.9 deg) |
| beam-2 contrast / balance | 2.8 dB / 18 dB | **11.2 dB / 14.5 dB** |
| SLL (with +2.5 dB) | -15 | -11.7 (OFF) / -11.5 (ON) |
Compared with the current LP board (HFSS 13.2 dBi, contrast 6.5 dB, SLL about -5): beam 1 +9 dB, circular, contrast +4.7 dB, SLL +6 dB better, beam 2 only 7.7 dBi (A/B board GA = PEC: 8.8) and balance 14.5 dB (targets: contrast >= 12 dB nearly met, balance <= 3 dB NOT met, SLL <= -13 NOT met, beam 1 >= 17 met).
Why: the GB cells alone radiate about 3 dBi into beam 2 (a 2.6 mS step on 12 % of the cells against tensor cells with alpha up to 4.3 Np/m); the 8 dBi at GB ON comes from the coherent doubling with the tensor background at (25, 180), which the joint optimisation sets to cancel the GB OFF field: 102 of 180 GB sites relocated, 50 tensor slots rotated by 90 degrees (sign flip of J), beam-1 amplitude loss <= 2 dB budget (0.1 dB realised). Balance <= 3 dB is out of reach with these cells: it needs a stronger switched cell (step 3).
## Diode reduction / GB placement (model; rel -1, loss budget 0.4 dB)
| GB cells / diodes | beam 2 ON / OFF [dBi] | contrast | balance | beam 1 AR OFF |
|---|---|---|---|---|
| 180 / 360 | 8.0 / -1.1 | 9.2 | 14.1 | 2.5 |
| 150 / 300 | 6.7 / -2.8 | 9.5 | 15.7 | 2.4 |
| 120 / 240 | 5.1 / -6.0 | 11.1 | 17.4 | 2.4 |
| 90 / 180 | 3.3 / -6.4 | 9.7 | 19.4 | 2.3 |
Contrast is flat (9-11 dB) down to 180 diodes; the beam-2 level falls 1.3 dB (300), 2.9 dB (240), 4.7 dB (180). The current GB regions are NOT optimal for a tensor background: the single-cell coherence ranking alone does not help (the background interference dominates); the joint site + slot-rotation optimisation does (contrast 2.8 -> 9.2-11.2 dB at 360 diodes).
## Risks
1. **Relative sign** of the tensor and GB cells (physical argument only): if wrong the design inverts (beam 2 +7.9 dBi OFF, -1.3 dBi ON). First HFSS check: compare S_w 1e6 vs 6.1 at (25, 180).
2. The GB cell fields carry the 2.5 dB contrast optimism seen on the A/B board (contrast 9 vs HFSS 6.5): expect 7-9 dB instead of 11.2.
3. Mutual coupling between rotated tensor slots and GB crossed-slit cells is not modelled; AR may differ either way.
## CP feasibility notes (not part of this package)
(ii) Switchable CP would need a cell whose two orthogonal axes are independently switchable (e.g. two diodes per cell on the slot axes with opposite-handed 90-degree rotated states, one line per handedness or per axis) with a susceptance step at least 5x the GB cell's 2.6 mS so that its radiated power per area is comparable to the tensor cells (alpha 4 Np/m); with one line only the cell must flip the sign of J (slot rotation by 90 degrees is the sign flip: a pair of diodes that selects between two orthogonally oriented slots).
