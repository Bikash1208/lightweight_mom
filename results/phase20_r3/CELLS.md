# Phase 20 R3 item 2 — the three starting cells (NOT designs)

**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.** Code: `src/ris/cells.py` (`Lp1Cell`, `Lp2Cell`, `Cp2Cell`, diode kits SMP1345 / MADP-000907 /
ideal), `src/polygon_geometry.py`, rotated loads in `src/loads.py`, circular basis in `src/ris/polarization.py` (convention: `docs/theory/CONVENTIONS.md`). Data: `cells_start.json`
(`scripts/research/phase20_r3_cells_start.py`), 10 GHz, RT5880 1.575 mm, PEC ground, mesh 48 (LP) / 54 (CP; 0.25 mm cells), normal and 30° incidence. All dimensions are starting values that
item 3 optimises. LP-1 and LP-2 defaults are the best of 70-point random scans for four distinct states only.

| cell | what | period | diodes | states | source |
|---|---|---|---|---|---|
| **LP-1** | three-section patch, a diode bridge across each gap (strips 0.5 and 2.0 mm wide) | 12 mm | 2 | (d1, d2) | roadmap baseline; literature A, B |
| **LP-2** | **new single-layer proposal of this project (no published precedent)**: patch 6.0 × 3.5 mm with two diode-switched stubs of different length | 12 mm | 2 | (d1, d2) | decision 2026-10-19 (LP must be 2-bit) |
| **CP-2** | regular octagon (circumradius 4.5 mm) + four radial stubs at 0/45/90/135°, one diode each, one-hot state; centre via omitted | 13.5 mm | 4 | which stub is ON | literature G (`PAPER_NOTES.md`) |

## Where they stand (10 GHz, normal incidence)
* **LP-1 (SMP1345):** TM phases 172 / 154 / 145 / 161° with |Γ| 0.88–1.0 — four *distinct* states but a ~27° spread: far from 0/90/180/270°. TE (the other polarisation) is unchanged (118° in all states).
* **LP-2 (SMP1345):** TM phases 101 / 90 / 72 / 45° (a ~55° spread, |Γ| ≥ 0.945); TE 137° in all states. Also far from 2-bit.
* **CP-2:** the handedness-preserving term |r_RR| is only 0.07–0.20 (an ideal half-wave plate gives 1): the starting stubs give a ~25° x′/y′ phase difference, not the 180° of the published element, so the cell is a **weak** PB element.
  With **ideal switches** the four states are rotations of one another: r_RR phases −6 / 85 / 175 / −94° (steps 91, 90, 91, 88°, staircasing error under ±3°). With **MADP-000907** (L_s = 0 assumed, C_T 0.025 pF) the ladder nearly survives:
  −23 / 68 / 167 / −101° (steps 91, 99, 92, 78°). With **SMP1345** it does **not**: −67 / −5 / −168 / −109° (steps 62, 197, 59, 106°) — the OFF diode is not an open at 10 GHz (0.19 pF in series with 0.7 nH is −40 Ω),
  so the three OFF stubs load the patch and the one-hot states stop being rotated copies. Consequence for item 3: the SMP1345 CP-2 cell needs OFF stubs that are insensitive to the OFF impedance (or a topology change), the MADP cell is the closer to the published element; both enter the tolerance objective.

## Modelling findings of this item
* A rotated (diagonal) diode region must **conduct transversely** (`transverse_conducts=True`): on a 45° strip a few cells wide the transverse block forces J_x = J_y on a staircase and kills the longitudinal current too (tested; MOM_FORMULATION Sec. 15).
* Rotated-load integrals use a composite Gauss rule (convergent, not exact); uniform-current dissipation is exact to ~1e-12 (conducting) / ~1e-6 (blocked, penalty rounding); mirror symmetry of the cell response holds to 1e-9.
* Polygon patches are staircased; `edges_aligned` is False for any non-axis-aligned edge. The 45° states are therefore rotated copies only to within the staircasing (about ±3° of r_RR phase here).
* The circular-basis convention differs from Phase 5's (propagation-relative handedness); `r_LR` of a plate is −1, a rotated half-wave plate gives `r_RR = exp(+j2α)` (tested).
* Not done in this item: CP-3 (four-slot circle, the alternative), the HFSS chain/pad representation of the CP-2 diodes (item 5), any optimisation (item 3).
