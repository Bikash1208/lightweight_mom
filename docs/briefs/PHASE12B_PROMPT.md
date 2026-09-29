# Phase 12B brief — mesh alignment, duplicate-geometry audit, grounded period study

You are working in `lightweight_mom`. Read `CLAUDE.md`, `CONVENTIONS.md`, `PHASE12A_GROUND_PLANE.md`,
`results/phase12a_ground_plane/phase12a_report.md` and `UNIT_CELL_LIBRARY.md` first.

## Why this phase exists (measured, reproduce before relying on it)

1. **Mesh/geometry misalignment.** A mesh cell is metal iff its centroid is metal (`CONVENTIONS.md`).
   Phase 7 and Phase 12A use 32 cells over a 4.8 mm period (0.15 mm), but patch sizes step by 0.2 mm,
   so adjacent sizes often select the same cells:
   - Phase 12A `sweep.csv`: L = 1.4/1.6, 2.0/2.2, 2.6/2.8, 3.2/3.4, 3.8/4.0 mm give identical S11 at every
     frequency.
   - Phase 7 `results/unit_cell_library/database/phase_lookup_12p5GHz.csv`: 324 geometries but only
     **169 distinct complex responses** (105 duplicate groups).
2. **Phase 12A's `mesh_error = 0.048` is probably a staircasing artefact, not under-resolution.** Its
   refinement study (L = 2.8 mm, 8/16/32/64 cells) never has a mesh edge on the patch edge; the phase
   oscillates (−26.4°, −35.8°, −31.0°, −33.8°) instead of converging monotonically. Test this hypothesis;
   don't assume it.
3. **The 4.8 mm period (~0.2 λ at 12.5 GHz) is too small for a grounded patch to traverse resonance.**
   At 12.5 GHz the phase moves −21° → −85° over L = 1.0–4.0 mm and the last two points (4.2, 4.4 mm; 0.4 mm
   gap) supply most of the 203.6°. Coverage at 10 GHz is only ~86°.

## Tasks

### A. Aligned meshes (core capability, no physics change)
- Add a helper (e.g. in `discretization.py`) that, for a lattice and a set of rectangle dimensions and
  centres, reports whether every feature edge lies on a mesh line, and proposes the smallest aligned
  `cells_x/cells_y` (e.g. for P = 4.8 mm, 0.2 mm steps centred at 0 → 48 cells, 0.1 mm). Edges sit at
  `−P/2 + k·P/N`; derive the condition, don't guess it, and test it.
- Add a **duplicate-geometry detector**: two configs whose metal masks are identical on the chosen mesh
  must be flagged (warning + record flag) before solving. Plan-time sweeps (`sweep.plan_sweep`) should
  report how many requested geometries are mesh-distinct.
- Additive only: existing configs and results must stay bit-identical; do not change the centroid rule.

### B. Re-assess Phase 12A with aligned meshes
- Redo the grounded refinement study for L = 2.8 mm on aligned meshes (e.g. 24 / 48 / 96 cells) with
  Floquet order scaled per the existing rule and accuracy guidance; also check Floquet convergence
  separately. Report whether the error now converges monotonically and meets the 0.02 tolerance.
- Re-run the Phase 12A sweep (same substrate, 4.8 mm, L = W = 1.0–4.4 mm, 10–15 GHz) on the aligned mesh,
  plus finer L steps near resonance (e.g. 0.1 mm for L ≥ 3.6 mm, on a mesh aligned to them). Report
  coverage, largest gap and |S11| next to the Phase 12A numbers.
- Update the Phase 12A status only if the evidence supports it: write a new gate record
  (e.g. `GROUNDED_STACK_VALIDATED` or still `GROUNDED_STACK_PARTIAL` with the reason). Do not edit the
  Phase 12A outputs in place; write to `results/phase12b_*`.

### C. Duplicate-geometry audit of earlier phases (read-only)
- Quantify the duplicates in the Phase 7 library (per frequency, via the metal masks, not just equal S
  values) and list them.
- Trace which geometries Phases 8C–11 actually used (Phase 9 contrast pairs, Phase 10 dataset
  train/validation/blind splits and held-out centre, Phase 11's 5 centre geometries and patterns) and
  report whether any "distinct" geometries were mesh-identical — in particular any train/test leakage
  in Phase 10. For each affected finding, state whether its conclusion changes, and why. Do not rebuild
  the Phase 7 library or re-run earlier phases; if a conclusion is actually at risk, say so and propose
  the minimal re-run.

### D. Grounded period study
- Same substrate (RT/Duroid 5880, 1.575 mm), grounded, rectangular PEC patch L = W, aligned meshes.
  Periods: 4.8 mm (baseline from B), 8.0 mm, and ~12 mm (≈ λ/2 at 12.5 GHz). L from ~0.2 P up to P minus a
  documented minimum gap; step small enough to resolve the resonance (refine near it).
- For each period, 10–15 GHz: phase vs L curves, phase coverage and largest gap at 12.5 GHz, coverage vs
  frequency, |S11| (loss) vs L, and the L at which resonance occurs.
- Checks per period: grating-order onset (`discretization.py`, including the oblique-incidence limit
  `P < λ/(1 + sin θ)` for a documented set of θ — report, don't assume a feed angle), and the
  surface-wave-pole margin from Phase 12A (it shrinks as P grows; for a 12 mm period at 15 GHz the lowest
  harmonic sits much closer to the TM0 band than at 4.8 mm — measure it). Mesh/Floquet convergence for one
  near-resonant patch per period.
- Report measured coverage honestly. A single-layer resonant patch typically cannot reach a full 360°;
  if it doesn't, say so and list the standard remedies (thinner substrate, stacked patches,
  multi-resonant elements) as open options — do not implement them here.
- Register nothing as production data; any new sweep is a benchmark, tagged grounded + period + mesh.

## Deliverables
- `examples/run_phase12b.py` → `results/phase12b_mesh_period_study/` (CSV/JSON, plots: aligned vs
  unaligned convergence, phase vs L per period and frequency, coverage vs frequency per period, |S11| vs L,
  duplicate-audit table). Deterministic.
- `PHASE12B_MESH_PERIOD_STUDY.md`, and updates to `README.md` (status, layout, "Still open"),
  `UNIT_CELL_LIBRARY.md` (duplicate finding), `CONVENTIONS.md` (alignment rule), `CLAUDE.md` if commands
  change.
- Explicit statuses from fixed sets, e.g. mesh alignment `ALIGNED_MESH_VALIDATED`; audit
  `NO_CONCLUSION_CHANGED` / `CONCLUSIONS_AT_RISK` (listed); period study `COVERAGE_SUFFICIENT` /
  `COVERAGE_INSUFFICIENT` with the numbers. HFSS status unchanged; geometry inversion and STEP export
  stay LOCKED.

## Rules
- Measured numbers only; no fabricated or digitised data without recorded provenance.
- Don't change earlier results in place, the ungrounded physics, or the centroid rule.
- Don't choose a period or geometry by phase coverage alone — report coverage alongside grating,
  surface-wave, loss and bandwidth findings and leave the choice to the user.
- Keep new tests property-based; run the full suite at the end and report counts and runtime.
