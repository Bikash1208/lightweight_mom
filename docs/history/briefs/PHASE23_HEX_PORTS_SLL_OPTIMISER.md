# Phase 23 — hexagonal multi-port sub-apertures + sidelobe-cancelling cell-state optimiser

**Status:** active brief, written 2026-10-12. Read `docs/START_HERE.md`, `docs/STATUS.md`, `CLAUDE.md` first. Branch
`phase23-hex-sll` from **`phase22-auto-cell-states`** (the current tip: it has Phase 21's multi-port coded layouts and
exact export, and Phase 22's automatic cell states). All repo rules apply: one commit per item, fast suite green before
each, physics tests, honest labels, verify UI changes in a real browser, no solver-core changes (none are needed here).

## Why
The 4-port 1-bit boards (300 mm square sub-apertures, beams θ = 30°, φ = 45/135/225/315°, 10.54 GHz) predict co-pol
SLL ≈ −11 to −12 dB, set by the square, near-uniform sub-apertures and 1-bit quantisation. Two changes:
1. **Hexagonal sub-apertures** — tile with no gap; an idealised uniform comparison of equal-area apertures gives highest
   sidelobe square −14.1, hexagon −16.7, octagon −17.5, circle −17.6 dB; with a radial cosine taper square −19.1,
   hexagon −23.9 dB (idealised, ±1 dB; re-measure with the real model).
2. **Optimise the 1-bit cell states directly to cancel sidelobes** instead of only quantising the hologram
   (Johnson et al., "Sidelobe canceling for optimization of reconfigurable holographic metamaterial antenna",
   arXiv 1408.3032 / IEEE TAP 2015 — ref [66] in the user's literature review; related 1-bit pre-phasing: IEEE TAP
   2022, prephase-based equivalent amplitude tailoring). Include the case the user cares about: **all four ports
   excited together**.

## Work items

**A. Polygon (hexagonal) port regions.**
- Generalise `src/hologram/ports.py: PortRegion` from axis-aligned rectangles to a convex polygon (keep the rectangle
  constructor/fields working unchanged — existing rectangle layouts must reproduce **bit-identically**). Cell
  membership by cell-centre-in-polygon; area; bounding box for cropping in `multi_port_synthesis.port_native_grid`;
  edge-to-edge gap between polygons for `isolation_estimate`; `check_no_overlap` for polygons (touching edges allowed).
- New preset `hex4`: four regular hexagons in the compact cluster (two offset rows), parameter `side_mm` = hexagon side = circumradius
  (**default 150 mm**: flat-to-flat 259.8 mm, corner-to-corner 300 mm, area 58,457 mm² ≈ 0.65× the 300.8 mm square;
  also accept `flat_to_flat_mm` as an alternative input), optional `gap_mm` (default 0), orientation (flat-top / pointy-top), cluster centred at the origin,
  monopole at each hexagon centre. Ports ordered so default beams (θ = 30°, φ = 45/135/225/315°) point away from the
  cluster centre; document the order. Cells outside every hexagon are state 0 (empty) and the board outline is the
  cluster's bounding rectangle (simplest for HFSS) — make "trim board to the hexagon outline" an option.
- Export (Phase 21 exact STEP, AEDT script, DXF, CSV, PNG) works unchanged per cell; draw the hexagon outlines in the
  layout PNG and DXF (separate layer) and list them in the README.
- Tests: rectangle regions bit-identical; hexagon area/membership against analytic; cluster has zero overlap and zero
  gap at gap_mm = 0; isolation estimate gives 0 λsw for touching hexagons; a hex4 board synthesises end to end.

**B. Pattern evaluation that can be updated cell by cell.**
- The radiated field is linear in the per-cell aperture contributions (`coded_synthesis` stores
  `aperture_field = E_ref·exp(j·achieved phase)`; the vector model in `src/pattern/vector_aperture_pattern.py` is a DFT
  of that). Build, per port, the contribution of every cell to the co- and cross-pol far field on an optimisation grid
  of directions (u–v grid over the visible region, denser near beams), chunked to fit memory, so flipping one cell's
  state updates the pattern in O(directions). Verify against `compute_vector_pattern` / `_vector_far_field` to ~1e-10.
- Combined excitation: total field = Σ_ports a_p·e^{jα_p}·E_p (user-set amplitudes, phases optimisable).

**C. Sidelobe-cancelling optimiser.**
- Variables: the 1-bit state of every cell (start from the current coded design) and, for the all-ports mode, the
  four port phases α_p.
- Objective (per mode): minimise the peak co-pol sidelobe outside every requested beam's own main-lobe region (reuse
  `beam_polarization_metrics`' first-null exclusion and nearest-beam partition, so one beam's lobe is never counted as
  another's sidelobe), subject to: each beam's peak gain drops ≤ `max_gain_loss_db` (default 0.5 dB), pointing stays
  within 0.3°, beam levels stay balanced within 0.5 dB, XPD not worse by > 1 dB.
- Modes: (1) each port alone (single-port excitation, as HFSS Edit Sources one at a time); (2) all four ports together
  (joint over cells and α_p); (3) both at once (weighted worst case) — report which mode each design was optimised for.
- Method: deterministic and seeded; e.g. greedy/ranked single-cell flips with the O(directions) update, then simulated
  annealing or Johnson-style auxiliary-pattern cancellation as an option. No ML (repo rule).
- **Guard against over-fitting the grid:** report every number on a held-out fine grid (0.5° θ × 1° φ, the grid
  Phase 18 used) that the optimiser never saw, and at 10.54 GHz ± 2% (re-evaluate, don't re-optimise).
- App: a "Reduce sidelobes" option in Step 3 for coded multi-port designs: target mode, max gain loss, iteration
  budget; shows before/after SLL, gain, pointing per beam and the cells changed; export uses the optimised states.
- Tests: optimiser never increases the objective; constraints respected; a toy aperture where the optimum is known
  (e.g. small 1-D array) is recovered; held-out-grid numbers reported; reproducible with a fixed seed.

**D. Comparison runs (model-only; the user runs HFSS).** Same 10.54 GHz design, RT/duroid 5880 1.575 mm, 6.4 mm cells,
Phase 22's automatic 1-bit pair, RHCP beams θ = 30° at φ = 45/135/225/315°:
| Board | Notes |
|---|---|
| square 300.8 mm, gap 0 (existing Phase 21 reproduction) | baseline |
| square of equal area to the hexagon (≈ 241.8 mm, snapped to whole 6.4 mm cells), gap 0 | separates the shape effect from the area effect |
| hex4, side 150 mm, gap 0 | shape + 35% less area (expect ≈ −1.9 dB directivity vs the square; report it, don't hide it) |
| hex4 + optimiser, single-port mode | |
| hex4 + optimiser, all-ports mode (cells + port phases) | report the chosen α_p |
For each: directivity per beam, pointing, co-pol SLL (global and on the board x/y cuts through each beam, Phase 18
metric), XPD, level balance, cells changed — single-port and all-ports-excited. Export each optimised board's exact STEP
+ AEDT script to `hfss_exports/four_port_hex150_*` (STEP gitignored as usual) with a README giving the port phases to
set in HFSS Edit Sources for the all-ports case.

**E. Close-out.** `docs/STATUS.md` (capability rows: polygon ports, SLL optimiser — model-only), `CHANGELOG.md`,
`docs/USER_GUIDE.md` (how to use hex4 and "Reduce sidelobes"), screenshots, move this brief to `docs/history/briefs/`.

## Do not
- Change any existing rectangle-layout result (bit-identical gate).
- Claim an SLL improvement measured only on the optimisation grid.
- Touch the RIS (Phase 20) or transmitarray (Phase 19) work.
