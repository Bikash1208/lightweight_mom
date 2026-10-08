# 28 GHz nylon dual-beam (+-30 deg) transmitarray lens -- model-only, awaiting your HFSS

Files: `lens.step` (ONE fused solid named **Lens**, mm), `lens.stl` (printing), `layout.csv` / `layout.json` (per-cell x, y, kind, octagon
circumradius, incidence angle, target phase, residual, power transmission), `verification.json` (STEP re-import checks).
Regenerate: `examples/export_nylon_lens_28ghz.py` (STEP/STL are gitignored).

## Geometry
Regular-octagon outline, 150 mm across flats (flats on the x and y axes), 5 mm cell pitch, 25 mm total thickness,
0.6 mm base plate; z = 0 is the **exit** face, z = 25 the **feed-side** face. 512 pillar cells (octagonal nylon pillar on the plate),
8 empty cells (plate only) and 236 hole cells (solid nylon with an octagonal through-hole cut through plate and lens).
Minimum printable feature 0.6 mm (your stated value). One material: nylon, eps_r 2.45, tan_d 0.0185.

## Horn
Phase centre on the axis, F = 180 mm (= 1.2 D) in front of the feed-side face, i.e. at **z = 205 mm in the STEP coordinates**
(the feed-side face is z = 25), pointing at the lens, **x-polarised** (E along x). Model horn: cos^q pattern, q = 14.38, 10 dB edge
taper at the lens edge. The beams point to theta = +-30 deg in the x-z plane (phi = 0 and 180 deg) on the exit side (z < 0).

## Model predictions to compare with HFSS (all model-only: HFSS comparison NOT_AVAILABLE_FROM_CURRENT_DATA)
- **ISOLATED-CELL prediction (deterministic)**: peak realised gain per beam **25.6 dBi** at theta = 29.87 deg (continuous peak search; the 0.25-deg grid of `summary.json` reports 30.00), SLL **-9.7 dB**, max cross-pol -19.3 dBi (off-plane). Horn power = 1; includes spillover, dissipation and the model's radiation efficiency.
- Aperture-averaged power transmission of the cells: 0.817 (-0.88 dB);
  fraction of horn power reaching the exit aperture: 0.752; radiated into propagating waves: 0.596.
- Per-cell phase residual vs requested: max 0.76 deg, rms 0.086 deg (cells above 1 deg: 0).
- **Inter-beam level difference 0.00 dB is an IDENTITY of the mirror-symmetric layout, not a performance result.** It cannot fail and must not be quoted as evidence.

## Error bar on the prediction (Monte Carlo of the measured local-periodicity error)
The isolated-cell figures above treat every cell as if it sat in an infinite array of itself. Exact supercell RCWA shows each cell's effective
exit field differs from that by a perturbation that grows with the phase step to its neighbours, and this layout needs steps with a median of
43 deg (90th percentile 172 deg) over all 2904 neighbour edges. It is
**set by the required aperture phase gradient**, so it is propagated, not just noted. 432 of 1452 internal edges (30 %) are pillar|hole, with a median step of 122 deg across them. Method: `scripts/research/phase19_error_bar.py`; 2000 draws per
scenario (seed 20261008); percentiles are stable -- between 250 and 2000 draws the 10/50/90 values move by at most 0.04 dB (gain), 0.14 dB (SLL), 0.05 deg (pointing).
E_c' = E_c (1 + eps_c) with eps_c = K * mean(T_n/T_c - 1) over the neighbours that differ (|T_n/T_c - 1| >= 0.6), K bootstrapped from exact-supercell
samples with the same number of differing neighbours (1 / 2 / 4 from strips, alternating strips, 2x2 checkerboards). **Single-neighbour coupling does not add
linearly** -- the 2- and 4-neighbour supercells give a per-neighbour coupling of 0.29 and 0.15 against 0.40 for one -- so naive superposition (scenario `Cn` in
`error_bar.json`) over-predicts and is not used.

Median (10th ... 90th percentile), isolated-cell values in the first row:

| model | mean gain / beam [dBi] | beam theta, phi=0 [deg] | worst-beam pointing error [deg] | SLL [dB] | abs inter-beam level difference [dB] |
|---|---|---|---|---|---|
| ISOLATED-CELL (deterministic) | 25.63 | 29.87 | 0.13 | -9.7 | 0.00 (identity, see below) |
| S: measured phase (primary) | 21.55 (21.12 ... 21.99) | 30.26 (30.00 ... 30.46) | 0.40 (0.24 ... 0.53) | -8.2 (-9.2 ... -7.1) | 0.41 (0.08 ... 0.94) |
| R: same magnitudes, random phase | 25.63 (25.29 ... 25.99) | 29.87 (29.67 ... 30.07) | 0.26 (0.13 ... 0.40) | -8.8 (-9.8 ... -7.8) | 0.39 (0.07 ... 0.90) |
| R0: order-0 relation, random phase | 25.63 (25.46 ... 25.78) | 29.87 (29.80 ... 29.93) | 0.13 (0.07 ... 0.20) | -9.3 (-9.9 ... -8.8) | 0.18 (0.03 ... 0.46) |

- **Phase assumption.** The supercell data DO constrain the perturbation phase (single-neighbour coupling is ~0.40 at -52 deg, not a random phase; two or four
  neighbours give a real ~0.5-0.65), so scenario S uses the measured complex K. Scenario R keeps the same magnitudes but randomises the phase: it leaves the
  coherent beam unchanged to first order and therefore shows **no** gain loss; the 4 dB between S and R is purely this phase assumption.
- **How far to trust S.** Checked on four real 9-cell windows of this layout solved exactly (`row_validation.json`, `row_models.json`): the isolated-cell order
  spectrum is off by 0.92, 1.42, 0.90, 1.44 (relative norm, i.e. as large as the signal); model S reduces that to 0.70, 1.36, 0.82, 1.32
  (random-phase R makes it worse: 1.25, 1.75, 1.21, 1.80). Total propagating power: exact/isolated = 1.03, 1.20, 1.05, 1.25, model S
  0.94, 0.88, 0.76, 0.91 -- S under-predicts power by roughly 0.4-1.4 dB there, so **S probably overstates the gain loss**; the true systematic
  effect is a gain reduction of unknown size between roughly 0 and the S value, and HFSS is the only thing that fixes it. The window test is at normal incidence,
  x-neighbours only.
- **Inter-beam level difference.** The nominal **0.00 dB is an identity of the mirror-symmetric layout** (the layout and the horn are symmetric under x -> -x, so the two
  beams are equal by construction). It cannot fail and is **not evidence** of anything. The Monte Carlo difference above is what a perturbation does to it.
- The error bar covers the local-periodicity error only. Not included: print tolerance, nylon eps_r spread, horn pattern/phase-centre error, truncation of the
  Huygens exit-aperture model, and the lens edge clipping effects.

## Cell-pitch lever (estimate only -- no re-design, no layout committed)
`scripts/research/phase19_pitch_estimate.py` runs only the phase synthesis on the finer grid. p = lambda/3 needs 42 cells across (D/p must be an even integer),
so **p = 150/42 = 3.5714285... mm, not an exact decimal** (current p = 5 mm = 0.467 lambda).

| | p = 5 mm (now) | p = 3.5714 mm (lambda/3) |
|---|---|---|
| cells | 756 | 1452 (1.92x) |
| minimum feature | 0.6 mm | 0.6 mm (unchanged); octagon R 0.32 ... 1.61 mm (was ... 2.38) |
| smallest pocket (fill cap 0.804) | 2.43 mm wide, 10:1 | 1.74 mm wide, 14:1 |
| median neighbour phase step (target) | 43 deg | 27 deg (0.62x) |
| 90th-percentile step; edges over 90 deg | 172 deg; 27 % | 75 deg; 0 % |
| coupling-error rms of eps (measured K, equal-amplitude targets) | 0.62 | 0.36 (0.57x) |

Expected error reduction from the measured relation: about **43 %** in the rms perturbation (the error is driven by the phase step). The cost is
2x the cells and narrower pockets (14:1); the pillar family shrinks (largest fill 0.57), which the library would have to be re-checked against.

## STEP coordinates are exact (Phase 19 Part D-bis)
- **Why**: the OpenCASCADE writer prints reals with 12 significant digits and its boolean kernel leaves 1e-16..1e-10 noise, so the first export carried literals such as
  `8.500000000001E-02` (intended 0.085), `6.169999999999E-02` (0.0617) and `1.110223024625E-16` (0) -- 34,764 of them. The earlier "no float noise" check had a blind spot
  below 1e-13 and missed them.
- **Fix**: after the kernel writes, every `CARTESIAN_POINT`, `DIRECTION` and `VECTOR` magnitude is snapped to its intended design value when within 1e-7 mm, then re-written
  with the shortest decimal that round-trips to the same double (`4.0` -> `4.`, `0.085` -> `0.085`; never a fixed `%.Nf`). Intended values come from the layout table
  (cell centres, +-pitch/2 corners, plate and total thickness, outline, octagon vertices rounded to 1e-4 mm, outline-flat x grid intersections), never from kernel arithmetic.
- **Result on this lens**: 409,872 coordinates walked, 34,764 inexact before, **0 after**; largest snap displacement
  4.5e-13 mm; no two previously distinct vertices merged (12,436 vertices, 0 violations,
  0 collapsed edges); one valid solid, self-intersection analyzer clean; volume before vs after differs by 1.7e-16 (limit 1e-9).
- **Uncertainty**: `UNCERTAINTY_MEASURE_WITH_UNIT` (it lives in the STEP data section, in the geometric representation context, not the header) is set to
  **1.E-6 mm** (OpenCASCADE's default is 1e-7). It must exceed the residual coordinate noise (<= 1e-9 mm after this export) and the snap tolerance (1e-7) by an order of
  magnitude so an importer welds duplicate vertices, and stay far below the smallest real vertex separation (0.0617 mm, a factor
  61,700) so it never merges distinct ones.
- **Cannot be exact in decimal (do not "fix")**: octagon outline vertices (75 tan 22.5 deg = 31.0660171779821...; stored as 31.066, a 1.7e-5 mm deliberate rounding);
  all pillar/pocket octagon vertex coordinates (cos/sin of 22.5 deg multiples times a solved radius; rounded to 1e-4 mm, up to 5e-5 mm); the solved radii themselves; the edge
  directions of the rounded octagons and the 1/sqrt(2) directions of the outline flats; intersections of octagon edges with the 45-degree outline flats. The pitch (150/30 = 5 mm)
  and the grid are exact; **a pitch such as 150/42 (lambda/3) would not be**.
- **STL**: an STL is a triangulated approximation; exactness is not a meaningful property of it, so none is claimed or chased there.
- **Short edges**: 1928 vertex pairs lie closer than the 0.6 mm minimum feature (smallest 0.0617 mm), where pillar/pocket octagons are cut by the diagonal outline flats. Review these
  before printing; this export does not change the geometry.

## Suggested HFSS setup
Material nylon (eps_r 2.45, tan_d 0.0185); horn at the phase centre above; radiation (or PML) boundary >= lambda/4 (2.7 mm)
beyond the structure on all sides; far-field sphere, Ludwig-3 co/cross for x polarisation; fine mesh in the nylon (lambda/(10 n) ~ 0.7 mm).
The lens has 756 cells; symmetry planes (x and y) can be used for the nominal geometry.

## Verification history
- Original Phase 19 series (Parts A-D): the fast test suite was run **once on the final tree** and not at each intermediate commit, so commit `af7b0aa` (the supercell-capable
  RCWA extension, the only follow-up commit that modified existing solver files) is verified by its outcome -- it reproduces the committed library to CSV rounding and every later
  test passed on top of it -- rather than individually. The Part D-bis commits (error bar, exact STEP coordinates) each ran the full fast suite before committing.
- Part A S21 reference plane: the "fix" was to a wrong *statement* (an empty cell was documented as S21 = 1) in the docstring and the tests, made before the first Phase 19 commit; the
  solver always returned S21 = exp(-j k0 cos(theta) h) for an empty cell. Generated after the first Phase 19 commit (20bf86d, 2026-10-08 13:01) that contains the correction: cell library 2026-10-08 14:53,
  layout and prediction 2026-10-08 15:28, this STEP 2026-10-08 18:40 -- so the correction is not just in the docs and tests.
