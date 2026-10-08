# Nylon octagon cell library, 28 GHz (Phase 19 Part B) — model-only

Regenerate: `.venv/bin/python examples/build_nylon_cell_library_28ghz.py` (~3 min, 4 processes).
Solver: `src/dielectric_cell` (RCWA, normal-vector factorisation, M = 8). HFSS comparison:
`NOT_AVAILABLE_FROM_CURRENT_DATA`.

Fixed: pitch 5.0 mm, total thickness 25.0 mm (1.0 mm base plate + 24.0 mm pillar layer), nylon eps_r 2.45,
tan_d 0.0185, minimum feature 0.6 mm (octagon flat-to-flat and every gap/wall >= 0.6 mm, i.e. R in [0.325, 2.381] mm).
Two families: `pillar` (octagon in air on the plate) and `hole` (solid nylon, octagonal pocket open to the feed side).
58 cells; theta in {0, 10, 20, 25} deg, phi in {0, 45} deg (C4 symmetry folds every other azimuth into [0, 45]).

## Results (normal incidence, p polarisation)
- Phase span of the whole library: **455°**. Best contiguous 360° window: **369°, 44 cells, worst-case |S21| = −1.99 dB**
  (full-library loss range −0.64 … −2.39 dB; the 1 mm plate alone costs ≈ −0.64 dB through Fabry–Perot reflection).
- Pillars alone cannot reach 360° under the 0.6 mm rule (~284° with the plate): the largest pillar (R = 2.381 mm)
  is only 64 % fill. Solid-with-hole cells supply the high-fill end.
- Oblique incidence changes phase strongly (e.g. +49° at 22.6° for the R = 1.8 mm pillar, p-pol, phi = 0) and azimuth-dependently (up to ≈ 1.5° between phi = 0 and 45°);
  Part C therefore uses each cell's response at its own (theta, phi). Cross-pol |t_sp| ≤ 0.011 (−39 dB) at 22.6°.
- Pocket aspect ratio is up to 24 mm / 0.6 mm = 40 : 1 for the smallest holes — check your printer/nylon grade can resolve this.

`library.csv`: complex `t_*`, `r_*` (2×2, s/p basis, project convention exp(+jωt), reference planes in
`docs/theory/DIELECTRIC_CELL_RCWA.md`). `summary.json`: parameters and the numbers above.

## Through-hole vs blind pockets, and plate thickness (Phase 19 follow-up items 2-3)
Restricted to the best-360 deg window (fill <= 0.804 for the chosen 0.6 mm plate; smallest pocket used has a flat-to-flat width of ~2.4 mm, about 10:1 aspect).
Normal-incidence library at 0.6 mm plate (`library_plate0.6mm_{blind,through}.csv`): both span 462.9 deg with loss
-2.39 ... -0.31 dB; best 360 deg window worst-case |S21| -1.88 dB blind vs
-1.83 dB through (window spans 365.5 / 362.7 deg).

`plate_study.json` (`scripts/research/phase19_plate_study.py`; theta-interpolated ranking, not the final solve): total height 25 mm, illumination-weighted
aperture-averaged transmission for the dual-beam phase map, both families over their full range:

| plate [mm] | pockets | fill cap | aperture-avg loss [dB] | worst cell [dB] |
|---|---|---|---|---|
| 0.6 | blind | 0.804 | -0.852 | -1.91 |
| 0.6 | through | 0.804 | -0.836 | -1.84 |
| 1 | blind | 0.804 | -1.060 | -2.00 |
| 1 | through | 0.826 | -0.990 | -1.87 |
| 2 | blind | 0.826 | -1.239 | -2.11 |
| 2 | through | 0.847 | -1.130 | -1.91 |
| 3.42 | blind | 0.901 | -1.001 | -1.96 |
| 3.42 | through | 0.917 | -0.923 | -1.95 |

Pick: **0.6 mm plate, through-holes** (best aperture-averaged transmission and the only one that keeps the pocket fill cap at 0.804). 3.42 mm (half-wave) is second but
needs a cap of 0.92 (narrower pockets, harder to print). A 0.6 mm floor is thin; check it against your printer.
