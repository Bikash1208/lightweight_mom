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
