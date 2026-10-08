# 28 GHz nylon dual-beam (+-30 deg) transmitarray lens -- model-only, awaiting your HFSS

Files: `lens.step` (ONE fused solid named **Lens**, mm), `lens.stl` (printing), `layout.csv` / `layout.json` (per-cell x, y, kind, octagon
circumradius, incidence angle, target phase, residual, power transmission), `verification.json` (STEP re-import checks).
Regenerate: `examples/export_nylon_lens_28ghz.py` (STEP/STL are gitignored).

## Geometry
Regular-octagon outline, 150 mm across flats (flats on the x and y axes), 5 mm cell pitch, 25 mm total thickness,
0.6 mm base plate; z = 0 is the **exit** face, z = 25 the **feed-side** face. 520 pillar
cells (octagonal nylon pillar on the plate) and 236 hole cells (solid nylon with an octagonal
through-hole cut through plate and lens). Min feature 0.6 mm. One material: nylon, eps_r 2.45, tan_d 0.0185.

## Horn
Phase centre on the axis at **z = +180 mm** (F = 1.2 D) above the feed-side face, pointing at the lens, **x-polarised**
(E along x). Model horn: cos^q pattern, q = 14.38, 10 dB edge taper at the lens edge. The beams point to
theta = +-30 deg in the x-z plane (phi = 0 and 180 deg) on the exit side (z < 0).

## Model predictions to compare with HFSS (all model-only: HFSS comparison NOT_AVAILABLE_FROM_CURRENT_DATA)
- Peak realised gain per beam (this model, horn power = 1, includes spillover, dissipation, model radiation efficiency): **25.6 dBi**
  at theta = 30.00 deg (phi 0) and 30.00 deg (phi 180); level difference **-0.00 dB**
  (zero by the layout's mirror symmetry, not an independent measurement).
- Co-pol side-lobe level relative to the beams: **-9.7 dB**; max cross-pol -19.3 dBi (off-plane).
- Aperture-averaged power transmission of the cells: 0.817 (-0.88 dB);
  fraction of horn power reaching the exit aperture: 0.752; radiated into propagating waves: 0.596.
- Per-cell phase residual vs requested: max 0.76 deg, rms 0.086 deg (cells above 1 deg: 0).
- Pillar|hole adjacencies in this layout: 432 of 1452 internal cell edges (29.8 %); median target-phase step across them 122 deg.
- **Local-periodicity caveat** (1x4 strips, order 0, normal incidence; `results/transmitarray_dual30/local_periodicity.json`): the model treats each cell as if it sat in an infinite array of itself. Against exact supercell RCWA the error vector, relative to the cell amplitude, is ~0.07 at a 60 deg neighbour step, ~0.2 at 110 deg and 0.4-0.5 at 150-175 deg, at pillar|hole and same-family boundaries alike; two dense-hole|thin-pillar pairs were worse (up to 0.70). Expect the HFSS pattern to depart from the numbers above (lower gain, higher side lobes) -- they are an optimistic, model-only bound.

## Suggested HFSS setup
Material nylon (eps_r 2.45, tan_d 0.0185); horn at the phase centre above; radiation (or PML) boundary >= lambda/4 (2.7 mm)
beyond the structure on all sides; far-field sphere, Ludwig-3 co/cross for x polarisation; fine mesh in the nylon (lambda/(10 n) ~ 0.7 mm).
Because the cell responses are local-periodicity predictions, expect the pillar|hole boundaries and the lens edge to be the first
places a full-wave run departs from the model. The lens has 756 cells; a symmetry plane (x and y) can be used.
