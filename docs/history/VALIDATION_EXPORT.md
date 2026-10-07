# Validation unit-cell export (Phase 6)

Exports physical, importable CAD models of representative unit cells **for manual comparison
against a full-wave solver (HFSS, CST, openEMS, ...)**, built directly from the same geometry
parameters the MoM solver consumes. This phase touches no EM formulation: `green_function.py`,
`mom_matrix.py`, `solver.py`, `sparameters.py` and every earlier phase are unmodified and are not
imported by `cad_export.py` in the write direction — `cad_export.py` reads `UnitCell`/`Geometry`/
`Substrate` objects (Phase 1) but never writes back to or influences the solve.

```
UnitCell (Phase 1, unmodified) ──┬──> cad_export.py       ──> unit_cell.step
                                  └──> runner.run_sweep     ──> Database (Phase 5, unmodified)
                                         (Phase 5, unmodified)     │
                                                                   ├──> mom_results.json / .csv
                                                                   ├──> plots/*.png
                                                                   └──> README.md
```

## 1. What gets exported, per case

One folder per case (`<case_id>/`, `src/validation_export.py::export_case`):

| file | contents | built from |
|---|---|---|
| `unit_cell.step` | `Substrate` + `Metal` solids, one STEP (ISO-10303-21) assembly | `cad_export.export_unit_cell_step` |
| `geometry.json` | every geometry/material parameter given to the solver, plus the CAD-only metal export thickness | `case.config.unit_cell` |
| `solver_config.json` | mesh, Floquet order, reference planes, frequency grid, excitation convention, `database.CONVENTIONS` | `case.config`, `Database.metadata` |
| `mom_results.json` / `.csv` | the Phase 5 database format (DATABASE.md Sec. 6), for this one cell across its own frequency grid | `Database.to_json` / `to_csv`, unmodified |
| `plots/s11.png`, `plots/s21.png` | magnitude (linear) and phase (deg) vs frequency, `_xx` and `_yy` | `database.series` / `phase_series_deg` |
| `plots/co_polarized.png` | `S11_xx, S21_xx, S11_yy, S21_yy` magnitude, dB, one axes | ″ |
| `plots/cross_polarized.png` | `S11_yx, S21_yx, S11_xy, S21_xy` magnitude, dB, one axes (present even for unrotated cases, where it is ≈ −∞ dB / the noise floor) | ″ |
| `README.md` | every fact needed to reproduce the model manually in HFSS (Sec. 4 below), rendered from the case + its own results | `render_readme` |

## 2. Coordinate / unit mapping (the one place besides `units.py` that converts)

The solver works internally in SI (metres, radians, hertz); STEP's native unit is **millimetres**.
`cad_export.py` converts with `units.to_mm`/`to_deg` in both build directions (export, and the
independent closed-form checks used by the tests) — no parameter is re-entered by hand.

* `x, y`: as in `geometry.py` — origin at the cell centre, lattice-aligned axes, unchanged.
* `z`: the solver's stack (`GREEN_FUNCTION.md`) is air (`z > 0`) | substrate (`−h < z < 0`) |
  zero-thickness PEC pattern at `z = 0` | air (`z < −h`). The exported **`Substrate`** solid is a
  `period_x × period_y × h` box occupying `z ∈ [−h, 0]`. The solver's metal is an idealized
  zero-thickness surface current; a CAD solid cannot have zero volume, so the exported
  **`Metal`** solid is a thin slab `z ∈ [0, t]`, flush with the substrate's top face, where `t`
  is `metal_export_thickness_m` (default `cad_export.DEFAULT_METAL_EXPORT_THICKNESS_M` = 35 µm,
  a nominal 1 oz/ft² copper thickness).
* **`t` has zero effect on any solved quantity.** It is not read anywhere outside `cad_export.py`
  and `validation_export.py`; it exists only so the STEP solid is a solid. In HFSS, assign the
  `Metal` object a **Perfect E** boundary (not a lossy copper volume) so the imported model
  matches the zero-thickness assumption the S-parameters were actually computed under — see each
  case's own `README.md` Sec. "Reproducing this model in HFSS" for the concrete steps.
* Rotation: `geometry.py`'s convention is counter-clockwise about the feature's own centroid,
  then translated to `(center_x, center_y)`. `cad_export._feature_solid` reproduces this exactly
  (`cq.Workplane.rotate` about `+Z` matches `Geometry.corners()`'s sign convention — checked by
  `tests/test_cad_export.py::test_rotation_changes_recovered_bbox` and the closed-form bbox
  cross-check below, not merely asserted).
* `polarity == "metal"` (`rectangle`, `patch`, `dipole`): the `Metal` solid *is* the feature.
* `polarity == "aperture"` (`slot`): the `Metal` solid is a full-period plate with the feature
  boolean-cut out (`cq.Workplane.cut`), matching `UnitCell.is_metal`'s complement rule exactly —
  the same rule `metal_area` and the solver's own `is_metal` use, not a separate one.

## 3. Case selection

Ten cases (`src/validation_cases.py::CASES`), chosen to cover every category the Phase 6 brief
asks for: non-resonant, resonant, narrow, wide, four distinct (length, width) pairs, two rotated
(one of them also off-centre), two slot/complement (one of them also rotated). All ten share one
lattice, substrate and reference-plane configuration (so only geometry, mesh and frequency grid
vary and the substrate solid is identical/comparable across cases); case03 (the narrow, 4:1
aspect-ratio dipole) uses a finer 32×32 mesh so the rooftop mesh resolves the narrow dimension at
all — a coarser mesh made that case's metal footprint fall exactly between mesh-cell boundaries
under floating-point rounding, yielding zero metal cells (see `validation_cases.py`'s module
docstring for the full account, including how the resonant case's centre frequency was found by
*measuring* the solver's own output, not assumed).

The full case table and the specific numbers behind "resonant"/"non-resonant"/etc. are documented
in `src/validation_cases.py`'s module docstring — this file does not repeat them, to avoid the two
falling out of sync; read that docstring for the case catalog itself.

## 4. What the README documents, per case

Every `README.md` states explicitly (checked field-by-field by `tests/test_validation_export.py`,
not just present by inspection): the frequency range and step; the unit-cell period; substrate
`eps_r`/`tan_delta`/thickness; the metal geometry (`type`, `length`, `width`, `center`, `rotation`,
polarity, fill fraction); the mesh (`cells_x`, `cells_y`, cell size) and number of rooftop
unknowns; the Floquet truncation order and the rank rule it satisfies; the excitation convention
(unit-amplitude normal incidence, both x and y solved regardless of any single `polarization`
value in the case config); the reference-plane locations (`z_top`, `z_bottom`) and that only the
specular (0,0) order is a port mode; the literal **`exp(+jωt)`** time convention and what to do if
comparing against an `exp(−jωt)` tool (conjugate the S-parameters); and an explicit 8-step manual
HFSS reproduction recipe (units, material, boundary, periodic box and boundaries, Floquet ports,
frequency sweep, and how to compare) ending in an explicit statement that this project performs
**no automatic comparison** against HFSS and applies **no tuning** of the MoM results to match it.

## 5. Determinism / reproducibility — and one honest caveat

* **STEP file bytes are not reproducible across runs**, even for the identical geometry: OpenCASCADE's STEP writer embeds a
  generation timestamp in the header (`FILE_NAME(...,'<timestamp>',...)`). The **geometric content**
  is reproducible: re-importing two independently exported files of the same case and comparing
  bounding boxes/volumes (`cad_export.import_and_summarize`) agrees to floating-point precision
  (`tests/test_validation_export.py::test_step_geometric_content_reproducible_despite_timestamp`).
* `mom_results.json`/`.csv` reuse Phase 5's `Database.to_json`/`to_csv`, which DATABASE.md Sec. 5
  documents as byte-identical for repeated runs *given* an identical `Database` object. In
  practice, two **independent solves** of the same case are not always bit-identical here: this
  project's `numpy` is linked against Apple's Accelerate BLAS/LAPACK (`numpy.show_config()`),
  which is multi-threaded and, for some geometries (case01 among them), does not reduce
  floating-point sums in the same order on every call. The effect is confined to the ~1e-9
  relative level in derived, already-near-zero diagnostics (`reciprocity_error`,
  `passivity_min_eig`); the S-matrix itself (what the plots and README show) is far more tightly
  reproducible because it comes directly from the linear solve rather than a subsequent
  eigendecomposition. This is a **pre-existing, environment-level property of the unmodified
  solver stack**, reproduced with plain `runner.run_sweep` calls with no Phase 6 code involved —
  not introduced by, or in scope to fix from, this phase. `tests/test_validation_export.py`
  checks numerical closeness across repeated exports rather than asserting false byte identity.
* CAD exports do not use randomness and do not depend on wall-clock time except for that one
  STEP header field.

## 6. Verification (`tests/test_cad_export.py`, `test_validation_cases.py`, `test_validation_export.py`, `test_hfss_comparison.py`)

| property | how it is checked |
|---|---|
| every exported STEP opens as a valid solid model | OpenCASCADE's own `BRepCheck_Analyzer`, on both the pre-export solid and an independent re-import from disk, for every ad hoc test cell and all 10 shipped cases |
| STEP parameters == MoM geometry parameters | substrate bbox == `(±Px/2, ±Py/2, −h, 0)` (from `units.to_mm`, not cadquery); metal bbox == closed-form rotated-rectangle bbox (`expected_feature_bbox_mm`, independent of cadquery); metal volume / thickness == `UnitCell.metal_area` (mm²); metal centroid == `(center_x, center_y)` for metal polarity, cell centre for aperture polarity — all checked against values recovered by *re-importing the STEP file from disk*, not just the pre-export in-memory solid |
| non-vacuousness | rotation, centre offset and export thickness are each perturbed independently and shown to change the corresponding recovered quantity (bbox, centroid, volume respectively) — the consistency checks are sensitive to the geometry, not passing regardless of input |
| case catalog | all 10 cases build valid configs; every requested category (non-resonant, resonant, narrow, wide, 4 distinct length/width pairs, rotated, slot) is present; rotated cases have nonzero rotation; slot cases are aperture polarity; case construction is deterministic |
| folder layout | exactly the 7 required entries (`unit_cell.step, geometry.json, solver_config.json, mom_results.csv, mom_results.json, plots/, README.md`) and exactly the 4 required plot files |
| `geometry.json` / `solver_config.json` content | every field cross-checked against `case.config`/`case.config.unit_cell` directly, for a plain, a rotated and a slot case |
| `mom_results.json`/`.csv` | independently re-solved with a fresh, un-orchestrated `runner.run_sweep` call (same pattern as the Phase 5 database tests) and compared to the exported file, not merely internally consistent with itself |
| plots | PNG signature and nonzero size; the underlying `_mag_phase`/`_db` helpers checked against `Database.series`/`phase_series_deg` ground truth directly |
| README content | every required fact (Sec. 4 above) checked by regex/substring, including the literal `exp(+j*omega*t)` string and the "does not perform this comparison automatically" / "does not tune" disclaimers |
| error handling | `metal_export_thickness_m <= 0` rejected; `cadquery` absent raises `ImportError` naming `requirements-cad.txt` (`pytest.importorskip("cadquery")` at module scope, so the base suite is unaffected either way) |
| Phase 6A CSV/plots (`test_hfss_comparison.py`) | exact requested column names/order; values cross-checked against an independent `Database.from_json` read (not the module's own intermediate table); unrotated case's cross-pol columns sit at the documented dB floor; `README.md`'s facts checked by substring; missing `mom_results.json` raises `FileNotFoundError` rather than solving; re-export from the same input file is numerically identical |

## 6a. Exact overlay-ready comparison export (Phase 6A, `src/hfss_comparison.py`)

For a single case (`case01_non_resonant_patch` by default), `scripts/research/export_hfss_comparison.py`
adds, into that case's already-existing folder:

* `hfss_comparison_reference.csv` — `frequency_GHz`, then `S11xx/S11yy/S21xx/S21yy_mag_dB` and
  `_phase_deg`, then `S11xy/S11yx/S21xy/S21yx_mag_dB` — the exact columns and order needed to
  overlay an HFSS export column-for-column.
* `plots/hfss_{s11,s21}_{magnitude_db,phase_deg}.png` and `plots/hfss_cross_polarization_db.png`
  — one curve pair (or quadruple, for cross-pol) per plot, `20*log10(|S|)` / principal phase in
  degrees, at exactly the case's stored frequency points.
* `HFSS_COMPARISON.md` — the MoM↔HFSS quantity/port/mode mapping, this project's `exp(+jωt)`
  convention and how to handle an `exp(−jωt)` HFSS export, the exact reference-plane
  z-coordinates, and an explicit statement that no unwrapping, normalization, smoothing,
  interpolation or empirical correction is applied anywhere in this export.

**This module reads only the already-exported `mom_results.json`/`geometry.json`/
`solver_config.json` of an existing case folder (`load_case_artifacts`) — it never calls
`runner.run_sweep` or otherwise re-solves.** Run `export_validation_cases.py` (or
`validation_export.export_case`) for that case first; `hfss_comparison.export_hfss_comparison`
raises `FileNotFoundError` rather than solving if `mom_results.json` is missing.

## 7. What this phase does **not** do

* Does not compare exported models against HFSS (or any other tool) automatically — every
  README's reproduction recipe is manual, by design (Phase 6 brief item 7).
* Does not tune, correct or otherwise adjust the MoM results toward any expected HFSS behaviour
  (item 8); the numbers in `mom_results.json`/`.csv` are exactly what `runner.run_sweep` (Phase 5,
  unmodified) produces for that case's config.
* Does not model finite metal thickness or conductivity; `metal_export_thickness_m` is a CAD
  solid-rendering convenience only (Sec. 2).
* Does not add oblique-incidence ports, grating-order ports, or any new EM capability — the
  frequency/geometry ranges of the 10 cases were chosen (by measuring the existing solver's own
  output) to stay specular-only and within the existing model's assumptions throughout.
