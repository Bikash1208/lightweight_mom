# Scientific feasibility gate (Phase 8A)

A read-only gate that checks, from the real Phase 7 production database only, whether the current
unit-cell library can physically implement a holographic-metasurface design *before* any GUI or
final geometry generation is attempted. It never modifies the solver, the mesh classifier, or the
Phase 7 results; it only reads them and reports.

```
results/unit_cell_library/ (Phase 7, unmodified)
        |
        v
feasibility_library_integrity.py --> feasibility_phase_source.py --> feasibility_coverage.py
        |                                    |                              |
        +------------------+-----------------+------------------------------+
                            v
                   feasibility_gate.py  --> results/phase8a_feasibility/
                     (also: feasibility_grating.py, feasibility_spatial.py,
                      reference_field.py, feasibility_required_phase.py,
                      feasibility_matching.py, feasibility_plots.py)
```

## The one fact that shapes this whole phase

This project has **no implemented holographic-synthesis formulation**
(`README.md`: "Holographic synthesis is not implemented"; `DATABASE.md`; `src/__init__.py`) and
**no reference-field / beam-steering / spatial-mapping module** of any kind (grep-verified: no
`aperture`/`beam`/`reference_field` module anywhere in `src/` before this phase). Consequently:

* `feasibility_phase_source.validate_holographic_phase_source()` cannot determine what quantity a
  holographic generator would require, so it always returns `status="PHASE_SOURCE_UNVALIDATED"`
  for this project as it stands — not a bug, the correct answer given the evidence.
* `feasibility_required_phase.compute_required_phase()` never computes a `required_phase(x, y)`
  array — brief instruction: "do NOT implement a second holographic formulation." It reports
  `status="HOLOGRAPHIC_FORMULATION_MISSING"` instead.
* `feasibility_spatial.check_mapping_method_defined()` always reports `defined=False`: no
  validated reference-grid → MTS-unit-cell mapping exists to check against.

Running the gate against the real library (`python scripts/research/run_phase8a_feasibility.py`) therefore
reports **`PHASE_SOURCE_UNVALIDATED`** (exit code 4) today. This is the expected, honest
conclusion (matches the brief's own Sec. 48), reached from real, computed data — not a
placeholder or a software defect.

## What is fully implemented and computed (Mode A, library-only)

* **Library integrity** (`feasibility_library_integrity.py`): 324 geometries, 201 frequency points
  per cell, one STEP file per successful cell, no duplicate `cell_id`, JSON/CSV consistency, and
  every manifest number cross-checked against the actual files (never trusted alone).
* **Phase coverage** (`feasibility_coverage.py`), reusing `phase_extraction.phase_coverage`
  exactly (Phase 7, unmodified) so the numbers are bit-identical/regression-tested against
  `LIBRARY_REPORT.md`'s own 24.84°/335.16°:
  * at the design frequency (12.5 GHz): span, largest gap, usable-geometry count, phase density
    (smallest/median/largest gap, distinct-phase-state count);
  * across all 201 stored frequency points (`phase_coverage_vs_frequency.csv`/`.png`): coverage
    ranges 23.24°-25.36° across the whole 10-15 GHz band — the narrow coverage is not specific to
    12.5 GHz;
  * resonance (`library_classifier.detect_resonance`, unmodified): 0 cells RESONANT in-band, so
    Sec. 16's coverage-vs-resonance correlation is honestly reported as "cannot be computed" from
    zero data points, not fabricated.
* **Grating/period check** (`feasibility_grating.py`): reuses
  `discretization.grating_onset_frequency`'s exact formula (cross-checked by test); onset ≈
  62.46 GHz, far above the 10-15 GHz band, `is_grating_risk=False`.
* **An incidental finding worth recording**: at the library's chosen 32×32 mesh (0.15 mm cells,
  the same mesh-coarseness trade-off `LIBRARY_REPORT.md` Sec. 16 already documents),
  `phase_density` finds only **169 distinct phase states among the 324 successful cells** at
  12.5 GHz — e.g. `(L=4.4, W=3.8)` and `(L=4.4, W=4.0)` mm produce a *bit-identical* stored S11xx
  value. This is consistent with, not contradictory to, the already-documented mesh-coarseness
  limitation (some distinct nominal geometries discretize to the same mesh pattern); it is not
  investigated further here (out of this phase's scope, and no solver code may be changed to chase
  it), but it further reduces the library's *effective* usable phase diversity below its nominal
  324-geometry count.

## What is implemented but structurally cannot pass yet (Mode B, design-specific)

Given no reference field exists to actually supply, Mode B's machinery is implemented and unit
tested with synthetic inputs, but is not exercised end-to-end against a real design:

* `reference_field.py`: reference-field CSV validation (format/coordinates/grid/NaN/duplicates)
  and beam-definition validation (theta/phi/weight/phase/polarization) — a minimal schema Phase 8A
  defines for its own input-validation purposes, explicitly **not** claimed as an existing project
  convention.
* `feasibility_matching.py`: circular phase matching and `is_phase_coverage_sufficient()` against
  *any* required-phase array (real or synthetic) — fully implemented, fully tested (partial-range
  feasible vs. full-range-insufficient synthetic cases, brief Sec. 44 Tests 10-11), never
  extrapolates and never rescues a mismatch with an arbitrary phase offset.

## Locked prohibitions actually enforced

No phase unwrapping/interpolation/extrapolation, no invented geometry, no arbitrary phase offset,
no silent nearest-neighbour spatial mapping, no claim of HFSS validation anywhere (every output —
`feasibility.json`, `feasibility_report.md`, the CLI banner — says `NOT PERFORMED`, tested
directly).

## Running it

```bash
.venv/bin/python scripts/research/run_phase8a_feasibility.py        # Mode A against results/unit_cell_library
```

Exit codes: `src/feasibility_schema.py::EXIT_CODES` (0 = `DESIGN_FEASIBLE`, 4 = the current
`PHASE_SOURCE_UNVALIDATED` result, ...). Output: `results/phase8a_feasibility/` (`feasibility.json`,
`feasibility_report.md`, phase-coverage CSVs, and the 3 Mode-A plots — the 3 design-specific plots
are correctly absent when no design-specific input was supplied).

## Verification

`tests/test_feasibility_schema.py`, `test_feasibility_phase_source.py`,
`test_feasibility_library_integrity.py`, `test_feasibility_coverage.py`, `test_reference_field.py`,
`test_feasibility_spatial.py`, `test_feasibility_matching.py`, `test_feasibility_grating.py`,
`test_feasibility_plots.py`, `test_feasibility_gate.py` — small deterministic fixtures plus
(skipped if absent) direct checks against the real 324×201 production library.
