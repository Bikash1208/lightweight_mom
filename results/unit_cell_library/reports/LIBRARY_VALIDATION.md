# LIBRARY_VALIDATION.md

## HFSS validation status

**No HFSS comparison has been performed.** Every number in this library is produced by this project's own lightweight periodic MoM solver (Phases 1-6A, unmodified). The subset below is a deterministic, documented *recommendation* of which cells to validate in HFSS first; it is not evidence that any of them agree with HFSS, because no HFSS simulation has been run or imported for this library (req. 22).

## Recommended HFSS validation subset

9 of 324 cells, selected deterministically (no random sampling, no manual picks) to cover: a small patch, a large patch, an approximately square medium patch, both extreme aspect ratios, the most resonant cell found, one non-resonant cell, and phase-quartile representatives at the design frequency. Full table: `HFSS_VALIDATION_SUBSET.csv`.

| cell_id | L_mm | W_mm | reason(s) |
|---|---|---|---|
| b92a9e8aeb43e660 | 1 | 1 | small patch (smallest L and W in the sweep); non-resonant cell (deterministic: smallest (L, W) with no detected resonance) |
| 1e445dd7687d4aba | 1 | 4.4 | small/large aspect ratio (narrowest L, widest W) |
| d595367f60c5b49e | 1.2 | 2.2 | phase-region representative (S11xx phase quartile sample, phase=83.6 deg) |
| cfa191715f326fc4 | 2.2 | 3 | phase-region representative (S11xx phase quartile sample, phase=84.9 deg) |
| 14a6d57a7536658e | 2.8 | 2.8 | approximately square, medium-sized patch (sweep midpoint) |
| 60fb0051afdf5074 | 3 | 1.4 | phase-region representative (S11xx phase quartile sample, phase=84.6 deg) |
| f48f9e9d071fcad8 | 3.2 | 3.6 | phase-region representative (S11xx phase quartile sample, phase=79.9 deg) |
| d4f002bd5d2b97cc | 4.4 | 1 | small/large aspect ratio (widest L, narrowest W) |
| 755ddf3509b3cc79 | 4.4 | 4.4 | large patch (largest L and W in the sweep); phase-region representative (S11xx phase quartile sample, phase=60.1 deg) |

## Internal consistency (JSON vs. CSV)

`unit_cell_library.json` (65124 records) vs. `unit_cell_library.csv` (65124 rows): **consistent**.

## Build completeness

- geometries: 324 expected, 324 solved OK, 0 invalid/skipped.
- computed (this run): 0; reused from checkpoint: 65124.
