# LIBRARY_REPORT.md

Production rectangular-PEC-patch unit-cell library (Phase 7), generated from the existing, unmodified lightweight periodic MoM solver (Phases 1-6A). No EM formulation, Green's function, MoM matrix, basis functions or S-parameter code was changed to build this library.

## 1. Physical configuration

- Period: 4.8 x 4.8 mm (square lattice).
- Substrate: `Rogers RT/Duroid 5880` (authoritative source: `config.yaml`, this project's own root config), eps_r = 2.2, tan_delta = 0.0009, thickness = 1.575 mm.
- Metal: PEC (`materials.py`, Phase 1's only supported model).
- Patch sweep: L, W in {1, ..., 4.4} mm, step 0.2 mm (18 x 18 = 324 geometries).
- Mesh: 32x32 cells, Floquet order 32 (accuracy trade-off documented in Sec. 16, "Limitations").

## 2. Number of geometries

324 enumerated (18x18); 324 solved OK, 0 invalid/skipped (a priori geometry check, before any solve). No geometry was silently dropped or substituted (req. 4/28).

## 3. Frequency range

10 GHz to 15 GHz, step 25 MHz, 201 points inclusive; design frequency 12.5 GHz is an exact grid point (verified, never interpolated).

## 4. Solver convention

Time convention: `exp(+j*omega*t)` (CONVENTIONS.md). S-parameters: 4x4 [[S11, S12], [S21, S22]] over (port1 x, port1 y, port2 x, port2 y); entry name S<out><in>_<received><incident>; port 1 = top reference plane (wave incident travelling -z), port 2 = bottom reference plane (wave incident travelling +z). Phase: arg(S), principal value in (-180, 180] deg; a delay is a negative phase. No phase unwrapping is performed anywhere in this library (req. 7, 9).

## 5. Polarization convention

global x and y unit vectors, identical at both ports; 4x4 [[S11, S12], [S21, S22]] over (port1 x, port1 y, port2 x, port2 y); entry name S<out><in>_<received><incident>; port 1 = top reference plane (wave incident travelling -z), port 2 = bottom reference plane (wave incident travelling +z)

## 6. Database structure

`database/unit_cell_library.json` (one entry per geometry, nested by frequency), `database/unit_cell_library.csv` (flat, one row per (cell, frequency), 24 columns), `database/phase_lookup_12p5GHz.csv` (req. 10's design-frequency lookup table). Both JSON and CSV re-export the same underlying `database.Database` (Phase 5, unmodified) this build produced; internal consistency between them: **consistent** (`reports/LIBRARY_VALIDATION.md`).

**HFSS comparison compatibility (req. 23):** this library uses exactly the conventions Phase 6A's `HFSS_COMPARISON.md` already documents per case, not a second incompatible one: `S11`/`S21`/`S12`/`S22` here are the same `S_ij[received, incident]`, index 0=x, 1=y, `exp(+j*omega*t)` quantities as `S_PARAMETERS.md`; `MoM S11xx <-> HFSS S(Port1:x, Port1:x)`, `MoM S21xx <-> HFSS S(Port2:x, Port1:x)` (Floquet port modes x/y = the (0,0)-order linear polarizations, `HFSS_COMPARISON.md` Sec. 2). If your HFSS export uses `exp(-j*omega*t)`, conjugate every HFSS S-parameter before comparing.

## 7. Geometry generation

324 STEP files written to `geometry/` (`UC_L<L>_W<W>.step`), one per successfully-solved geometry, built directly from the same `UnitCell` the solver used (`cad_export.py`, Phase 6, unmodified). Every export passed its independent bounds check (req. 17: 0 < L < period_x, 0 < W < period_y, and the patch bounding box strictly inside the cell). STEP cannot represent a periodic boundary condition itself; the unit-cell footprint (substrate/metal solids) stands in for it, per `cad_export.py`.

## 8. Phase extraction method

Component: `S11_xx` (status: `SPECIFICATION_MANDATED`). Req. 10's lookup-table schema names S11 explicitly; this project's existing solver has no ground plane, so S21 is not negligible here (unlike a typical grounded reflectarray element) and whether S11 or S21 phase is the one the eventual holographic generator needs is an architecture decision this project has not made. Both S11 and S21 (all four blocks) are stored in full in the library database regardless of this choice; only the phase-vs-geometry table and the coverage/collision/resonance analyses below are pinned to S11_xx.

## 9. Phase coverage

324 successful cells at 12.5 GHz: phase range [60.08, 84.91] deg (naive, non-circular span 24.84 deg); circularly-aware coverage = **24.84 deg of 360 deg**. See `plots/phase_coverage_map.png` and `plots/phase_vs_geometry.png`.

## 10. Largest phase gap

**335.16 deg**, between 84.91 and 60.08 deg (wrapping through +-180 deg if the gap straddles it). This is the single largest arc of phase this library cannot represent with any geometry in the sweep; 3 phase region(s) (`library_schema.PHASE_COLLISION_BIN_DEG` = 10 deg bins) have more than one candidate geometry (no winner chosen automatically, req. 11).

## 11. Resonance analysis

0 cell(s) classified RESONANT, 324 NON_RESONANT, out of 324 solved cells (criterion: local |S11xx(f)| maximum at a strictly interior frequency sample with prominence >= 0.1, `library_classifier.detect_resonance`). See `plots/resonance_map.png`.

## 12. Classification rules

- Phase usability: PHASE_USABLE if |S11xx| at design frequency >= -40 dB (else PHASE_UNUSABLE); 324 usable, 0 unusable (`library_schema.PHASE_USABLE_MIN_MAG_DB`).
- Grating risk: GRATING_RISK if a non-specular Floquet order propagates at any frequency in the sweep (`discretization.is_specular_only`, unmodified); onset frequency for this period = 62.46 GHz. 0 of 324 cells flagged.
- Every threshold above is a named constant in `library_schema.py`, documented there and exercised at its boundary in `tests/test_library_classifier.py`.

## 13. Failed cells

None.

## 14. STEP export status

324 of 324 successfully-solved geometries exported as STEP; every export's substrate/metal solids independently re-validated as topologically valid solids on re-import (`cad_export.solid_validity`, `BRepCheck_Analyzer`) in `tests/test_geometry_export.py`.

## 15. HFSS validation subset

See `reports/HFSS_VALIDATION_SUBSET.csv` and `reports/LIBRARY_VALIDATION.md` -- a deterministic recommendation only. **No HFSS simulation has been run or imported for this library; no HFSS agreement is claimed anywhere in this report** (req. 22, 31.18).

## 16. Limitations

- No ground plane in the underlying solver (README "Assumptions"); S21 is not negligible for these cells, unlike a typical grounded reflectarray element -- see Sec. 8's phase-source discussion.
- Mesh (32x32 cells, Floquet order 32) is below this project's own "2-4x cells" practical accuracy guideline (`MOM_FORMULATION.md`) purely for build-time cost; every cell's a priori `metal_area_error`/`area_ok` (Phase 5, unmodified) is stored in `database/unit_cell_library.json`'s source `database.json` (and reachable via `database.Database.cell(...).validation`) rather than hidden.
- STEP export cannot represent a periodic lattice-pair boundary condition; only the physical solids are exported (Sec. 7).
- No interpolation or extrapolation anywhere in this library or its query API (req. 26); `find_phase_candidates`/`invert_phase` only ever return geometries actually simulated.

## 17. Reproducibility

`manifest.json`'s `configuration_hash` (fa2708bd1eee4c80...) identifies every physical input to this build. `library_builder.py` reuses `runner.run_sweep`'s checkpointed, resumable, deterministic solve (Phase 5, unmodified); `tests/test_library_builder.py` verifies numerically identical output across two builds of the same small configuration. `creation_timestamp` is the one field that is not reproducible run to run, by design (same convention as every earlier phase's file headers).

## 18. HFSS validation status

**HFSS validation is pending.** This report and this library are MoM-only results. No claim of agreement with HFSS (or any other full-wave tool) is made anywhere in this project until actual HFSS results are imported and quantitatively compared -- see `reports/LIBRARY_VALIDATION.md` and `VALIDATION_EXPORT.md`.
