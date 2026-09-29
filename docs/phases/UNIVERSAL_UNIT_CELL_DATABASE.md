# UNIVERSAL_UNIT_CELL_DATABASE.md

Phase 8D-REV: a period-agnostic, topology-aware unit-cell electromagnetic database architecture, built as a query/registry layer over the existing, unmodified `database.Database` objects Phases 1-7 already produce. Nothing in the MoM solver, mesh classifier, Green's function, or Phase 7 database/geometry was modified or re-run.

## 1. Why 4.8 mm is not universal

`library_schema.LibraryConfig.default()`'s `period_m = 4.8e-3` (`_PERIOD_M`) is a locked *Phase 7 build choice*, not a solver constant: the underlying `periodic_bc.RectangularLattice(period_x, period_y)` (Phase 1, unmodified) already accepts any positive period -- it was never hard-coded at the solver level. Phase 8D-REV's job is to stop that one build's choice from leaking into the design layer as an assumption.

## 2. Database hierarchy

`UniversalDatabase` (a registry) holds any number of `LibraryEntry` objects, each one an independent solved `database.Database` tagged with its own topology/period/material/frequency-range/polarization/incidence identity and provenance. `family.py` aggregates every entry sharing one topology into a `UnitCellFamily` (Sec. 6, 9); the existing 4.8 mm build is registered as exactly one entry, role `reference_validation_library`, library_id `phase7_reference_validation_library_4p8mm` -- never treated as *the* database.

## 3. Unit-cell family concept

`UnitCellFamily` for topology `rectangular_patch` reports `parameter_schema=('length_m', 'width_m')` (the geometry parameters that actually vary, computed from the data, not hand-listed) and `fixed_parameters={center_x_m, center_y_m, rotation_deg}`. A future family (slot-loaded patch, arm-length dipole, ...) needs no schema change: `record_schema.py`'s `extra_geometry_parameters` dict already carries any `database.PARAMETER_NAMES` beyond length/width generically.

## 4. Period parameterization

`LibraryEntry.p_x_mm`/`p_y_mm` are per-entry fields; `UniversalDatabase.periods_for_topology(topology)` currently returns `((4.8, 4.8),)` -- one validated period, from one registered build. Registering a second build at a different period adds a second tuple; nothing here assumes there is only one.

## 5. Complex scattering storage

`record_schema.UnitCellRecord.s_parameters` stores all 16 linear-basis entries (`database.component_names('linear')`) as `{re, im, mag, phase_deg}` -- the complex value is primary (Sec. 8); magnitude/phase are derived from it, never the reverse. Cross-polarized terms (`S11_xy`, `S11_yx`, `S21_xy`, `S21_yx`, ...) are read from the internal Phase 7 checkpoint (`phase8c.data_source.load_full_database`, itself read-only and self-cross-checking against the delivered library), not discarded the way the delivered `unit_cell_library.json`'s `_xx`-only export does.

## 6. Material parameterization

`LibraryEntry.substrate` (`{'name': 'Rogers RT/Duroid 5880', 'eps_r': 2.2, 'tan_delta': 0.0009, 'thickness_mm': 1.575}`) is a per-entry dict (name, eps_r, tan_delta, thickness_mm); `family.material_stacks` collects every distinct stack observed across all registered entries of one topology -- never a single hard-coded substrate.

## 7. Frequency parameterization

`LibraryEntry.frequency_min_hz`/`frequency_max_hz`/`n_frequency_points` (10-15 GHz, 201 points for this entry) are per-entry; `query.UnitCellQuery.frequency_hz` matches only an exact stored grid point of whichever entry is being searched -- never interpolated.

## 8. Polarization parameterization

`LibraryEntry.polarizations_supported = ('x', 'y')` (both linear incident polarizations, per `database.CONVENTIONS`'s linear basis) is a per-entry field, queried via `UnitCellQuery.polarization`.

## 9. Incidence parameterization

`LibraryEntry.incidence_supported = ('normal',)` -- today every registered entry is normal-incidence only (`phase8b.conventions.NORMAL_INCIDENCE_ONLY`, solver-wide, unmodified); the field exists specifically so a future oblique-incidence build is representable without a schema change.

## 10. Query mechanism

`query.search_unit_cells` applies topology -> material -> frequency -> period -> polarization -> incidence -> response-availability -> validity-status filters, in that fixed order (Sec. 14), and never ranks by phase coverage (see `query_spec.md`).

## 11. Project-library extraction

`project_library.extract_project_library(universal_db, query)` returns a `ProjectLibrary` referencing (never copying) the `LibraryEntry` objects a `MATCH`ed query selected -- the small, task-specific subset a real HMSA design would actually use (Sec. 17).

## 12. Missing-data behavior

A query that reaches zero candidates at any filter stage returns `DATABASE_COVERAGE_MISSING` with the exact stage named (`eliminated_at`) and the full per-stage candidate counts (`filter_trail`); `build_request.generate_build_request` turns that into a machine-readable `LibraryBuildRequest` -- nothing is rescaled, interpolated, or silently generated (Sec. 18, 21).

## 13. Library-building workflow

Unchanged: `library_schema.LibraryConfig` + `library_builder.build_library_database` + `unit_cell_library.build_full_library` (Phase 7, unmodified) build one `database.Database`; `library_entry.from_database` then registers it into a `UniversalDatabase`. The existing 4.8 mm build is simply the first (and, at this point, only) library built this way.

## 14. Validation levels

`validation_levels.VALIDATION_LEVELS` = DATA_PRESENT < NUMERICALLY_VALIDATED < CROSS_SOLVER_VALIDATED < HFSS_VALIDATED < EXPERIMENTALLY_VALIDATED. The benchmark entry's actual, measured status is `DATA_PRESENT` (every stored record's own `checks_passed` diagnostic -- Phase 5, unmodified -- determines this; nothing stronger is claimed: no cross-solver, HFSS, or experimental comparison exists anywhere in this repository).

## 15. Provenance

Every `UnitCellRecord.provenance` carries `source_type` (`"MoM"`), `solver_version` (`0.5.0`, from `manifest.json`), `configuration_hash` (`fa2708bd1eee4c80809cc8c7569b7acf445441e3e7301f9ec7f65b71eeac2f18`), and its own `validation_status` -- never presented as validated merely because it exists in the database (Sec. 22).

## 16. Integrity rules

`integrity.run_integrity_checks` (Sec. 24-25, tests UDB-01..UDB-12): unique record identity, valid parameter ranges, complete frequency grids, finite complex-valued responses, preserved cross-polarization, no fabricated values for failed cells, complete metadata, a content hash, unique library ids, and JSON round-trip exactness.

