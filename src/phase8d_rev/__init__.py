"""Phase 8D-REV: universal, period-agnostic, topology-aware unit-cell database architecture.

This package does not solve anything: it is a query/registry/schema layer over
`database.Database` objects that Phases 1-7's existing, unmodified solver already produces (one
per "library build"). Its job is to stop the project from treating the Phase 7 build's own
choices -- period = 4.8x4.8 mm, one substrate, one topology, 10-15 GHz -- as global constants.

Central objects (see the individual modules for detail):

- `record_schema.UnitCellRecord` -- one unit-cell/frequency response, the schema of Sec. 7.
- `family.UnitCellFamily` -- what a topology's parameter/period/frequency/polarization/incidence
  coverage actually is, derived from the data that exists, never assumed.
- `library_entry.LibraryEntry` -- one registered `database.Database` plus its identity/provenance.
- `universal_database.UniversalDatabase` -- a registry of `LibraryEntry` objects; `benchmark()`
  wraps the existing, byte-unmodified Phase 7 324-cell/4.8 mm library as exactly one entry among
  possibly many, tagged `reference_validation_library` (Sec. 3).
- `query.search_unit_cells` -- Sec. 14's physical-filter-first query engine; returns either
  matches or an explicit `DATABASE_COVERAGE_MISSING` (never fabricated/interpolated data, Sec. 21).
- `build_request.LibraryBuildRequest` -- the machine-readable "build this" specification a missed
  query produces (Sec. 18-19); nothing is silently generated from it.
- `response_matching.nearest_response_match` -- Sec. 15's complex-response nearest-neighbour
  search, reporting amplitude/phase/complex error separately, never collapsed into one score.
- `candidate_period.py` -- Sec. 12-13: CANDIDATE vs VALIDATED period, never a universal
  0.2*lambda-style threshold invented here.
- `integrity.py` -- Sec. 24-25's UDB-01..UDB-12 checks.

Frozen (never modified, never re-run): the MoM solver, mesh classifier, Green's function, Phase 7
database/geometry/STEP files, `results/unit_cell_library/` in its entirety. The existing 4.8 mm
library is read wholesale, unmodified, as one `LibraryEntry` -- it is a `benchmark_library`
proving the architecture works, not the final or only database (Sec. 3, Sec. 20).

No GUI, no final MTS STEP, no aperture calculation, and no bypass of the still-open Phase 8C
`CONTROL_VARIABLE_UNVALIDATED` result (Sec. 28-29): this package is agnostic about which
scattering quantity the eventual holographic formulation will use.
"""
