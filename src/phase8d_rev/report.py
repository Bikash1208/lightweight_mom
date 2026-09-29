"""Sec. 27, 32: `UNIVERSAL_UNIT_CELL_DATABASE.md` (the 16-topic architecture explainer) and
`validation_report.md` (Sec. 31's acceptance checklist + Sec. 32's mandatory final answer).
"""

# Sec. 31: every acceptance-criteria line, in order, each mapped to the boolean this run actually
# measured -- never asserted without a corresponding check.
ACCEPTANCE_CRITERIA = (
    "4.8 mm is no longer a global design assumption",
    "Existing 4.8 mm Phase-7 library remains untouched",
    "Period is represented as a database parameter",
    "Multiple periods are supported by the schema",
    "Multiple topologies are supported",
    "Multiple substrates/material stacks are supported",
    "Frequency is parameterized",
    "Polarization is parameterized",
    "Incidence conditions are parameterized",
    "Complete available complex S-matrix data are preserved",
    "Cross-polarization data are preserved",
    "Database provenance is preserved",
    "Missing combinations return an explicit failure/status",
    "No response is interpolated or extrapolated between unsupported periods",
    "No geometry is rescaled to fake another period",
    "Query engine can recover the existing 324-cell 4.8 mm benchmark",
    "Query engine does not assume 4.8 mm",
    "Project-specific libraries can be extracted from the universal database",
    "Missing project libraries produce a precise build request",
    "Phase-8D control-variable uncertainty remains explicit",
    "Full regression suite passes",
)


def render_universal_unit_cell_database_md(family, entry):
    return "\n".join([
        "# UNIVERSAL_UNIT_CELL_DATABASE.md", "",
        "Phase 8D-REV: a period-agnostic, topology-aware unit-cell electromagnetic database "
        "architecture, built as a query/registry layer over the existing, unmodified "
        "`database.Database` objects Phases 1-7 already produce. Nothing in the MoM solver, mesh "
        "classifier, Green's function, or Phase 7 database/geometry was modified or re-run.", "",
        "## 1. Why 4.8 mm is not universal", "",
        "`library_schema.LibraryConfig.default()`'s `period_m = 4.8e-3` (`_PERIOD_M`) is a locked "
        "*Phase 7 build choice*, not a solver constant: the underlying `periodic_bc."
        "RectangularLattice(period_x, period_y)` (Phase 1, unmodified) already accepts any "
        "positive period -- it was never hard-coded at the solver level. Phase 8D-REV's job is to "
        "stop that one build's choice from leaking into the design layer as an assumption.", "",
        "## 2. Database hierarchy", "",
        "`UniversalDatabase` (a registry) holds any number of `LibraryEntry` objects, each one an "
        "independent solved `database.Database` tagged with its own topology/period/material/"
        "frequency-range/polarization/incidence identity and provenance. `family.py` aggregates "
        "every entry sharing one topology into a `UnitCellFamily` (Sec. 6, 9); the existing "
        f"4.8 mm build is registered as exactly one entry, role `{entry.role}`, library_id "
        f"`{entry.library_id}` -- never treated as *the* database.", "",
        "## 3. Unit-cell family concept", "",
        f"`UnitCellFamily` for topology `{family.topology}` reports `parameter_schema="
        f"{family.parameter_schema}` (the geometry parameters that actually vary, computed from "
        "the data, not hand-listed) and `fixed_parameters={" + ", ".join(family.fixed_parameters) +
        "}`. A future family (slot-loaded patch, arm-length dipole, ...) needs no schema change: "
        "`record_schema.py`'s `extra_geometry_parameters` dict already carries any "
        "`database.PARAMETER_NAMES` beyond length/width generically.", "",
        "## 4. Period parameterization", "",
        f"`LibraryEntry.p_x_mm`/`p_y_mm` are per-entry fields; `UniversalDatabase."
        f"periods_for_topology(topology)` currently returns `{family.supported_periods_mm}` -- "
        "one validated period, from one registered build. Registering a second build at a "
        "different period adds a second tuple; nothing here assumes there is only one.", "",
        "## 5. Complex scattering storage", "",
        "`record_schema.UnitCellRecord.s_parameters` stores all 16 linear-basis entries "
        "(`database.component_names('linear')`) as `{re, im, mag, phase_deg}` -- the complex "
        "value is primary (Sec. 8); magnitude/phase are derived from it, never the reverse. Cross-"
        "polarized terms (`S11_xy`, `S11_yx`, `S21_xy`, `S21_yx`, ...) are read from the internal "
        "Phase 7 checkpoint (`phase8c.data_source.load_full_database`, itself read-only and self-"
        "cross-checking against the delivered library), not discarded the way the delivered "
        "`unit_cell_library.json`'s `_xx`-only export does.", "",
        "## 6. Material parameterization", "",
        f"`LibraryEntry.substrate` (`{entry.substrate}`) is a per-entry dict (name, eps_r, "
        "tan_delta, thickness_mm); `family.material_stacks` collects every distinct stack observed "
        "across all registered entries of one topology -- never a single hard-coded substrate.", "",
        "## 7. Frequency parameterization", "",
        f"`LibraryEntry.frequency_min_hz`/`frequency_max_hz`/`n_frequency_points` "
        f"({entry.frequency_min_hz/1e9:g}-{entry.frequency_max_hz/1e9:g} GHz, "
        f"{entry.n_frequency_points} points for this entry) are per-entry; `query.UnitCellQuery."
        "frequency_hz` matches only an exact stored grid point of whichever entry is being "
        "searched -- never interpolated.", "",
        "## 8. Polarization parameterization", "",
        f"`LibraryEntry.polarizations_supported = {entry.polarizations_supported}` (both linear "
        "incident polarizations, per `database.CONVENTIONS`'s linear basis) is a per-entry field, "
        "queried via `UnitCellQuery.polarization`.", "",
        "## 9. Incidence parameterization", "",
        f"`LibraryEntry.incidence_supported = {entry.incidence_supported}` -- today every "
        "registered entry is normal-incidence only (`phase8b.conventions.NORMAL_INCIDENCE_ONLY`, "
        "solver-wide, unmodified); the field exists specifically so a future oblique-incidence "
        "build is representable without a schema change.", "",
        "## 10. Query mechanism", "",
        "`query.search_unit_cells` applies topology -> material -> frequency -> period -> "
        "polarization -> incidence -> response-availability -> validity-status filters, in that "
        "fixed order (Sec. 14), and never ranks by phase coverage (see `query_spec.md`).", "",
        "## 11. Project-library extraction", "",
        "`project_library.extract_project_library(universal_db, query)` returns a `ProjectLibrary` "
        "referencing (never copying) the `LibraryEntry` objects a `MATCH`ed query selected -- the "
        "small, task-specific subset a real HMSA design would actually use (Sec. 17).", "",
        "## 12. Missing-data behavior", "",
        "A query that reaches zero candidates at any filter stage returns "
        "`DATABASE_COVERAGE_MISSING` with the exact stage named (`eliminated_at`) and the full "
        "per-stage candidate counts (`filter_trail`); `build_request.generate_build_request` turns "
        "that into a machine-readable `LibraryBuildRequest` -- nothing is rescaled, interpolated, "
        "or silently generated (Sec. 18, 21).", "",
        "## 13. Library-building workflow", "",
        "Unchanged: `library_schema.LibraryConfig` + `library_builder.build_library_database` + "
        "`unit_cell_library.build_full_library` (Phase 7, unmodified) build one `database.Database`; "
        "`library_entry.from_database` then registers it into a `UniversalDatabase`. The existing "
        "4.8 mm build is simply the first (and, at this point, only) library built this way.", "",
        "## 14. Validation levels", "",
        "`validation_levels.VALIDATION_LEVELS` = DATA_PRESENT < NUMERICALLY_VALIDATED < "
        "CROSS_SOLVER_VALIDATED < HFSS_VALIDATED < EXPERIMENTALLY_VALIDATED. The benchmark entry's "
        f"actual, measured status is `{entry.validation_status}` (every stored record's own "
        "`checks_passed` diagnostic -- Phase 5, unmodified -- determines this; nothing stronger is "
        "claimed: no cross-solver, HFSS, or experimental comparison exists anywhere in this "
        "repository).", "",
        "## 15. Provenance", "",
        f"Every `UnitCellRecord.provenance` carries `source_type` (`\"MoM\"`), `solver_version` "
        f"(`{entry.solver_version}`, from `manifest.json`), `configuration_hash` "
        f"(`{entry.configuration_hash}`), and its own `validation_status` -- never presented as "
        "validated merely because it exists in the database (Sec. 22).", "",
        "## 16. Integrity rules", "",
        "`integrity.run_integrity_checks` (Sec. 24-25, tests UDB-01..UDB-12): unique record "
        "identity, valid parameter ranges, complete frequency grids, finite complex-valued "
        "responses, preserved cross-polarization, no fabricated values for failed cells, complete "
        "metadata, a content hash, unique library ids, and JSON round-trip exactness.", "",
    ]) + "\n"


def render_validation_report_md(criteria_results, benchmark_summary, integrity_summary, final_answer,
                                 control_variable_status, current_benchmark_library, period_is_global_constraint):
    lines = [
        "# validation_report.md", "",
        "## Acceptance criteria (Sec. 31)", "",
    ]
    for name, passed in criteria_results:
        lines.append(f"- [{'x' if passed else ' '}] {name}")
    lines += [
        "", "## Benchmark query results (Sec. 20)", "",
        benchmark_summary, "",
        "## Integrity check summary (Sec. 24-25)", "",
        integrity_summary, "",
        "## Mandatory final answer (Sec. 32)", "",
        "> Can the project now maintain a universal, period-agnostic, topology-aware unit-cell "
        "electromagnetic database from which an HMSA-specific library can be selected or requested "
        "according to the actual antenna requirements, without forcing the design to use the "
        "existing 4.8-mm library?", "",
        f"**{final_answer}**", "",
        f"- `CONTROL_VARIABLE_STATUS` = `{control_variable_status}` (Phase 8C's own result; this "
        "phase does not touch it, Sec. 28)", "",
        f"- `CURRENT_BENCHMARK_LIBRARY` = `{current_benchmark_library}`", "",
        f"- `PERIOD_AS_GLOBAL_CONSTRAINT` = `{period_is_global_constraint}`", "",
        "No future period is claimed validated unless an actual electromagnetic dataset for it is "
        "registered in the `UniversalDatabase` (Sec. 32's explicit rule).", "",
    ]
    return "\n".join(lines) + "\n"
