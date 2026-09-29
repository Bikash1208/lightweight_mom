"""Sec. 22-23: provenance and validation-level taxonomy.

Five explicit levels, most to least demanding. A `LibraryEntry`/record must never be reported as
a stronger level than what has actually been done in this repository -- and today, nothing in
this repository has ever been cross-checked against an independent solver, HFSS, or a physical
measurement (grep-confirmed: no CST/HFSS/measurement data file exists anywhere in the project;
`hfss_comparison.py`/`library_validation.py` only *export* MoM data for a future HFSS comparison
that has not happened). Only the two weakest levels are therefore ever actually assigned by this
package; the other three are represented so the schema does not need to change the day a real
cross-solver or experimental comparison is finally run.
"""

from dataclasses import dataclass

VALIDATION_LEVELS = (
    "DATA_PRESENT",
    "NUMERICALLY_VALIDATED",
    "CROSS_SOLVER_VALIDATED",
    "HFSS_VALIDATED",
    "EXPERIMENTALLY_VALIDATED",
)


@dataclass(frozen=True)
class Provenance:
    source_type: str  # "MoM" | "HFSS" | "CST" | "measurement" | "other"
    source_version: str | None
    solver_version: str | None
    configuration_hash: str | None
    date: str | None
    validation_status: str

    def __post_init__(self):
        if self.validation_status not in VALIDATION_LEVELS:
            raise ValueError(f"validation_status {self.validation_status!r} not in {VALIDATION_LEVELS}")


def record_validation_status(diagnostics):
    """Sec. 23: NUMERICALLY_VALIDATED iff this specific record's own self-consistency diagnostics
    (`database.DIAGNOSTIC_TYPES`: reciprocity_error, passivity_min_eig, incoming_residual,
    cells_per_wavelength -- Phase 5, unmodified) report ``checks_passed`` True; DATA_PRESENT
    otherwise. This is a real, per-record, data-driven determination (the diagnostic already
    exists in every stored record), not a blanket label assigned to the whole library."""
    if diagnostics.get("checks_passed"):
        return "NUMERICALLY_VALIDATED"
    return "DATA_PRESENT"


def library_validation_status(records):
    """The validation level a whole `LibraryEntry` may honestly claim: NUMERICALLY_VALIDATED only
    if *every* record in it passed its own checks, else DATA_PRESENT -- never averaged or rounded
    up. Raises ValueError on an empty record set (there is nothing to validate)."""
    records = list(records)
    if not records:
        raise ValueError("library_validation_status: no records to evaluate")
    statuses = {record_validation_status(r.diagnostics) for r in records}
    if statuses == {"NUMERICALLY_VALIDATED"}:
        return "NUMERICALLY_VALIDATED"
    return "DATA_PRESENT"
