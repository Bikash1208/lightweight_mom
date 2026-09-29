"""Phase 9A: physical control-variable validation. Reuses Phase 7's own solved data (never
re-solved) and Phase 8B/8C's architecture findings (never re-derived); compares against HFSS data
that -- as of this run -- does not exist anywhere in this project (`hfss_reference.py`), so every
production comparison honestly reports `NOT_AVAILABLE_FROM_CURRENT_DATA` rather than a fabricated
number (Rule 1, Sec. 45).
"""

from dataclasses import dataclass

from ..phase8b.conventions import COORDINATE_SYSTEM, PORT_DEFINITION
from ..phase8c import topology_audit as topo_mod
from .error_metrics import compute_error_stats

# Sec. 7's minimum validation subset (mm).
REFERENCE_CELLS_MM = (
    ("CELL-01", 1.0, 1.0),
    ("CELL-02", 1.0, 4.4),
    ("CELL-03", 4.4, 1.0),
    ("CELL-04", 4.4, 4.4),
    ("CELL-05", 2.0, 2.0),
    ("CELL-06", 2.8, 2.8),
    ("CELL-07", 3.6, 2.0),
    ("CELL-08", 2.0, 3.6),
    ("CELL-09", 3.2, 2.4),
    ("CELL-10", 2.4, 3.2),
)

# Sec. 11's minimum frequency set -- every value below is an exact point of Phase 7's 10-15 GHz /
# 25 MHz grid (verified by test), so no interpolation is ever needed for this specific set.
VALIDATION_FREQUENCIES_GHZ = (10.0, 11.0, 12.0, 12.5, 13.0, 14.0, 15.0)

CONTROL_VARIABLE_STATES = (
    "VALIDATED",
    "PARTIALLY_VALIDATED",
    "INCONSISTENT",
    "UNVALIDATED",
    "ARCHITECTURE_INCOMPATIBLE",
)


@dataclass(frozen=True)
class SolverReferenceRow:
    cell_id: str
    l_mm: float
    w_mm: float
    frequency_ghz: float
    s11: complex
    s21: complex


def extract_solver_reference_data(database, cells_mm=REFERENCE_CELLS_MM, frequencies_ghz=VALIDATION_FREQUENCIES_GHZ):
    """Sec. 7, 11: solver S11xx/S21xx for the 10 reference cells at the 7 validation frequencies,
    read directly from the real, already-solved Phase 7 checkpoint (`database.Database`,
    unmodified) -- 'use the existing database records rather than creating duplicate records'
    (Sec. 7). Raises KeyError if a requested (L, W) or frequency is not in the database (never
    silently substituted or interpolated)."""
    freqs_hz = database.frequencies
    rows = []
    for cell_id, l_mm, w_mm in cells_mm:
        cell = next(
            (
                c
                for c in database.cells
                if abs(c.parameters["length_m"] * 1e3 - l_mm) < 1e-6
                and abs(c.parameters["width_m"] * 1e3 - w_mm) < 1e-6
            ),
            None,
        )
        if cell is None or cell.status != "ok":
            raise KeyError(f"{cell_id}: (L={l_mm}, W={w_mm}) mm not found as an OK cell in the database")
        for f_ghz in frequencies_ghz:
            f_hz = f_ghz * 1e9
            matches = [i for i, f in enumerate(freqs_hz) if abs(f - f_hz) <= 1e-6 * f_hz]
            if not matches:
                raise KeyError(f"{f_ghz} GHz is not an exact grid point of this database")
            record = database.record(cell.index, matches[0])
            rows.append(
                SolverReferenceRow(
                    cell_id=cell_id,
                    l_mm=l_mm,
                    w_mm=w_mm,
                    frequency_ghz=f_ghz,
                    s11=record.component("S11_xx"),
                    s21=record.component("S21_xx"),
                )
            )
    return rows


def architecture_audit(database, mean_power):
    """Sec. 8, 14: documents this project's OWN configuration/architecture (reusing
    `phase8c.topology_audit`, `phase8b.conventions`, unmodified) and explicitly marks every
    HFSS-specific configuration field `UNKNOWN` -- no HFSS project exists in this repository, so
    none of those fields can be determined (Rule 1: 'do not guess')."""
    audit = topo_mod.audit_topology(database, mean_power)
    return {
        "hfss_solution_type": "UNKNOWN (no HFSS project exists in this repository)",
        "hfss_boundary_conditions": "UNKNOWN",
        "hfss_periodic_boundaries": "UNKNOWN",
        "hfss_floquet_ports": "UNKNOWN",
        "hfss_number_of_floquet_modes": "UNKNOWN",
        "hfss_mode_polarization": "UNKNOWN",
        "hfss_reference_impedance": "UNKNOWN",
        "hfss_deembedding_reference_plane": "UNKNOWN",
        "hfss_coordinate_system": "UNKNOWN",
        "hfss_mesh_convergence": "UNKNOWN",
        "solver_ground_plane": audit.ground_plane,
        "solver_ground_plane_evidence": audit.ground_plane_evidence,
        "solver_architecture_classification": audit.architecture_classification,
        "solver_architecture_evidence": audit.architecture_evidence,
        "solver_incident_wave": PORT_DEFINITION.port1_incident_direction + " at " + PORT_DEFINITION.port1_plane,
        "solver_reflected_wave": "specular (0,0) Floquet order returning through port 1 (top reference plane)",
        "solver_transmitted_wave": "specular (0,0) Floquet order continuing through port 2 (bottom reference plane, "
        + PORT_DEFINITION.port2_plane
        + ")",
        "solver_holographic_surface_location": COORDINATE_SYSTEM.pattern_plane,
        "solver_port_modes": PORT_DEFINITION.modes,
        "solver_coordinate_system": COORDINATE_SYSTEM.z_convention,
    }


def classify_control_variable(stats):
    """Sec. 13: one of `CONTROL_VARIABLE_STATES`, never a numerical winner score.
    ``stats``: `error_metrics.ComplexErrorStats`, or `None` if no HFSS data was available."""
    if stats is None:
        return "UNVALIDATED", "NOT_AVAILABLE_FROM_CURRENT_DATA: no HFSS reference file found (hfss_reference.py)"
    if abs(stats.mean_abs_phase_error_deg - 180.0) < 5.0:
        return "INCONSISTENT", (
            "mean phase error is within 5 deg of 180 deg across every compared point -- an "
            "objective sign/convention mismatch (e.g. time-convention or port-direction), not a "
            "matter of numerical agreement, and not something this module corrects automatically "
            "(Rule 2)."
        )
    return "UNVALIDATED", (
        "solver-vs-HFSS error statistics were computed in full (see the accompanying "
        "ComplexErrorStats), but no project document or cited scientific basis defines an "
        "acceptable cross-solver agreement tolerance for this application (Sec. 24); "
        "VALIDATED/PARTIALLY_VALIDATED cannot be honestly asserted without one."
    )


def compare_against_hfss(hfss_records, solver_rows, component):
    """Sec. 9-11: pairs solver rows with HFSS records by (cell_id, frequency_ghz) and computes
    `error_metrics.compute_error_stats` over every matched point. Raises ValueError if a solver
    row has no matching HFSS record (never silently drops or substitutes a point)."""
    hfss_by_key = {(r.cell_id, round(r.frequency_ghz, 6)): r for r in hfss_records}
    solver_values, hfss_values = [], []
    for row in solver_rows:
        key = (row.cell_id, round(row.frequency_ghz, 6))
        if key not in hfss_by_key:
            raise ValueError(f"no HFSS record for {key}")
        hfss_row = hfss_by_key[key]
        solver_values.append(getattr(row, component))
        hfss_values.append(getattr(hfss_row, component))
    return compute_error_stats(solver_values, hfss_values)
