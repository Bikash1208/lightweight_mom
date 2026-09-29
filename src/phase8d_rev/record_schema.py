"""Sec. 7-8: the mandatory per-response record schema, built directly from an already-solved
`database.Database`'s `CellEntry`/`Record` (Phase 5, unmodified) -- nothing here re-solves,
re-derives an S-parameter, or invents a value. The complex S-matrix entry (Sec. 8's "primary
electromagnetic datum") is copied verbatim from `database.Record.s`; magnitude/phase/power are
derived exactly as `database.py`/`phase8c.power_balance` already compute them.
"""

from dataclasses import dataclass

from ..database import PARAMETER_NAMES, component_names
from ..phase8c.power_balance import power_balance_for_record
from .validation_levels import Provenance, record_validation_status

# The two PARAMETER_NAMES this schema promotes to first-class geometry fields (Sec. 7's own
# example); the rest are carried generically in `extra_geometry_parameters` so a future unit-cell
# family (slot_length, arm_length, ...) never needs a schema change (Sec. 9's explicit requirement).
_PRIMARY_GEOMETRY = ("length_m", "width_m")
_EXTRA_GEOMETRY = tuple(n for n in PARAMETER_NAMES if n not in _PRIMARY_GEOMETRY)

LINEAR_COMPONENT_NAMES = component_names("linear")   # 16 names, e.g. "S11_xx"


@dataclass(frozen=True)
class UnitCellRecord:
    # --- identity ---------------------------------------------------------
    record_id: str
    library_id: str
    topology_id: str
    simulation_id: str

    # --- geometry -----------------------------------------------------------
    topology: str
    p_x_mm: float
    p_y_mm: float
    L_mm: float
    W_mm: float
    extra_geometry_parameters: dict     # e.g. {"center_x_m": ..., "center_y_m": ..., "rotation_deg": ...}

    # --- materials ------------------------------------------------------------
    substrate_name: str
    epsilon_r: float
    tan_delta: float
    substrate_thickness_mm: float
    additional_layers: tuple            # empty in this solver; never fabricated

    # --- electromagnetic conditions --------------------------------------------
    frequency_hz: float
    incidence_theta_deg: float
    incidence_phi_deg: float
    polarization_basis: str             # "linear (x, y)" -- both incident polarizations are stored
    propagation_direction: str

    # --- boundary conditions ----------------------------------------------------
    periodic_x: bool
    periodic_y: bool
    port_type: str
    floquet_mode_definition: str

    # --- scattering response (Sec. 8: complex value is primary) -----------------
    s_parameters: dict                  # name -> {re, im, mag, phase_deg} for all 16 linear entries

    # --- derived quantities (Sec. 7, computed from s_parameters, never stored instead of it) -----
    reflected_power_co: float
    reflected_power_cross: float
    transmitted_power_co: float
    transmitted_power_cross: float
    absorbed_power_fraction: float
    polarization_conversion_note: str
    phase_coverage_metadata: None       # library-level quantity; never fabricated per-record (Sec. 8)

    # --- provenance / validity (Sec. 22-23) --------------------------------------
    provenance: Provenance


def record_from_database(database, library_id, topology, cell, record, extra=None):
    """Build one `UnitCellRecord` from an existing `database.Database`'s `CellEntry`/`Record`
    (Phase 5/7, unmodified). ``extra``: an optional dict of additional identity fields
    (``configuration_hash``, ``solver_version``, ``source_type``) supplied by the caller when
    known (e.g. from `manifest.json`); left ``None`` where the caller has not supplied one --
    never guessed here."""
    extra = extra or {}
    ctx = database.metadata.get("context") or {}
    lattice = ctx.get("lattice", {})
    substrate = ctx.get("substrate", {})

    s_entries = {}
    for name in LINEAR_COMPONENT_NAMES:
        z = record.component(name)
        s_entries[name] = {
            "re": float(z.real), "im": float(z.imag), "mag": float(abs(z)),
            "phase_deg": float(record.phase_deg(name)),
        }

    balance = power_balance_for_record(record)
    extra_geometry = {
        name: cell.parameters[name] for name in _EXTRA_GEOMETRY
    }

    record_id = f"{library_id}:{cell.cell_id}:{record.frequency_index}"
    provenance = Provenance(
        source_type=extra.get("source_type", "MoM"),
        source_version=extra.get("source_version"),
        solver_version=extra.get("solver_version"),
        configuration_hash=extra.get("configuration_hash"),
        date=extra.get("date"),
        validation_status=record_validation_status(record.diagnostics),
    )

    return UnitCellRecord(
        record_id=record_id, library_id=library_id, topology_id=topology, simulation_id=record_id,
        topology=topology,
        p_x_mm=float(lattice.get("period_x_m", float("nan")) * 1e3) if "period_x_m" in lattice else float("nan"),
        p_y_mm=float(lattice.get("period_y_m", float("nan")) * 1e3) if "period_y_m" in lattice else float("nan"),
        L_mm=round(cell.parameters["length_m"] * 1e3, 9),
        W_mm=round(cell.parameters["width_m"] * 1e3, 9),
        extra_geometry_parameters=extra_geometry,
        substrate_name=substrate.get("name", "unavailable"),
        epsilon_r=substrate.get("eps_r", float("nan")),
        tan_delta=substrate.get("tan_delta", float("nan")),
        substrate_thickness_mm=float(substrate["thickness_m"] * 1e3) if "thickness_m" in substrate else float("nan"),
        additional_layers=(),
        frequency_hz=record.frequency,
        incidence_theta_deg=0.0,        # NORMAL_INCIDENCE_ONLY (phase8b.conventions) -- solver-wide limit
        incidence_phi_deg=0.0,
        polarization_basis="linear (x, y); both incident polarizations stored (database.CONVENTIONS)",
        propagation_direction="port 1 (-z) and port 2 (+z) incidence, two-port S-matrix",
        periodic_x=True, periodic_y=True,
        port_type="two reference planes, specular (0,0) Floquet order only (phase8b.conventions.PORT_DEFINITION)",
        floquet_mode_definition="kxm = kx0 + 2*pi*m/Px, kyn = ky0 + 2*pi*n/Py (periodic_bc.RectangularLattice)",
        s_parameters=s_entries,
        reflected_power_co=balance.R_co, reflected_power_cross=balance.R_cross,
        transmitted_power_co=balance.T_co, transmitted_power_cross=balance.T_cross,
        absorbed_power_fraction=balance.absorbed_fraction,
        polarization_conversion_note=balance.note,
        phase_coverage_metadata=None,
        provenance=provenance,
    )
