"""Sec. 9: the unit-cell-family abstraction, derived from what the registered `LibraryEntry`
objects actually contain -- never a hand-written list of "the" parameters for a topology. A future
family (slot-loaded patch, arm-length dipole, ...) needs no schema change: `parameter_schema`
below is computed from `database.PARAMETER_NAMES`'s actual variance across the entry's own cells.
"""

from dataclasses import dataclass

import numpy as np

from ..database import PARAMETER_NAMES
from ..phase8b.conventions import (
    NO_GROUND_PLANE,
    NORMAL_INCIDENCE_ONLY,
    OBLIQUE_INCIDENCE_PORTS_IMPLEMENTED,
    SURFACE_WAVE_DISPERSION_IMPLEMENTED,
)


@dataclass(frozen=True)
class UnitCellFamily:
    family_id: str
    topology: str
    parameter_schema: tuple  # PARAMETER_NAMES that actually vary across this family's cells
    fixed_parameters: dict  # PARAMETER_NAMES held constant (value, not swept) in this family
    material_stacks: tuple  # distinct {name, eps_r, tan_delta, thickness_mm} dicts observed
    supported_periods_mm: tuple  # distinct (p_x_mm, p_y_mm) pairs actually present
    supported_frequency_range_hz: tuple  # (min, max) across all entries of this topology
    supported_polarizations: tuple
    supported_incidence: tuple
    solver_method: str
    validity_constraints: dict
    geometry_export_method: str
    library_ids: tuple  # entries contributing to this family


def _varying_parameters(database):
    varying, fixed = [], {}
    for name in PARAMETER_NAMES:
        values = np.array([c.parameters[name] for c in database.cells])
        if values.size and np.ptp(values) > 0.0:
            varying.append(name)
        elif values.size:
            fixed[name] = float(values[0])
    return tuple(varying), fixed


def family_from_entries(topology, entries):
    """Aggregate every `LibraryEntry` sharing ``topology`` into one `UnitCellFamily`. Raises
    ValueError if ``entries`` is empty or mixes topologies (a family is topology-scoped by
    definition, Sec. 6-9)."""
    entries = list(entries)
    if not entries:
        raise ValueError("family_from_entries: no entries")
    if any(e.topology != topology for e in entries):
        raise ValueError("family_from_entries: entries must all share the same topology")

    varying_sets, fixed_union = set(), {}
    material_stacks = set()
    periods = set()
    for e in entries:
        v, f = _varying_parameters(e.database)
        varying_sets |= set(v)
        fixed_union.update(f)
        material_stacks.add(
            (e.substrate["name"], e.substrate["eps_r"], e.substrate["tan_delta"], e.substrate["thickness_mm"])
        )
        periods.add((e.p_x_mm, e.p_y_mm))

    freq_min = min(e.frequency_min_hz for e in entries)
    freq_max = max(e.frequency_max_hz for e in entries)
    pols = tuple(sorted(set().union(*(set(e.polarizations_supported) for e in entries))))
    incidence = tuple(sorted(set().union(*(set(e.incidence_supported) for e in entries))))

    return UnitCellFamily(
        family_id=f"family_{topology}",
        topology=topology,
        parameter_schema=tuple(sorted(varying_sets)),
        fixed_parameters=fixed_union,
        material_stacks=tuple(
            {"name": n, "eps_r": e, "tan_delta": t, "thickness_mm": h} for (n, e, t, h) in sorted(material_stacks)
        ),
        supported_periods_mm=tuple(sorted(periods)),
        supported_frequency_range_hz=(freq_min, freq_max),
        supported_polarizations=pols,
        supported_incidence=incidence,
        solver_method="lightweight_mom periodic MoM (Phases 1-7), unmodified",
        validity_constraints={
            "normal_incidence_only": NORMAL_INCIDENCE_ONLY,
            "ground_plane": "ABSENT" if NO_GROUND_PLANE else "PRESENT",
            "oblique_incidence_ports_implemented": OBLIQUE_INCIDENCE_PORTS_IMPLEMENTED,
            "surface_wave_dispersion_implemented": SURFACE_WAVE_DISPERSION_IMPLEMENTED,
        },
        geometry_export_method="cad_export.py STEP export (Phase 6, unmodified)",
        library_ids=tuple(e.library_id for e in entries),
    )
