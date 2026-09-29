"""Sec. 16-17: period-independent analysis. Reuses `phase8d_rev.candidate_period` (the Nyquist-
sampling candidate-period heuristic and the CANDIDATE-vs-VALIDATED distinction) and
`phase8d_rev.build_request`/`query` (unmodified) -- this module adds no new period logic of its
own, only the Sec. 16 "CURRENT_DATASET_LIMITATION vs REQUIRED_LIBRARY_SPECIFICATION" framing.
"""

from dataclasses import dataclass

from ..phase8d_rev.build_request import generate_build_request
from ..phase8d_rev.candidate_period import candidate_periods_from_sampling, validate_period
from ..phase8d_rev.query import UnitCellQuery, search_unit_cells


@dataclass(frozen=True)
class PeriodRequirementAssessment:
    registered_periods_mm: tuple
    candidate_periods_mm: tuple            # from the sampling heuristic, all CANDIDATE_HEURISTIC
    registered_period_matches_candidates: bool | None
    status: str                             # "CURRENT_DATASET_LIMITATION" | "REQUIRED_LIBRARY_MISSING"
    required_library_specification: dict | None


def assess_period_requirement(universal_db, topology, k_required_max_per_m, frequency_hz, substrate_name=None):
    registered = universal_db.periods_for_topology(topology)
    candidates = candidate_periods_from_sampling(k_required_max_per_m, n_candidates=5) if k_required_max_per_m else ()
    candidate_p_mm = tuple(round(c["p_mm"], 6) for c in candidates)

    matches = None
    if registered and candidate_p_mm:
        p_min, p_max = min(candidate_p_mm), max(candidate_p_mm)
        matches = any(p_min * 0.5 <= p_x <= p_max * 2.0 for p_x, _ in registered)   # order-of-magnitude check only

    required_spec = None
    status = "CURRENT_DATASET_LIMITATION"
    query = UnitCellQuery(topology=topology, substrate_name=substrate_name, frequency_hz=frequency_hz)
    if registered:
        # is the registered period actually VALIDATED for this exact frequency/substrate?
        p_x, p_y = registered[0]
        validated = validate_period(universal_db, topology, p_x, p_y, frequency_hz=frequency_hz, substrate_name=substrate_name)
        if validated.status != "VALIDATED" or matches is False:
            status = "REQUIRED_LIBRARY_MISSING"
    else:
        status = "REQUIRED_LIBRARY_MISSING"

    if status == "REQUIRED_LIBRARY_MISSING":
        result = search_unit_cells(universal_db, query)
        if result.status == "DATABASE_COVERAGE_MISSING":
            req = generate_build_request(result)
            required_spec = {
                "topology": req.topology, "p_x_mm": candidate_p_mm[0] if candidate_p_mm else None,
                "p_y_mm": candidate_p_mm[0] if candidate_p_mm else None,
                "frequency_hz": req.frequency_start_hz, "substrate": req.material_stack,
                "polarization": req.polarization, "incidence": req.incidence,
                "geometry_parameters": req.geometry_parameter_ranges,
                "response_required": "see control_response.py / complex_controllability.py for this run",
                "reason": req.reason,
            }

    return PeriodRequirementAssessment(
        registered_periods_mm=registered, candidate_periods_mm=candidate_p_mm,
        registered_period_matches_candidates=matches, status=status,
        required_library_specification=required_spec,
    )
