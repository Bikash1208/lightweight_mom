"""Sec. 14, 20: the physical-filter-first query engine. Filters are applied strictly in the order
Sec. 14 lists (topology, material stack, frequency, period, polarization, incidence, response
availability, validity status); nothing here ranks by phase coverage or any other controllability
metric (Sec. 14's explicit "do not rank cells merely by phase coverage").
"""

from dataclasses import dataclass

_RTOL = 1e-6


@dataclass(frozen=True)
class UnitCellQuery:
    topology: str | None = None
    substrate_name: str | None = None
    eps_r: float | None = None
    tan_delta: float | None = None
    substrate_thickness_mm: float | None = None
    backing: str | None = None           # "air" | "pec" (Phase 12A: never mix grounded/ungrounded)
    p_x_mm: float | None = None
    p_y_mm: float | None = None
    frequency_hz: float | None = None
    polarization: str | None = None      # "x" | "y"
    incidence: str | None = None         # "normal"
    min_validation_status: str | None = None   # e.g. "NUMERICALLY_VALIDATED"


@dataclass(frozen=True)
class QueryResult:
    status: str                          # "MATCH" | "DATABASE_COVERAGE_MISSING"
    query: UnitCellQuery
    matching_entries: tuple              # LibraryEntry objects that survived every filter
    filter_trail: tuple                  # (filter_name, n_candidates_remaining) per stage, in order
    eliminated_at: str | None            # which filter first reduced candidates to zero, if any


def _close(a, b, rtol=_RTOL):
    return abs(a - b) <= rtol * max(abs(a), abs(b), 1.0)


_FILTERS = (
    ("topology", lambda e, q: q.topology is None or e.topology == q.topology),
    ("material_stack", lambda e, q: (
        (q.substrate_name is None or e.substrate["name"] == q.substrate_name)
        and (q.eps_r is None or _close(e.substrate["eps_r"], q.eps_r))
        and (q.tan_delta is None or _close(e.substrate["tan_delta"], q.tan_delta))
        and (q.substrate_thickness_mm is None or _close(e.substrate["thickness_mm"], q.substrate_thickness_mm))
        # Phase 12A: an explicit backing request is never satisfied by a library built for the
        # other backing -- a grounded (backing="pec") and an ungrounded (backing="air") library
        # are never silently mixed in a search result.
        and (q.backing is None or e.backing == q.backing)
    )),
    ("frequency", lambda e, q: q.frequency_hz is None or e.has_frequency(q.frequency_hz)),
    ("period", lambda e, q: (q.p_x_mm is None and q.p_y_mm is None) or e.matches_period(
        q.p_x_mm if q.p_x_mm is not None else e.p_x_mm,
        q.p_y_mm if q.p_y_mm is not None else e.p_y_mm,
    )),
    ("polarization", lambda e, q: q.polarization is None or q.polarization in e.polarizations_supported),
    ("incidence", lambda e, q: q.incidence is None or q.incidence in e.incidence_supported),
    ("response_availability", lambda e, q: e.n_ok_cells > 0),
    ("validity_status", lambda e, q: q.min_validation_status is None or _meets_validation(e, q.min_validation_status)),
)

_VALIDATION_ORDER = (
    "DATA_PRESENT", "NUMERICALLY_VALIDATED", "CROSS_SOLVER_VALIDATED", "HFSS_VALIDATED",
    "EXPERIMENTALLY_VALIDATED",
)


def _meets_validation(entry, minimum):
    return _VALIDATION_ORDER.index(entry.validation_status) >= _VALIDATION_ORDER.index(minimum)


def search_unit_cells(universal_db, query: UnitCellQuery):
    """Sec. 14, 20-21: apply every filter in `_FILTERS`' fixed order over
    ``universal_db.entries()``. Returns a `QueryResult`; if candidates reach zero at any stage,
    the result is `DATABASE_COVERAGE_MISSING` and names the exact stage (Sec. 21: never silently
    substitute, rescale or interpolate a near match)."""
    candidates = list(universal_db.entries())
    trail = [("initial", len(candidates))]
    eliminated_at = None
    for name, predicate in _FILTERS:
        candidates = [e for e in candidates if predicate(e, query)]
        trail.append((name, len(candidates)))
        if not candidates and eliminated_at is None:
            eliminated_at = name
    status = "MATCH" if candidates else "DATABASE_COVERAGE_MISSING"
    return QueryResult(
        status=status, query=query, matching_entries=tuple(candidates),
        filter_trail=tuple(trail), eliminated_at=eliminated_at,
    )
