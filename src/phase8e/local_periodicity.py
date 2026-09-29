"""Sec. 15: the local-periodic approximation (unit-cell response computed in isolation under a
plane wave stands in for its response under the actual, spatially-varying HMSA field). This
project's solver has no coupled-multi-cell or near-field-array solve capability -- checked, not
assumed -- so this module can never report the approximation as validated; it always reports
`LOCAL_PERIODICITY_UNVALIDATED`, alongside the one supporting number it *can* honestly compute:
how many unit cells span one fringe of the fastest required spatial variation (Sec. 15's "period
relative to spatial variation").
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class LocalPeriodicityAssessment:
    period_m: float
    k_t_dominant_per_m: float | None  # dominant transverse wavenumber driving the modulation
    cells_per_fringe: float | None  # (2*pi/k_t_dominant) / period, None if k_t_dominant is 0/None
    status: str  # always "LOCAL_PERIODICITY_UNVALIDATED"
    note: str


def assess_local_periodicity(period_m, k_t_dominant_per_m=None):
    """Sec. 15: never claims validity merely because the cell is electrically small (this
    project's period, 4.8 mm at 12.5 GHz, is ~0.2 lambda -- small, but smallness alone is not
    evidence, Sec. 15's explicit rule). No coupled-cell/near-field-array solver exists anywhere in
    `src/` (grep-verified: no 'mutual coupling'/'array factor'/'finite array' module), so this
    function can never report anything but `LOCAL_PERIODICITY_UNVALIDATED`."""
    cells_per_fringe = None
    if k_t_dominant_per_m:
        import math

        fringe_period_m = 2.0 * math.pi / k_t_dominant_per_m
        cells_per_fringe = fringe_period_m / period_m
    return LocalPeriodicityAssessment(
        period_m=period_m,
        k_t_dominant_per_m=k_t_dominant_per_m,
        cells_per_fringe=cells_per_fringe,
        status="LOCAL_PERIODICITY_UNVALIDATED",
        note=(
            "No coupled-multi-cell or near-field-array solve exists in this project to validate "
            "the local-periodic approximation; "
            + (
                f"the fastest required spatial variation spans approximately "
                f"{cells_per_fringe:.2f} unit cells per fringe (informative only, not a validity test)."
                if cells_per_fringe is not None
                else "no dominant transverse wavenumber was supplied to even compute the sampling-density number."
            )
        ),
    )
