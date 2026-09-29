"""Phase 8A Sec. 29-31: reference-grid / MTS-grid sampling compatibility and spatial-mapping
availability.

This project has no validated reference-field-coordinate -> MTS-unit-cell-centre mapping method
(grep-verified: no `resample`/`mapping`/`aperture` module anywhere in `src/`). `Section 31`
explicitly forbids assuming nearest-neighbour (or any other) mapping without implementing and
validating it as its own feature, so `check_mapping_method_defined` always reports "not defined"
for this project as it stands -- a real limitation, not a placeholder bug.
"""

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class SamplingCompatibility:
    reference_dx_mm: float | None
    reference_dy_mm: float | None
    period_x_mm: float
    period_y_mm: float
    grids_identical: bool
    note: str


def check_sampling_compatibility(reference_dx_mm, reference_dy_mm, period_x_mm, period_y_mm, rtol=1e-6):
    """Sec. 30: identifies whether the reference-field sample spacing equals the unit-cell period
    -- never silently resamples either way."""
    identical = (
        reference_dx_mm is not None and reference_dy_mm is not None
        and abs(reference_dx_mm - period_x_mm) <= rtol * period_x_mm
        and abs(reference_dy_mm - period_y_mm) <= rtol * period_y_mm
    )
    note = (
        "Reference-field sampling matches the unit-cell period exactly."
        if identical else
        "Reference grid and MTS grid are non-identical. A validated spatial resampling/mapping "
        "method is required before final geometry generation."
    )
    return SamplingCompatibility(reference_dx_mm, reference_dy_mm, period_x_mm, period_y_mm, identical, note)


@dataclass(frozen=True)
class MappingAvailability:
    defined: bool
    note: str


def check_mapping_method_defined():
    """Sec. 31: this project's actual state -- no mapping method exists, checked by inspection
    (see module docstring), not assumed."""
    return MappingAvailability(
        defined=False,
        note=(
            "No validated reference-field-coordinate -> MTS-unit-cell-centre mapping method "
            "exists in this project (no resampling/aperture-mapping module in src/). "
            "Nearest-neighbour or any other mapping must not be assumed without being "
            "implemented and validated as its own feature (Sec. 31)."
        ),
    )


@dataclass(frozen=True)
class ApertureEstimate:
    max_nx: int | None
    max_ny: int | None
    note: str


def estimate_max_aperture(reference_extent_x_mm, reference_extent_y_mm, period_x_mm, period_y_mm):
    """Sec. 29: an integer floor-division upper bound only -- Phase 8A does not choose the final
    aperture (that belongs to a later aperture-selection phase)."""
    if reference_extent_x_mm is None or reference_extent_y_mm is None:
        return ApertureEstimate(None, None, "reference-field extent unknown; cannot estimate aperture")
    nx = int(math.floor(reference_extent_x_mm / period_x_mm))
    ny = int(math.floor(reference_extent_y_mm / period_y_mm))
    return ApertureEstimate(
        nx, ny,
        f"floor(extent/period): Nx={nx}, Ny={ny} (upper bound only; final aperture selection is a later phase)",
    )
