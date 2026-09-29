"""Phase 8A Sec. 28: grating-lobe / period feasibility check.

Reuses exactly the formula `discretization.grating_onset_frequency` uses (f = c0 / max(Px, Py)),
inlined here only because that function needs a solved `UnitCell` object and this gate works from
the library's own stored `period_x_mm`/`period_y_mm` (no re-solve, no new criterion) --
`tests/test_feasibility_grating.py` cross-checks the two give identical numbers for the locked
Phase 7 period, so this is verified reuse, not an independent re-derivation (brief: "do not invent
a new grating criterion").
"""

from dataclasses import dataclass

from .constants import C0


@dataclass(frozen=True)
class GratingCheckResult:
    period_x_m: float
    period_y_m: float
    frequency_hz: float
    wavelength_m: float
    period_over_wavelength: float
    grating_onset_hz: float
    is_grating_risk: bool


def grating_check(period_x_m, period_y_m, frequency_hz):
    if period_x_m <= 0.0 or period_y_m <= 0.0:
        raise ValueError(f"periods must be > 0, got ({period_x_m!r}, {period_y_m!r})")
    if frequency_hz <= 0.0:
        raise ValueError(f"frequency_hz must be > 0, got {frequency_hz!r}")
    wavelength_m = C0 / frequency_hz
    onset_hz = C0 / max(period_x_m, period_y_m)      # discretization.grating_onset_frequency's formula
    return GratingCheckResult(
        period_x_m=period_x_m, period_y_m=period_y_m, frequency_hz=frequency_hz,
        wavelength_m=wavelength_m,
        period_over_wavelength=max(period_x_m, period_y_m) / wavelength_m,
        grating_onset_hz=onset_hz,
        is_grating_risk=bool(frequency_hz >= onset_hz),
    )
