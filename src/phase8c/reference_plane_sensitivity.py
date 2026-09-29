"""Phase 8C Sec. 29: reference-plane sensitivity of the candidate scattering-parameter phase.

`CONVENTIONS.md` ("Waves, ports and S-parameters") already documents the exact effect of moving
the S-parameter reference planes:

    S11' = S11 * exp(-2j k0 D1),   S22' = S22 * exp(+2j k0 D2),
    S21' = S12' = S21 * exp(+j k0 (D2 - D1))

(D1, D2 = displacement of the new top/bottom planes from the stored ones). This module reuses that
exact formula (not a new derivation) purely to *demonstrate* the sensitivity -- it never applies it
to the stored database (Sec. 29's explicit "do not arbitrarily remove a linear phase term / do not
introduce an arbitrary phase offset").
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber
from ..phase import wrapped_difference_deg


@dataclass(frozen=True)
class ReferencePlaneSensitivity:
    frequency_hz: float
    shift_m: float
    phase_before_deg: float
    phase_after_deg_S11: float
    phase_after_deg_S21: float
    delta_deg_S11: float
    delta_deg_S21: float
    note: str


def demonstrate_sensitivity(s11_value, s21_value, frequency_hz, shift_m=1.0e-3):
    """Applies the `CONVENTIONS.md` moving-plane formula to a *copy* of one already-stored (S11,
    S21) pair, for a purely illustrative top-plane shift D1=``shift_m`` (D2=0) -- proves phase is
    reference-plane dependent; changes nothing in any stored database."""
    k0 = free_space_wavenumber(frequency_hz)
    s11_after = s11_value * np.exp(-2j * k0 * shift_m)
    s21_after = s21_value * np.exp(1j * k0 * (0.0 - shift_m))
    phase_before = float(np.degrees(np.angle(s11_value)))
    phase_after_s11 = float(np.degrees(np.angle(s11_after)))
    phase_after_s21 = float(np.degrees(np.angle(s21_after)))
    delta_s11 = float(wrapped_difference_deg(phase_after_s11, phase_before))
    delta_s21 = float(wrapped_difference_deg(phase_after_s21, float(np.degrees(np.angle(s21_value)))))
    return ReferencePlaneSensitivity(
        frequency_hz=frequency_hz,
        shift_m=shift_m,
        phase_before_deg=phase_before,
        phase_after_deg_S11=phase_after_s11,
        phase_after_deg_S21=phase_after_s21,
        delta_deg_S11=delta_s11,
        delta_deg_S21=delta_s21,
        note=(
            f"A {shift_m * 1e3:g} mm top-reference-plane shift alone changes S11 phase by "
            f"{delta_s11:.2f} deg and S21 phase by {delta_s21:.2f} deg at "
            f"{frequency_hz / 1e9:g} GHz (`CONVENTIONS.md` moving-plane formula, reused, not "
            "applied to any stored data). Phase is therefore only meaningful relative to the "
            "specific reference planes recorded in the database metadata "
            "(`context['reference_planes']`); it cannot be used as a geometry-control law without "
            "holding that reference fixed -- which every Phase 7 build already does consistently."
        ),
    )
