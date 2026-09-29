"""Phase 8B Sec. 19-20: the holographic modulation/interference quantity.

**Two distinct, literature-sourced candidate formulations are implemented and reported honestly
side by side -- neither is presented as "the" project formulation, because no project document
authorizes either one (Sec. 1, 19's explicit "do not select an equation merely because it looks
like a standard holographic formula").**

## Formulation A: interference-modulation index (Fong et al. 2010)

    Psi(x, y) = Re{ E_ref(x, y) . conj(E_obj(x, y)) }

the interference-pattern kernel of holographic *impedance-modulated leaky-wave* antennas
(`source_registry.json` `fong2010`), where the physical control variable is a **grounded**
surface's reactance, `X(x,y) = X_avg * (1 + M * Psi(x,y))` (M = modulation index, X_avg = average
reactance -- both are *fabrication calibration constants* obtained from a unit-cell library that
characterizes a grounded structure's normal-incidence reflection phase; Phase 7's library does not,
see `control_variable.py`). This module implements `interference_modulation` (the Psi kernel
itself, dimensionless, numerically testable) and `reactance_from_modulation` (the X(x,y) mapping,
implemented as a pure, parameterized function -- explicitly marked unvalidated, since M/X_avg
cannot be calibrated from this project's data).

## Formulation B: phase-difference hologram (reflectarray/transmitarray-style)

    required_phase(x, y) = wrap( arg(E_obj(x, y)) - arg(E_ref(x, y)) )

the simpler, also well-established technique used by passive reflectarray/transmitarray "phase
hologram" designs (each unit cell imposes, via its normal-incidence reflection or transmission
phase, exactly the phase needed to turn the locally-incident reference wave into the locally-
desired object wave). Its physical control variable is a **normal-incidence co-polarized
reflection or transmission phase** -- the same *kind* of quantity Phase 7's library actually
stores (S11xx or S21xx), which is why `control_variable.py` investigates this formulation as the
more promising candidate. This kernel needs only `arg(E_ref)`/`arg(E_obj)`, so it is the one case
where a phase-only reference field is sufficient (`reference_wave.phase_only_reference_requires_amplitude_assumption`).

Status of both: **MATHEMATICALLY DEFINED, NUMERICALLY IMPLEMENTED, NUMERICALLY TESTED** (standard,
convention-consistent physics). **NOT physically validated as this project's intended
formulation** -- no project document picks between them, and neither's physical control variable
has been shown to correspond to what Phase 7's library actually measured (that is
`control_variable.py`'s job, not this module's).
"""

import numpy as np

from ..phase import wrap_deg


def interference_modulation(e_ref, e_obj):
    """Formulation A's kernel: `Psi = Re{E_ref * conj(E_obj)}`, real-valued, same shape as the
    (broadcastable) complex inputs. Dimensionless only if `E_ref`/`E_obj` are both normalized to
    unit amplitude elsewhere; this function performs no implicit normalization (Sec. 42)."""
    e_ref, e_obj = np.asarray(e_ref, dtype=complex), np.asarray(e_obj, dtype=complex)
    return (e_ref * np.conj(e_obj)).real


def reactance_from_modulation(psi, x_avg_ohm, modulation_index):
    """Formulation A's X(x,y) = X_avg (1 + M Psi) mapping -- a pure, documented function of
    already-known M/X_avg, **not** a claim that M/X_avg are known for this project (they are not;
    see `control_variable.py`'s surface-impedance investigation, Sec. 24)."""
    if not (0.0 < modulation_index < 1.0):
        raise ValueError(f"modulation_index must be in (0, 1), got {modulation_index!r}")
    psi = np.asarray(psi, dtype=float)
    return x_avg_ohm * (1.0 + modulation_index * psi)


def phase_difference_hologram(e_ref, e_obj):
    """Formulation B's kernel: `required_phase = wrap(arg(E_obj) - arg(E_ref))` [deg], the
    standard phase-conjugate/phase-difference hologram equation. Accepts complex fields (uses only
    their phase) so it works identically whether the reference is phase-only or complex."""
    e_ref, e_obj = np.asarray(e_ref, dtype=complex), np.asarray(e_obj, dtype=complex)
    return wrap_deg(np.degrees(np.angle(e_obj)) - np.degrees(np.angle(e_ref)))


def phase_difference_hologram_from_phases_deg(ref_phase_deg, obj_phase_deg):
    """Same as `phase_difference_hologram` but for already-extracted phase-only inputs (Sec. 9.B
    reference field, or any other phase-only source) -- avoids constructing a fake unit-amplitude
    complex number just to re-extract its phase."""
    return wrap_deg(np.asarray(obj_phase_deg, dtype=float) - np.asarray(ref_phase_deg, dtype=float))
