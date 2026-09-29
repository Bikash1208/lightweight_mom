"""Task A: the per-cell dipole-moment deviation Phase 9/12C measured is a proxy (Why Sec. 2 of
the brief) -- what the design actually needs is the reflection-phase/amplitude error of the centre
cell. This module defines that conversion and validates it.

Conversion. ``sparameters.scattering_from_fields`` already accepts an arbitrary *supplied*
mean-current array (it does not require that array to have come from actually solving that exact
cell -- Sec. "mean_currents" of its own docstring just says "the Floquet coefficient of the solved
current for each excitation"). So: given any local (0,0) dipole moment/mean-current for a cell of
the design's own period, plugging it into that same, UNMODIFIED function together with the design's
own ``LayeredGreenFunction``/``ReferencePlanes`` yields the "equivalent reflection coefficient" that
cell would present if its own environment were exactly periodic -- i.e. exactly what
``compute_s_parameters`` itself would report for an ISOLATED solve of that same current. This is
not a new physical model, just the existing S-parameter definition applied to a current that did
not come from an isolated solve.

Validation (required by the brief, "must reproduce S11 exactly"): for a uniform array, the
embedded centre-cell dipole moment IS the isolated solve's own dipole moment (Phase 9/12C's
homogeneous correctness proof), so this conversion applied to it must equal (and does, by
construction: it is the same function call with the same input) the true isolated S11 to machine
precision. See ``tests/test_phase12d_reflection_metric.py``.
"""

from dataclasses import dataclass

import numpy as np

from ..green_function import LayeredGreenFunction
from ..sparameters import scattering_from_fields


@dataclass(frozen=True)
class EquivalentReflection:
    s11_xx: complex
    magnitude: float
    phase_deg: float


def equivalent_reflection_xx(mean_current_x, substrate, frequency_hz, reference_planes, backing="pec"):
    """The S11_xx an isolated periodic cell with this exact (0,0) x-polarized current Floquet
    coefficient (``mean_current_x`` [A/m], e.g. ``basis_functions.RooftopBasis.
    current_coefficients(...)[0]``, NOT a physical dipole moment [A.m] -- see
    ``equivalent_reflection_from_dipole_moment`` for that conversion) would present -- Sec. above.
    Only the x-polarization is used (this design's own excitation, ``design_point.py``:
    ``excitation.polarization: x``); y and cross-pol are set to zero, matching a square patch's own
    S11_yx ~ 0 (Phase 9's own numeric mirror-symmetry finding). Exact (machine precision) when
    ``mean_current_x`` is the true Floquet coefficient -- verified directly against
    ``compute_s_parameters`` (``tests/test_phase12d_reflection_metric.py``)."""
    green = LayeredGreenFunction(substrate, frequency_hz, backing=backing)
    mean_currents = np.zeros((1, 2, 2), dtype=complex)
    mean_currents[0, 0, 0] = mean_current_x   # port 0 (only port, grounded), pol x, component x
    sparams = scattering_from_fields(green, reference_planes, mean_currents)
    s11xx = sparams.component("S11_xx")
    return EquivalentReflection(s11xx, abs(s11xx), float(np.degrees(np.angle(s11xx))))


def equivalent_reflection_from_dipole_moment(dipole_px, period_m, substrate, frequency_hz,
                                              reference_planes, backing="pec"):
    """As ``equivalent_reflection_xx``, but from a physical dipole moment [A.m] over the FULL cell
    footprint (Phase 9's ``dipole_extraction.local_dipole_moment`` with ``half_width = period/2``,
    the quantity Phase 9/12C's own gradient-case results (``isolated_px``/``embedded_px``) store):
    ``mean_current_x = dipole_px / (period_m ** 2)`` (a dipole moment is the Floquet coefficient
    times the cell area it was integrated over). This introduces a small (~0.1-0.2 deg phase,
    <0.1% magnitude, measured) residual relative to the exact conversion above, because
    ``local_dipole_moment`` is itself a finite (48x48) Riemann-sum approximation of the same
    integral ``current_coefficients`` evaluates in closed form -- negligible next to this phase's
    5-10 deg / 0.5-1 dB tolerances, and reported as a known artifact, not hidden."""
    return equivalent_reflection_xx(dipole_px / (period_m ** 2), substrate, frequency_hz,
                                     reference_planes, backing)


def reflection_phase_amplitude_error(isolated_px, embedded_px, period_m, substrate, frequency_hz,
                                      reference_planes, backing="pec"):
    """(phase_error_deg, magnitude_error_db) between the isolated and embedded centre-cell
    equivalent reflection coefficients -- the design-relevant error Task A reports alongside
    Phase 9's raw dipole-moment metric. Inputs are dipole moments [A.m] (Phase 9/12C convention),
    converted via ``equivalent_reflection_from_dipole_moment``."""
    from ..phase import wrapped_difference_deg

    iso = equivalent_reflection_from_dipole_moment(isolated_px, period_m, substrate, frequency_hz,
                                                    reference_planes, backing)
    emb = equivalent_reflection_from_dipole_moment(embedded_px, period_m, substrate, frequency_hz,
                                                    reference_planes, backing)
    phase_error = abs(wrapped_difference_deg(emb.phase_deg, iso.phase_deg))
    mag_db_error = abs(20.0 * np.log10(emb.magnitude / iso.magnitude))
    return phase_error, mag_db_error, iso, emb
