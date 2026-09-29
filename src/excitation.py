"""Plane-wave excitation.

Geometry/direction convention (documented in CONVENTIONS.md):

* Incident half-space is z > 0; the wave travels toward -z.
* Incident wavevector k_inc = (kx0, ky0, -kz0), field ~ exp(-j*k_inc.r).
* theta = angle of incidence from the surface normal (+z axis reversed),
  phi = azimuth of the in-plane part of k_inc measured from +x toward +y:

      kx0 = k0*sin(theta)*cos(phi),  ky0 = k0*sin(theta)*sin(phi),
      kz0 = k0*cos(theta).

* Polarization is the direction of the incident E-field: 'x' or 'y'. It must be
  perpendicular to k_inc (checked). TE/TM labels are not defined yet.
* Amplitude: the incident field has unit amplitude, |E_inc| = 1 V/m, *evaluated at
  the pattern plane z = 0* (phase reference). Moving the phase reference to the
  S-parameter reference planes is a later (de-embedding) step.

Also here: ``ReferencePlanes`` (explicit S-parameter plane locations; they do not
enter Z or V) and ``stack_tangential_field`` (field of the *unpatterned* stack on
the pattern plane, the "incident" field the PEC metal must cancel).
"""

import math
from dataclasses import dataclass

import numpy as np

from ._validate import finite_real
from .green_function import LayeredGreenFunction, te_tm_unit_vectors
from .frequency import (
    angular_frequency,
    free_space_wavelength,
    free_space_wavenumber,
    validate_frequency,
)

_POLARIZATIONS = {"x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0)}


@dataclass(frozen=True)
class PlaneWaveExcitation:
    frequency: float     # [Hz]
    theta: float = 0.0   # [rad], 0 <= theta < pi/2
    phi: float = 0.0     # [rad]
    polarization: str = "x"

    def __post_init__(self):
        object.__setattr__(self, "frequency", validate_frequency(self.frequency))
        for name in ("theta", "phi"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if not 0.0 <= self.theta < math.pi / 2:
            raise ValueError(f"theta must satisfy 0 <= theta < pi/2, got {self.theta}")
        if self.polarization not in _POLARIZATIONS:
            raise ValueError(
                f"polarization must be one of {sorted(_POLARIZATIONS)}, got {self.polarization!r}"
            )
        k_hat = np.asarray(self.wavevector) / self.k0
        if abs(np.dot(self.e_field_direction, k_hat)) > 1e-12:
            raise ValueError(
                f"polarization {self.polarization!r} is not perpendicular to the incident "
                f"wavevector for theta={self.theta}, phi={self.phi}"
            )

    @property
    def omega(self):
        return angular_frequency(self.frequency)

    @property
    def k0(self):
        return free_space_wavenumber(self.frequency)

    @property
    def wavelength(self):
        return free_space_wavelength(self.frequency)

    @property
    def transverse_wavenumber(self):
        """(kx0, ky0) [rad/m]."""
        s = self.k0 * math.sin(self.theta)
        return s * math.cos(self.phi), s * math.sin(self.phi)

    @property
    def wavevector(self):
        """Incident propagation vector (kx0, ky0, -kz0) [rad/m]."""
        kx0, ky0 = self.transverse_wavenumber
        return (kx0, ky0, -self.k0 * math.cos(self.theta))

    @property
    def e_field_direction(self):
        """Unit vector of the incident E-field."""
        return np.array(_POLARIZATIONS[self.polarization])


@dataclass(frozen=True)
class ReferencePlanes:
    """Explicit S-parameter reference-plane locations (used by ``sparameters.py``).

    The planes sit at the outer boundaries of the modelled air regions:

        top plane    z_top    = +top_air_thickness             (above the pattern at z = 0)
        bottom plane z_bottom = -(substrate_thickness + bottom_air_thickness)

    They do NOT enter the impedance matrix or excitation vector: the Green's
    function treats both air regions as semi-infinite (or, for a grounded stack, the substrate
    as terminated by the ground). They only fix where reflected/transmitted phases are
    referenced when de-embedding is added.

    ``bottom_air_thickness=None`` (Phase 12A): there is no bottom port (a grounded stack,
    ``UnitCell.backing="pec"``); ``z_bottom`` then raises rather than returning a meaningless
    location. Every existing call passes a real ``bottom_air_thickness`` and is unaffected.
    """

    top_air_thickness: float             # [m] >= 0
    bottom_air_thickness: float | None   # [m] >= 0, or None: no bottom port (grounded stack)
    substrate_thickness: float           # [m] > 0

    def __post_init__(self):
        object.__setattr__(
            self, "top_air_thickness", finite_real("top_air_thickness", self.top_air_thickness)
        )
        object.__setattr__(
            self, "substrate_thickness", finite_real("substrate_thickness", self.substrate_thickness)
        )
        if self.bottom_air_thickness is not None:
            object.__setattr__(
                self, "bottom_air_thickness",
                finite_real("bottom_air_thickness", self.bottom_air_thickness),
            )
            if self.bottom_air_thickness < 0.0:
                raise ValueError("air-region thicknesses must be >= 0")
        if self.top_air_thickness < 0.0:
            raise ValueError("air-region thicknesses must be >= 0")
        if self.substrate_thickness <= 0.0:
            raise ValueError("substrate_thickness must be > 0")

    @property
    def z_top(self):
        return self.top_air_thickness

    @property
    def z_bottom(self):
        if self.bottom_air_thickness is None:
            raise ValueError(
                "no bottom reference plane exists for a grounded stack "
                "(ReferencePlanes.bottom_air_thickness is None)"
            )
        return -(self.substrate_thickness + self.bottom_air_thickness)

    @property
    def grounded(self):
        """True iff this plane set has no bottom port (``bottom_air_thickness is None``)."""
        return self.bottom_air_thickness is None


def stack_tangential_field(green: LayeredGreenFunction, excitation: PlaneWaveExcitation, port="top"):
    """Tangential E (Ex, Ey) at z = 0 of the stack WITHOUT metal, for a unit incident wave.

    port 'top': the wave of ``excitation`` (travelling -z) with unit amplitude at z = 0.
    port 'bottom': the same transverse E, incident from z < -h (travelling +z) with unit
    amplitude at z = -h; the field at z = 0 is then the wave transmitted through the substrate.

    Transmission-line model of each mode (derivation: MOM_FORMULATION.md Sec. 4 and
    S_PARAMETERS.md Sec. 3); for the top port

        V(0) = V+ (1 + Gamma),      1 + Gamma = (1 + r)(1 - r u) / (1 - r^2 u),
        u = exp(-2j b h),  r_TE = (a - b)/(a + b),  r_TM = (B - A)/(B + A),
        A = a, B = b/eps_c,  a, b = kz in air / substrate,

    and for the bottom port V(0) = V+ (1 - r^2) exp(-j b h) / (1 - r^2 u). Evaluated through
    ``LayeredGreenFunction.stack_fields`` at z = 0 (verified against independent field matching).

    Returns a complex array of shape (2,). At normal incidence TE and TM coincide.

    ``port='bottom'`` raises ``ValueError`` if ``green.backing == "pec"`` (Phase 12A): a grounded
    stack has no bottom port (``LayeredGreenFunction.stack_fields``).
    """
    if not np.isclose(green.frequency, excitation.frequency, rtol=1e-12, atol=0.0):
        raise ValueError(
            f"Green function frequency {green.frequency} != excitation frequency "
            f"{excitation.frequency}"
        )
    kx0, ky0 = excitation.transverse_wavenumber
    e_t = excitation.e_field_direction[:2].astype(complex)
    return green.stack_fields(kx0, ky0, e_t, 0.0, port)[0]
