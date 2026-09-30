"""Phase 8B Sec. 11-18: the object-wave (desired-beam) field, derived entirely from
`excitation.PlaneWaveExcitation` (Phase 1, unmodified) and `CONVENTIONS.md`'s own general plane-
wave rule -- no new sign or coordinate convention is introduced here.

## Derivation (Sec. 12-13's required "show the exact equation, derive the sign")

`excitation.PlaneWaveExcitation` already implements this project's one, tested transverse-
wavenumber convention for a wave characterized by (theta, phi):

    kx0 = k0 sin(theta) cos(phi),   ky0 = k0 sin(theta) sin(phi)          (`.transverse_wavenumber`)

and the project's one plane-wave phase convention, fixed for *any* propagation direction under
`exp(+j*omega*t)` (`CONVENTIONS.md` line 19, `phase8b.conventions.PLANE_WAVE_PHASE_CONVENTION`):

    E(r) = A * exp(-j * k . r)

An *object* (desired-beam) wave radiates away from the surface, i.e. its propagation direction has
a positive z-component (unlike `PlaneWaveExcitation`'s own incident wave, which is defined to
travel toward -z, `wavevector = (kx0, ky0, -kz0)`). The holographic modulation pattern, however, is
recorded on the surface plane z = 0 -- the same plane `CONVENTIONS.md` already fixes for the metal
pattern and the S-parameter phase reference (`phase8b.conventions.COORDINATE_SYSTEM`). At z = 0 the
`kz . z` term of `k . r` vanishes identically, so **the object field's spatial phase pattern on the
surface depends only on the transverse wavenumber (kx0, ky0), which has the identical formula for
an incident or an outgoing wave** -- the kz-sign ambiguity Sec. 13 warns about therefore does not
need to be resolved to compute a hologram at z = 0, and this module says so explicitly rather than
silently picking a kz sign it does not need.

    E_obj(x, y; theta, phi) = A * exp(j*phi0) * exp(-j*(kx0*x + ky0*y))     [at z = 0]

`transverse_wavenumber_deg` below computes kx0, ky0 with this identical formula, cross-checked
against `PlaneWaveExcitation.transverse_wavenumber` by test wherever that class is constructible
(see its own docstring for why it cannot always be reused by direct construction).
"""

import math
from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber

_POLARIZATION_VECTORS = {"x": np.array([1.0, 0.0]), "y": np.array([0.0, 1.0])}

# PHASE15 Part A: the polarization-compensation label set for `desired_polarization`, distinct from
# the legacy `polarization` field above (that one splits THIS module's own scalar object field into
# an (Ex, Ey) pair, Sec. 18 -- a different physical quantity from `hologram.polarization_compensation`'s
# p_hat_b, which projects onto the APERTURE CURRENT's local radial direction to build phi_ext,b,
# see that module's docstring). "natural" (no compensation, phi_ext=1) is the default so a `Beam`
# built without this field reproduces pre-Phase-15 hologram synthesis bit-identically.
DESIRED_POLARIZATIONS = frozenset({"natural", "lhcp", "rhcp", "hp", "vp", "slant_p45", "slant_m45"})

# `excitation.PlaneWaveExcitation` couples its transverse-wavenumber formula to a
# polarization-perpendicularity check (it raises unless the requested x/y polarization is exactly
# transverse to the propagation direction -- only true at theta=0 or at the two azimuths where
# that polarization axis already lies in the plane of incidence's normal). An object-wave beam
# needs the (kx0, ky0) formula at an *arbitrary* (theta, phi) independent of its own polarization
# label (Sec. 18: "arbitrary polarization" combined with "arbitrary propagation direction" is
# exactly the case `PlaneWaveExcitation` does not support -- `CONVENTIONS.md`: "TE/TM and
# circular-polarization labels are not yet defined"). `transverse_wavenumber_deg` therefore
# reimplements only the two-line formula directly (cross-checked byte-for-byte against
# `PlaneWaveExcitation.transverse_wavenumber` in every case where the latter *is* constructible,
# `tests/test_phase8b_object_wave.py::test_transverse_wavenumber_matches_plane_wave_excitation`)
# rather than silently working around the perpendicularity check.


@dataclass(frozen=True)
class Beam:
    """One desired object-wave beam (Sec. 11, 14). ``weight``: complex field-amplitude
    coefficient (Sec. 16 -- see `combine_beams`'s docstring for why this interpretation, not
    power or a normalized coefficient, was chosen). ``initial_phase_deg``: Sec. 17 -- enters as a
    constant additive phase on this beam's field, nothing else; it has no other role in this
    formulation, so it is not exposed as a free parameter anywhere else."""

    theta_deg: float
    phi_deg: float
    weight: complex = 1.0
    initial_phase_deg: float = 0.0
    polarization: str = "x"
    # PHASE15 Part A: independent of `polarization` above -- see `DESIRED_POLARIZATIONS`'s comment.
    desired_polarization: str = "natural"
    oam_charge: int = 0  # l_b: requested vortex OAM charge on top of `desired_polarization`'s own
    # compensation (0 = "Pencil" beam type in the app's terms; nonzero = "Vortex, charge l_b").

    def __post_init__(self):
        if not (0.0 <= self.theta_deg <= 90.0):
            raise ValueError(f"theta_deg must be in [0, 90], got {self.theta_deg!r}")
        if not (0.0 <= self.phi_deg < 360.0):
            raise ValueError(f"phi_deg must be in [0, 360), got {self.phi_deg!r}")
        if self.polarization not in _POLARIZATION_VECTORS:
            raise ValueError(f"polarization must be one of {sorted(_POLARIZATION_VECTORS)}, got {self.polarization!r}")
        if self.desired_polarization not in DESIRED_POLARIZATIONS:
            raise ValueError(
                f"desired_polarization must be one of {sorted(DESIRED_POLARIZATIONS)}, "
                f"got {self.desired_polarization!r}"
            )
        if not isinstance(self.oam_charge, (int, np.integer)):
            raise ValueError(f"oam_charge must be an integer, got {self.oam_charge!r}")


def transverse_wavenumber_deg(frequency_hz, theta_deg, phi_deg):
    """(kx0, ky0) [rad/m] = (k0 sin(theta) cos(phi), k0 sin(theta) sin(phi)) -- the same formula
    `excitation.PlaneWaveExcitation.transverse_wavenumber` implements (see module-level comment
    for why this cannot always be obtained *by constructing* that class for an object wave), using
    `frequency.free_space_wavenumber` (Phase 1, unmodified) for k0."""
    k0 = free_space_wavenumber(frequency_hz)
    theta, phi = math.radians(theta_deg), math.radians(phi_deg)
    s = k0 * math.sin(theta)
    return s * math.cos(phi), s * math.sin(phi)


def single_beam_object_field(x_m, y_m, beam: Beam, frequency_hz):
    """Sec. 14: `E_obj(x, y)` for one beam, at z = 0 -- see module docstring for the derivation.
    Returns a complex array (same shape as ``x_m``/``y_m``) of the *scalar* field along
    ``beam.polarization``; ``vector_field`` below gives the full (Ex, Ey) if needed."""
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    kx0, ky0 = transverse_wavenumber_deg(frequency_hz, beam.theta_deg, beam.phi_deg)
    phase0 = np.radians(beam.initial_phase_deg)
    return beam.weight * np.exp(1j * phase0) * np.exp(-1j * (kx0 * x_m + ky0 * y_m))


def single_beam_vector_field(x_m, y_m, beam: Beam, frequency_hz):
    """(Ex, Ey) arrays: the scalar field of `single_beam_object_field` times the polarization
    unit vector (Sec. 18) -- never collapsed to a scalar silently."""
    scalar = single_beam_object_field(x_m, y_m, beam, frequency_hz)
    pol = _POLARIZATION_VECTORS[beam.polarization]
    return scalar * pol[0], scalar * pol[1]


def combine_beams(x_m, y_m, beams, frequency_hz):
    """Sec. 15: multi-beam object field, combined at the **complex-field level**
    (`E_obj = sum_n w_n E_obj,n`), never by averaging phase maps or angles (Sec. 15's explicit
    prohibition) -- the field-level sum is what a real aperture radiating N simultaneous beams
    physically produces (superposition of the actual phasors), and is the only combination for
    which "the resulting field's phase" is even a well-defined single-valued quantity.

    ``weight`` (Sec. 16) is interpreted here as a **complex field-amplitude coefficient**
    (|w_n| scales that beam's field amplitude, arg(w_n) is an additional phase on top of
    ``initial_phase_deg``) -- the natural quantity a linear field superposition needs. Power
    weighting or a separately-normalized coefficient are different physical quantities that would
    require a different (undocumented, unauthorized) combination rule, so they are not used here
    without that being made explicit by the caller via a complex ``weight``."""
    if not beams:
        raise ValueError("combine_beams: at least one beam is required")
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    ex = np.zeros(np.broadcast(x_m, y_m).shape, dtype=complex)
    ey = np.zeros_like(ex)
    for beam in beams:
        bx, by = single_beam_vector_field(x_m, y_m, beam, frequency_hz)
        ex = ex + bx
        ey = ey + by
    return ex, ey


def combine_beams_scalar(x_m, y_m, beams, frequency_hz):
    """Scalar-field convenience wrapper of `combine_beams` when every beam shares one
    polarization (raises otherwise, rather than silently projecting onto one component)."""
    pols = {b.polarization for b in beams}
    if len(pols) != 1:
        raise ValueError(f"combine_beams_scalar needs all beams to share one polarization, got {pols}")
    ex, ey = combine_beams(x_m, y_m, beams, frequency_hz)
    return ex if beams[0].polarization == "x" else ey
