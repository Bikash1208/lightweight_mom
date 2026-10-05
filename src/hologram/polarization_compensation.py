"""Phase 15 Part A: the per-beam polarization-compensation term phi_ext,b (methodology source:
Ghosh, Ghosh, Saha, IEEE OJAP 2026, Table 1 and eqs. 13-15; `docs/briefs/
PHASE15_POLARIZATION_VORTEX_MULTIFEED.md`).

## Why

A monopole-fed HMSA's surface current is radial about the feed, `J = rho_hat(x,y) * J_s(x,y)`,
`rho_hat = x_hat*cos(phi') + y_hat*sin(phi')` (`phi'` = azimuth about the feed). The field a small
patch of that current radiates into a chosen far-field polarization `p_hat` is proportional to the
projection `conj(p_hat) . rho_hat` -- so to make the WHOLE aperture add constructively into `p_hat`
at the beam direction, each cell needs an extra hologram phase that cancels that projection's own
phase before the ordinary object-wave phase is applied:

    phi_ext,b(x,y) = exp(-j*arg(conj(p_hat_b) . rho_hat(x,y))) * exp(+j*l_b*phi'(x,y))

(the brief's eq., l_b = 0 for a "Pencil" beam; the second factor is an independently-requested extra
OAM winding, e.g. for a deliberately-designed vortex beam). This is ONE generic formula -- rather
than branching per polarization label -- because every special case in the brief falls out of it
algebraically: `p_hat` real and along x_hat (Horizontal) gives `conj(p_hat).rho_hat = cos(phi')`,
a real number whose `arg` is 0 or pi depending on its sign, i.e. exactly the brief's "0/pi flip
across the line through the feed perpendicular to x"; `p_hat` along y_hat (Vertical) gives
`sin(phi')`, the brief's "sign of sin(phi')"; `p_hat = (x_hat+y_hat)/sqrt(2)` (Slant +45) gives
`cos(phi'-45deg)`, the brief's "sign of cos(phi'-45deg)"; a circular `p_hat = (x_hat -+ j*y_hat)/sqrt(2)`
gives `conj(p_hat).rho_hat = exp(+-j*phi')/sqrt(2)`, a spiral -- matching the brief's "CP: p_hat*.rho_hat
= e^{+-jphi'}/sqrt(2) -> phi_ext is a spiral e^{-+jphi'}".

## Handedness (derived from this project's own convention, not the paper's)

Which of the two circular basis vectors `polarization.E_R`/`polarization.E_L` is "RHCP" depends on
the wave's propagation direction (`CONVENTIONS.md`'s `exp(+j*omega*t)`, `src/polarization.py`'s own
handedness derivation). An object (desired-beam) wave radiates AWAY from the surface, i.e. toward
`+z` (`phase8b.object_wave`'s own docstring) -- so `p_hat_b` for "RHCP"/"LHCP" is picked by calling
`polarization.handedness(label, "+z")` and choosing whichever of `E_R`/`E_L` matches, never by
hardcoding one of them (the paper's own sign convention may, and need not, agree with this
project's -- `CONVENTIONS.md`: "Not yet fixed... which handedness naming... the holographic
synthesis... will use"; Phase 15 fixes it here, for the object wave, by that same derivation).

## The "keep the pi/4 configurable" note

The paper's own Table 1 writes the CP case as `exp(-+j(atan2(...) + pi/4))` -- an extra fixed phase
on top of the spiral. Multiplying `p_hat` itself by `exp(+j*cp_phase_offset_deg)` is algebraically
equivalent to adding that same constant to `phi_ext,b`'s own phase (conjugating a Jones vector by an
extra phase subtracts that phase from `arg(conj(p_hat).rho_hat)`, so `phi_ext = exp(-j*arg(...))`
picks it back up as a positive additive phase) -- so `cp_phase_offset_deg` is exposed as exactly
that knob, applied only to the CP labels (the brief calls it "a constant... keep it configurable"),
default 0 (no offset; the derivation above needs none of its own tests to add one to pass).
"""

import numpy as np

from ..phase8b.object_wave import single_beam_object_field
from ..polarization import E_L, E_R, handedness

# PHASE15 Part A: real-valued Jones vectors for the linear labels (HP/VP/slant), matching
# `phase8b.object_wave.DESIRED_POLARIZATIONS` minus "natural" (handled separately, below) and minus
# the circular labels (handedness-dependent, resolved by `_circular_jones_vector`).
_SQRT_HALF = 2.0**-0.5
_LINEAR_JONES_VECTORS = {
    "hp": np.array([1.0, 0.0], dtype=complex),
    "vp": np.array([0.0, 1.0], dtype=complex),
    "slant_p45": np.array([_SQRT_HALF, _SQRT_HALF], dtype=complex),
    "slant_m45": np.array([_SQRT_HALF, -_SQRT_HALF], dtype=complex),
}

# Object-wave propagation direction (`phase8b.object_wave`'s own docstring: "radiates away from the
# surface, i.e. its propagation direction has a positive z-component") -- fixed once here, reused
# by every circular-label lookup rather than repeated at each call site.
_OBJECT_WAVE_DIRECTION = "+z"


def _circular_jones_vector(label, cp_phase_offset_deg):
    """"rhcp"/"lhcp" -> whichever of `polarization.E_R`/`E_L` is that handedness for a wave
    travelling `_OBJECT_WAVE_DIRECTION`, per `polarization.handedness` (derived, not hardcoded --
    see module docstring)."""
    wanted = "RHCP" if label == "rhcp" else "LHCP"
    for circ_label, basis_vector in (("R", E_R), ("L", E_L)):
        if handedness(circ_label, _OBJECT_WAVE_DIRECTION) == wanted:
            return basis_vector * np.exp(1j * np.radians(cp_phase_offset_deg))
    raise AssertionError(f"unreachable: neither E_R nor E_L is {wanted} travelling {_OBJECT_WAVE_DIRECTION}")


def polarization_jones_vector(desired_polarization, cp_phase_offset_deg=0.0):
    """The unit Jones vector `p_hat_b` for a `Beam.desired_polarization` label, or `None` for
    "natural" (no compensation -- the caller must special-case this, since there is no p_hat to
    project onto `rho_hat` at all, not a degenerate p_hat)."""
    if desired_polarization == "natural":
        return None
    if desired_polarization in ("rhcp", "lhcp"):
        return _circular_jones_vector(desired_polarization, cp_phase_offset_deg)
    if desired_polarization in _LINEAR_JONES_VECTORS:
        return _LINEAR_JONES_VECTORS[desired_polarization]
    raise ValueError(f"unknown desired_polarization {desired_polarization!r}")


def polarization_compensation_term(x_m, y_m, feed_xy_m, desired_polarization, oam_charge=0, cp_phase_offset_deg=0.0):
    """phi_ext,b(x,y), the module docstring's formula, broadcast over `x_m`/`y_m` (typically the
    hologram's own cell-centre grid). `feed_xy_m` = (x0, y0) [m], the feed this beam's current is
    radial about (Part C: the beam's ASSIGNED feed in a multi-feed design, not necessarily the
    aperture centre). Returns a unit-modulus complex array (a pure phase factor, by construction:
    both factors of the product have |.| = 1 everywhere) the same shape as `x_m`."""
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    x0, y0 = feed_xy_m
    phi_prime = np.arctan2(y_m - y0, x_m - x0)

    p_hat = polarization_jones_vector(desired_polarization, cp_phase_offset_deg)
    if p_hat is None:  # "natural": no projection factor at all, phi_ext = 1 from this factor
        projection_term = np.ones_like(phi_prime, dtype=complex)
    else:
        rho_hat_x, rho_hat_y = np.cos(phi_prime), np.sin(phi_prime)
        projection = np.conj(p_hat[0]) * rho_hat_x + np.conj(p_hat[1]) * rho_hat_y
        projection_term = np.exp(-1j * np.angle(projection))

    vortex_term = np.exp(1j * oam_charge * phi_prime)  # oam_charge=0 -> identically 1
    return projection_term * vortex_term


def _feed_xy_for_beam(feed_xy_m, i):
    """`feed_xy_m` is either one (x0, y0) shared by every beam (single-feed, the default), or a
    sequence of one (x0, y0) per beam (Part C: each beam radiates about its own ASSIGNED feed)."""
    x0_or_seq, _y0_or_seq = feed_xy_m
    if np.isscalar(x0_or_seq):
        return feed_xy_m
    return feed_xy_m[0][i], feed_xy_m[1][i]


def combine_beams_with_polarization(
    x_m, y_m, beams, frequency_hz, feed_xy_m=(0.0, 0.0), cp_phase_offset_deg=0.0, sector_owner_index=None
):
    """Sum_b weight_b * phi_ext,b(x,y) * single_beam_object_field(beam_b) -- the polarization-
    compensated generalization of `phase8b.object_wave.combine_beams_scalar`, reusing
    `single_beam_object_field` (Sec. 14, unmodified) for each beam's own scalar spatial pattern.

    **Bit-identical invariant (tested)**: when every beam has `desired_polarization="natural"` and
    `oam_charge=0`, `phi_ext,b` is identically 1 for every beam (see `polarization_compensation_term`),
    so this sum reduces EXACTLY to `combine_beams_scalar`'s own sum -- same `Beam.weight` semantics,
    no renormalization -- which is what makes pre-Phase-15 designs left at these defaults reproduce
    bit-identically.

    `sector_owner_index` (Part A's "sector multiplexing" option, paper eq. 17): `None` (default) is
    today's full-aperture combination (every beam contributes everywhere). Given an integer array
    the shape of `x_m`, holding for each cell the index (into `beams`) of the one beam that cell is
    assigned to, each cell instead uses ONLY its owning beam's own (single-beam) contribution -- a
    spatial partition of the aperture, not a different per-beam formula.
    """
    if not beams:
        raise ValueError("combine_beams_with_polarization: at least one beam is required")
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    e_obj = np.zeros(np.broadcast(x_m, y_m).shape, dtype=complex)
    for i, beam in enumerate(beams):
        scalar = single_beam_object_field(x_m, y_m, beam, frequency_hz)
        feed_xy_b = _feed_xy_for_beam(feed_xy_m, i)
        phi_ext = polarization_compensation_term(
            x_m, y_m, feed_xy_b, beam.desired_polarization, beam.oam_charge, cp_phase_offset_deg
        )
        contribution = scalar * phi_ext
        if sector_owner_index is not None:
            contribution = np.where(sector_owner_index == i, contribution, 0.0)
        e_obj = e_obj + contribution
    return e_obj


def oam_sign_mismatch(desired_polarization, oam_charge):
    """True when the requested OAM charge has the sign that does NOT give a clean circular vortex
    with a monopole feed: in this code's convention RHCP needs l<0 and LHCP needs l>0 (measured:
    RHCP with l=+2 puts its power into the LHCP pencil component). Such a request is a cross-polar
    pencil rather than the vortex asked for, so the app warns."""
    if oam_charge == 0 or desired_polarization not in ("rhcp", "lhcp"):
        return False
    return (desired_polarization == "rhcp" and oam_charge > 0) or (desired_polarization == "lhcp" and oam_charge < 0)
