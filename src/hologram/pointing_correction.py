"""Optional pointing correction (Phase 16 follow-ups, item 5). A small aperture squints toward broadside:
the wide array factor (about 0.9 lambda / (L cos theta)) and the falling element pattern shift the
peak, and the error grows as 1/L^2. The correction measures each beam's achieved direction on the
computed vector pattern, then re-aims the requested beam by the predicted miss: the new direction is
the requested one plus the miss, in direction-cosine (u, v) space. Off by default."""

from dataclasses import dataclass, replace

import numpy as np

from ..pattern.vector_aperture_pattern import aperture_current_from_scalar, beam_local_peak, compute_vector_pattern

MAX_ITERATIONS = 2
_GRID_THETA_MAX_DEG = 60.0
_GRID_N_THETA = 241
_GRID_N_PHI = 73


@dataclass(frozen=True)
class BeamAim:
    requested_theta_deg: float
    requested_phi_deg: float
    achieved_theta_deg: float
    achieved_phi_deg: float
    level_db: float

    @property
    def miss_deg(self):
        """Angular distance between the requested and achieved directions, in degrees."""
        t0, p0 = np.radians(self.requested_theta_deg), np.radians(self.requested_phi_deg)
        t1, p1 = np.radians(self.achieved_theta_deg), np.radians(self.achieved_phi_deg)
        cos_angle = np.cos(t0) * np.cos(t1) + np.sin(t0) * np.sin(t1) * np.cos(p0 - p1)
        return float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))


def _uv(theta_deg, phi_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return np.sin(t) * np.cos(p), np.sin(t) * np.sin(p)


def _theta_phi(u, v):
    return float(np.degrees(np.arcsin(np.clip(np.hypot(u, v), 0.0, 1.0)))), float(np.degrees(np.arctan2(v, u)) % 360.0)


def measure_beam_aims(design, beams, green, period_m, frequency_hz, feed_xy_m):
    """Achieved direction of every requested beam on the design's own vector pattern."""
    jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=feed_xy_m)
    pattern = compute_vector_pattern(
        design.x_m, design.y_m, jx, jy, frequency_hz, green, period_m, period_m,
        _GRID_THETA_MAX_DEG, _GRID_N_THETA, _GRID_N_PHI,
    )
    aims = []
    for beam in beams:
        theta, phi, level = beam_local_peak(pattern, beam.theta_deg, beam.phi_deg)
        aims.append(BeamAim(beam.theta_deg, beam.phi_deg, theta, phi, level))
    return aims


def corrected_beams(beams, aims, targets):
    """Each beam re-aimed by its predicted miss, in (u, v): the direction used this pass (`beams`)
    plus the gap between the target and what it achieved, `d_next = d_used + (target - achieved)`."""
    out = []
    for beam, aim, target in zip(beams, aims, targets, strict=True):
        ud, vd = _uv(beam.theta_deg, beam.phi_deg)
        ut, vt = _uv(target.theta_deg, target.phi_deg)
        ua, va = _uv(aim.achieved_theta_deg, aim.achieved_phi_deg)
        theta, phi = _theta_phi(ud + (ut - ua), vd + (vt - va))
        out.append(replace(beam, theta_deg=theta, phi_deg=phi))
    return out


def _against_targets(aims, targets):
    """The same achieved directions, with each requested direction replaced by the TARGET."""
    return [
        replace(aim, requested_theta_deg=t.theta_deg, requested_phi_deg=t.phi_deg)
        for aim, t in zip(aims, targets, strict=True)
    ]


def correct_beam_aim(beams, synthesize, measure, iterations=1):
    """`synthesize(beams) -> design`, `measure(design, beams) -> [BeamAim]`. Returns
    `(design, beams_used, history)` where `history[k]` is the list of aims at pass k (history[0]
    is the uncorrected design). The loop stops early once every miss is below 0.05 deg."""
    if not 1 <= iterations <= MAX_ITERATIONS:
        raise ValueError(f"iterations must be 1..{MAX_ITERATIONS}, got {iterations!r}")
    targets = list(beams)
    used = list(beams)
    design = synthesize(used)
    aims = _against_targets(measure(design, used), targets)
    history = [aims]
    for _ in range(iterations):
        if max(aim.miss_deg for aim in aims) < 0.05:
            break
        used = corrected_beams(used, aims, targets)
        design = synthesize(used)
        aims = _against_targets(measure(design, used), targets)
        history.append(aims)
    return design, used, history
