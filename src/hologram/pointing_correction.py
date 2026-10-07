"""Optional pointing correction (Phase 16 follow-ups, item 5). A small aperture squints toward broadside:
the wide array factor (about 0.9 lambda / (L cos theta)) and the falling element pattern shift the
peak, and the error grows as 1/L^2. The correction measures each beam's achieved direction on the
computed vector pattern, then re-aims the requested beam by the predicted miss: the new direction is
the requested one plus the miss, in direction-cosine (u, v) space. Off by default."""

from dataclasses import dataclass, replace

import numpy as np

from ..pattern.vector_aperture_pattern import aperture_current_from_scalar, beam_local_peak_fine

MAX_ITERATIONS = 2
FINE_GRID_STEP_DEG = 0.1
_FINE_GRID_RADIUS_DEG = 3.0


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
    """Achieved direction of every requested beam, from a fine local evaluation of the vector far
    field around each beam's own request (`beam_local_peak_fine`, resolution `FINE_GRID_STEP_DEG`
    deg) rather than a shared hemisphere grid -- see that function's docstring for why (it was the
    cause of pointing correction stalling, PHASE16 follow-ups). All beams in one call are
    normalised to the SAME peak so their relative levels stay meaningful."""
    jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=feed_xy_m)
    raw = [
        beam_local_peak_fine(
            design.x_m, design.y_m, jx, jy, frequency_hz, green, period_m, period_m,
            beam.theta_deg, beam.phi_deg, radius_deg=_FINE_GRID_RADIUS_DEG, step_deg=FINE_GRID_STEP_DEG,
        )
        for beam in beams
    ]
    peak_power = max(power for _, _, power in raw)
    aims = []
    for beam, (theta, phi, power) in zip(beams, raw, strict=True):
        level_db = 10.0 * np.log10(power / peak_power) if peak_power > 0.0 else float("-inf")
        aims.append(BeamAim(beam.theta_deg, beam.phi_deg, theta, phi, level_db))
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
    """`synthesize(beams) -> design`, `measure(design, beams) -> [BeamAim]` (requested = target).
    Returns `(design, beams_used, history)`: `history[k]` holds the aims of pass k (`history[0]` is the
    uncorrected design). The returned design is the pass with the smallest worst-case miss, so a
    correction never ships a worse design than the one it started from. Stops early once every miss is
    below 0.05 deg."""
    if not 1 <= iterations <= MAX_ITERATIONS:
        raise ValueError(f"iterations must be 1..{MAX_ITERATIONS}, got {iterations!r}")
    targets = list(beams)
    used = list(beams)
    design = synthesize(used)
    aims = _against_targets(measure(design, used), targets)
    history = [aims]
    best = (max(aim.miss_deg for aim in aims), design, used)
    for _ in range(iterations):
        if max(aim.miss_deg for aim in aims) < 0.05:
            break
        used = corrected_beams(used, aims, targets)
        design = synthesize(used)
        aims = _against_targets(measure(design, used), targets)
        history.append(aims)
        worst = max(aim.miss_deg for aim in aims)
        if worst < best[0]:
            best = (worst, design, used)
    return best[1], best[2], history


BALANCE_PASSES = 3
_BALANCE_TOLERANCE_DB = 0.1


def aim_dict(aim):
    return {
        "requested_theta_deg": aim.requested_theta_deg, "requested_phi_deg": aim.requested_phi_deg,
        "achieved_theta_deg": aim.achieved_theta_deg, "achieved_phi_deg": aim.achieved_phi_deg,
        "miss_deg": aim.miss_deg, "level_db": aim.level_db,
    }


def aim_report(history):
    """Before/after of a pointing correction, with `after` the best pass (fewest worst-case miss)."""
    best = min(range(len(history)), key=lambda k: max(a.miss_deg for a in history[k]))
    return {
        "passes": best,
        "before": [aim_dict(a) for a in history[0]],
        "after": [aim_dict(a) for a in history[best]],
        "measurement_grid_step_deg": FINE_GRID_STEP_DEG,
    }


def requested_levels_db(beams):
    """The requested relative power of each beam, from its own weight, strongest beam = 0 dB
    (a UI "Relative power" of 0 dB for every beam is therefore equal power)."""
    mags = np.array([abs(b.weight) for b in beams], dtype=float)
    return 20.0 * np.log10(mags / mags.max())


def balance_beam_levels(beams, targets_db, synthesize, measure, passes=BALANCE_PASSES):
    """Iterates `w <- w * 10^((target - achieved) / 20)` on each beam's internal weight so the achieved
    levels (relative to the strongest beam) approach `targets_db`. Keeps the pass with the smallest
    worst-case error. Returns `(design, beams_used, report)`; the report holds the relative levels
    before and after, the targets and the errors."""
    targets = np.asarray(targets_db, dtype=float)

    def relative_levels(aims):
        levels = np.array([a.level_db for a in aims], dtype=float)
        return levels - levels.max()

    used = list(beams)
    design = synthesize(used)
    levels = relative_levels(measure(design, used))
    err = float(np.max(np.abs(levels - targets)))
    history = [levels]
    best = (err, design, used, levels)
    best_pass = 0
    for _ in range(passes):
        if err < _BALANCE_TOLERANCE_DB:
            break
        used = [replace(b, weight=b.weight * 10.0 ** ((t - lv) / 20.0)) for b, t, lv in zip(used, targets, levels, strict=True)]
        design = synthesize(used)
        levels = relative_levels(measure(design, used))
        err = float(np.max(np.abs(levels - targets)))
        history.append(levels)
        if err < best[0]:
            best = (err, design, used, levels)
            best_pass = len(history) - 1
    err_best, design_best, used_best, levels_best = best
    report = {
        "passes": best_pass,
        "targets_db": [float(t) for t in targets],
        "before_db": [float(v) for v in history[0]],
        "after_db": [float(v) for v in levels_best],
        "max_error_before_db": float(np.max(np.abs(history[0] - targets))),
        "max_error_after_db": float(err_best),
    }
    return design_best, used_best, report


def refine_design(beams, synthesize, measure, *, pointing, balance, pointing_iterations=2):
    """The optional stages on one design: pointing correction (when `pointing`), then level matching
    (when `balance` and there are 2+ beams). Returns `(design, beams_used, report)`; `report` has a
    "pointing" and/or "balance" entry for each stage that ran."""
    report = {}
    if pointing:
        design, used, history = correct_beam_aim(beams, synthesize, measure, pointing_iterations)
        report["pointing"] = aim_report(history)
    else:
        used = list(beams)
        design = synthesize(used)
    if balance and len(beams) >= 2:
        design, used, balance_report = balance_beam_levels(used, requested_levels_db(beams), synthesize, measure)
        report["balance"] = balance_report
    return design, used, report
