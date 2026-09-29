"""Sec. 28-29: does a unit-cell-level surrogate error actually matter for the beam it is meant to
help synthesize? Builds ONE environment-aware aperture (the same 7x7 gradient grid
`generalization.py` already solves once), reads out each interior sub-cell's local dipole moment
as its aperture-field contribution (physically: the induced-current source term responsible for
re-radiation) for both the true heterogeneous-MoM values and the surrogate's predictions, and
feeds both into the SAME unmodified planar-aperture far-field integral Phase 8E already validated
(`phase8e.beam_reconstruction.far_field_pattern`).
"""

import numpy as np

from ..phase8b.object_wave import Beam
from ..phase8e.beam_reconstruction import (
    compare_continuous_vs_discrete,
    evaluate_beam_reconstruction,
    far_field_pattern,
)
from . import generalization as gen

DEFAULT_INTENDED_BEAM = Beam(theta_deg=20.0, phi_deg=0.0)
THETA_GRID_DEG, PHI_GRID_DEG = np.meshgrid(np.linspace(0.0, 70.0, 36), np.linspace(0.0, 359.0, 72))


def build_aperture_from_cases(cases, period_m, use_surrogate):
    """Sec. 28: (x, y) positions of each interior sub-cell (physical supercell coordinates) and
    its complex aperture value -- either the true MoM dipole moment or the surrogate's prediction
    (whichever is requested); samples where the surrogate returned ``OUT_OF_VALIDATED_DOMAIN`` are
    reported and excluded (never silently zeroed -- Sec. 36)."""
    n = gen.GRID_N
    xs, ys, values, excluded = [], [], [], 0
    for c in cases:
        cx = -n * period_m / 2 + (c.j + 0.5) * period_m
        cy = -n * period_m / 2 + (c.i + 0.5) * period_m
        value = c.px_surrogate if use_surrogate else c.px_mom
        if not isinstance(value, complex):
            excluded += 1
            continue
        xs.append(cx)
        ys.append(cy)
        values.append(value)
    return np.array(xs), np.array(ys), np.array(values, dtype=complex), excluded


def run_beam_validation(cases, period_m, frequency_hz, intended_beams=(DEFAULT_INTENDED_BEAM,)):
    """``cases``: the interior `generalization.SpatialCase` list from an ALREADY-solved
    environment-aware aperture (Sec. 28's own aperture, typically the gradient test's -- reused
    here, not re-solved, since it is the same one self-consistent MoM solution either way)."""
    x_mom, y_mom, a_mom, _ = build_aperture_from_cases(cases, period_m, use_surrogate=False)
    x_sur, y_sur, a_sur, n_excluded = build_aperture_from_cases(cases, period_m, use_surrogate=True)

    ff_mom = far_field_pattern(x_mom, y_mom, a_mom, frequency_hz, THETA_GRID_DEG, PHI_GRID_DEG, period_m, period_m)
    if len(a_sur) >= 4:
        ff_sur = far_field_pattern(x_sur, y_sur, a_sur, frequency_hz, THETA_GRID_DEG, PHI_GRID_DEG, period_m, period_m)
    else:
        ff_sur = None

    metrics_mom = evaluate_beam_reconstruction(ff_mom, THETA_GRID_DEG, PHI_GRID_DEG, intended_beams)
    metrics_sur = (
        evaluate_beam_reconstruction(ff_sur, THETA_GRID_DEG, PHI_GRID_DEG, intended_beams)
        if ff_sur is not None
        else None
    )
    comparison = compare_continuous_vs_discrete(ff_mom, ff_sur) if ff_sur is not None else None

    return {
        "cases": cases,
        "n_surrogate_excluded": n_excluded,
        "far_field_mom": ff_mom,
        "far_field_surrogate": ff_sur,
        "metrics_mom": metrics_mom,
        "metrics_surrogate": metrics_sur,
        "comparison": comparison,
        "theta_grid_deg": THETA_GRID_DEG,
        "phi_grid_deg": PHI_GRID_DEG,
    }
