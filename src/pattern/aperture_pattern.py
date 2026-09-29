"""Far-field pattern, beam pointing, side-lobe level (SLL) and a directivity estimate from an
aperture field. Built entirely on `phase8e.beam_reconstruction.far_field_pattern` and
`detect_local_maxima` (unmodified) -- this module adds the summary statistics (pointing/SLL/
directivity) an app page needs, not new radiation-integral physics.
"""

from dataclasses import dataclass

import numpy as np

from ..phase8e.beam_reconstruction import detect_local_maxima, far_field_pattern


@dataclass(frozen=True)
class PatternResult:
    theta_deg: np.ndarray       # 2-D angular grid
    phi_deg: np.ndarray
    power_db: np.ndarray         # relative to the pattern's own peak, 0 dB at peak
    main_beam_theta_deg: float
    main_beam_phi_deg: float
    sll_db: "float | None"       # highest sidelobe relative to the main beam, None if only 1 peak found
    directivity_estimate_dbi: float


def compute_pattern(x_m, y_m, aperture_field, frequency_hz, dx_m, dy_m,
                     theta_max_deg=90.0, n_theta=91, n_phi=73):
    """Far field over a hemisphere (theta in [0, theta_max_deg], phi in [0, 360)), on an
    `n_theta` x `n_phi` grid. Returns a `PatternResult`."""
    theta = np.linspace(0.0, theta_max_deg, n_theta)
    phi = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")

    field = far_field_pattern(x_m, y_m, aperture_field, frequency_hz, theta_grid, phi_grid, dx_m, dy_m)
    power = np.abs(field) ** 2
    peak = power.max()
    if peak <= 0.0:
        raise ValueError("compute_pattern: aperture field radiates zero power (all-zero input?)")
    power_db = 10.0 * np.log10(power / peak)

    peaks = detect_local_maxima(power, theta_grid, phi_grid, min_relative_db=-40.0)
    peaks = sorted(peaks, key=lambda p: -p.relative_power_db)
    main = peaks[0]
    sll_db = peaks[1].relative_power_db if len(peaks) > 1 else None

    directivity = _directivity_estimate(power, theta_grid)
    return PatternResult(
        theta_deg=theta_grid, phi_deg=phi_grid, power_db=power_db,
        main_beam_theta_deg=main.theta_deg, main_beam_phi_deg=main.phi_deg,
        sll_db=sll_db, directivity_estimate_dbi=directivity,
    )


def _directivity_estimate(power, theta_grid_deg):
    """D = 4*pi*P_max / integral(P * sin(theta) dtheta dphi) -- the standard directivity
    definition, evaluated by the same discrete grid the pattern was computed on (a hemisphere-only
    approximation for a reflectarray/HMSA that radiates into the forward half-space only: exact
    for a pattern with negligible back-radiation, which this project's own aperture-integration
    model cannot represent in the first place -- it has no notion of a backing structure's
    back-lobe)."""
    theta_rad = np.radians(theta_grid_deg)
    dtheta = np.radians(theta_grid_deg[1, 0] - theta_grid_deg[0, 0]) if theta_grid_deg.shape[0] > 1 else 1.0
    dphi = np.radians(360.0 / theta_grid_deg.shape[1])
    solid_angle_weight = np.sin(theta_rad) * dtheta * dphi
    integral = float(np.sum(power * solid_angle_weight))
    if integral <= 0.0:
        return float("nan")
    return 10.0 * np.log10(4.0 * np.pi * power.max() / integral)
