"""Horn model: field pattern cos^q(theta_f), x-polarised, phase centre on axis (ideal horn, no cross-pol)."""

from __future__ import annotations

import math

import numpy as np


def q_for_edge_taper(edge_theta_deg: float, edge_taper_db: float) -> float:
    """cos^q E-field exponent giving ``edge_taper_db`` (positive dB below peak) at the lens edge angle."""
    return -edge_taper_db / 20.0 * math.log(10.0) / math.log(math.cos(math.radians(edge_theta_deg)))


def field_amplitude(theta_f_deg, r_m, q: float):
    """Relative incident |E| at a cell: cos^q(theta_f)/r (spherical spreading included)."""
    return np.cos(np.radians(theta_f_deg)) ** q / r_m


def incident_power_density_cos(theta_f_deg, r_m, q: float):
    """Relative power flux through the lens plane per unit area: |E|^2 cos(theta_f)  (eta0 factor dropped)."""
    return field_amplitude(theta_f_deg, r_m, q) ** 2 * np.cos(np.radians(theta_f_deg))


def total_feed_power(q: float) -> float:
    """Same normalisation: integral of cos^{2q}(theta) over the forward hemisphere = 2 pi / (2q + 1)."""
    return 2.0 * math.pi / (2.0 * q + 1.0)
