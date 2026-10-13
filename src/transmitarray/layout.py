"""Lens outline and cell grid: regular octagon (across flats D) on an even square grid of pitch p."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LensLayout:
    d_flats_m: float
    pitch_m: float
    focal_m: float

    def __post_init__(self):
        n = self.d_flats_m / self.pitch_m
        if abs(n - round(n)) > 1e-9 or round(n) % 2:
            raise ValueError("D / pitch must be an even integer so the octagon flats sit on cell borders")

    @property
    def n_grid(self) -> int:
        return round(self.d_flats_m / self.pitch_m)

    def grid_centers(self):
        """Full square grid of cell centres, arrays (n, n) indexed [ix, iy]."""
        n, p = self.n_grid, self.pitch_m
        c = (np.arange(n) - (n - 1) / 2) * p
        return np.meshgrid(c, c, indexing="ij")

    def inside(self, x, y):
        """Cell-centre test against the regular octagon (apothem D/2): |x|,|y| <= a and |x|+|y| <= a*sqrt(2)."""
        a = self.d_flats_m / 2
        return (np.abs(x) <= a) & (np.abs(y) <= a) & (np.abs(x) + np.abs(y) <= a * math.sqrt(2) + 1e-12)

    def cells(self):
        """(ix, iy, x, y) of the cells whose centre lies inside the outline, row-major."""
        X, Y = self.grid_centers()
        ix, iy = np.nonzero(self.inside(X, Y))
        return ix, iy, X[ix, iy], Y[ix, iy]

    def octagon_vertices(self):
        a = self.d_flats_m / 2
        r = a / math.cos(math.radians(22.5))
        ang = np.deg2rad(22.5 + 45.0 * np.arange(8))
        return np.stack([r * np.cos(ang), r * np.sin(ang)], axis=1)

    def incidence(self, x, y):
        """(theta_deg, phi_deg, path length m) of the ray from the horn phase centre to cell (x, y)."""
        rho = np.hypot(x, y)
        return (
            np.degrees(np.arctan2(rho, self.focal_m)),
            np.degrees(np.arctan2(y, x)),
            np.sqrt(rho**2 + self.focal_m**2),
        )
