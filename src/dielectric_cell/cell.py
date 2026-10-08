"""Octagonal dielectric pillar unit cell (3D-printed lens element) on top of the RCWA engine.

Geometry: square lattice of pitch ``period_m``; one regular octagon of circumradius ``circumradius_m``
(0 = no pillar, air gap) extruded through ``height_m``; optional uniform base plate of
``base_plate_m`` below it (exit side).  Nylon is lossy: eps = eps_r (1 - j tan_delta), exp(+jwt).

Reference planes: S11 at the top face of the pillar layer (z = 0, incident side) and S21 at the bottom
face of the last layer (z = -(height + base plate)); the incident phase is referenced to the top face and
the transmitted phase to the bottom face, with no air propagation de-embedded.  An empty cell
(circumradius 0, no plate) is an air layer of the full height and so has S11 = 0 and
S21 = exp(-j k0 cos(theta) h), the free-space path.  Lens design uses the *same* total
thickness for every cell, so the absolute phase offset is a common constant.

Polarisation basis: s = E perpendicular to the plane of incidence, p = in the plane.  At normal incidence
phi = 0 deg, s = y and p = x.  Transmission uses identical unit vectors in and out; reflection uses
the tangential-E sign convention (r_pp = r_ss at normal incidence).  Matrices are [out, in] with
index 0 = s, 1 = p.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .fourier import octagon_vertices
from .rcwa import PolygonLayer, UniformLayer, solve_stack


@dataclass(frozen=True)
class OctagonPillarCell:
    period_m: float
    circumradius_m: float
    height_m: float
    eps_r: float
    tan_delta: float
    base_plate_m: float
    rotation_deg: float  # 22.5 = flats on the lattice axes, 0 = vertices on the axes

    def __post_init__(self):
        if not (self.period_m > 0 and self.height_m > 0 and self.eps_r >= 1 and self.tan_delta >= 0):
            raise ValueError("period, height must be > 0, eps_r >= 1, tan_delta >= 0")
        if self.circumradius_m < 0 or self.base_plate_m < 0:
            raise ValueError("circumradius and base plate thickness must be >= 0")
        if self.max_extent_m > self.period_m / 2 * (1 + 1e-9):
            raise ValueError(
                f"octagon (extent {self.max_extent_m * 1e3:.3f} mm) does not fit in a {self.period_m * 1e3:.3f} mm cell"
            )

    @property
    def max_extent_m(self) -> float:
        """Largest |x| or |y| reached by the octagon."""
        a = np.deg2rad(self.rotation_deg + 45.0 * np.arange(8))
        return float(self.circumradius_m * max(np.abs(np.cos(a)).max(), np.abs(np.sin(a)).max()))

    @property
    def total_thickness_m(self) -> float:
        return self.height_m + self.base_plate_m

    @property
    def fill_fraction(self) -> float:
        """Pillar area / cell area  (octagon area = 2*sqrt(2)*R^2)."""
        return 2.0 * math.sqrt(2.0) * self.circumradius_m**2 / self.period_m**2

    def layers(self):
        out = []
        if self.circumradius_m > 0:
            v = tuple(map(tuple, octagon_vertices(self.circumradius_m, self.rotation_deg)))
            out.append(PolygonLayer(self.height_m, v, self.eps_r, self.tan_delta))
        else:
            out.append(UniformLayer(self.height_m, 1.0, 0.0))
        if self.base_plate_m > 0:
            out.append(UniformLayer(self.base_plate_m, self.eps_r, self.tan_delta))
        return out


@dataclass(frozen=True)
class CellResponse:
    r: np.ndarray  # 2x2 [out, in], (s, p)
    t: np.ndarray
    theta_deg: float
    phi_deg: float
    freq_hz: float

    def t_pol(self, pol: str) -> complex:
        return complex(self.t[0, 0] if pol == "TE" else self.t[1, 1])

    def loss_db(self, pol: str) -> float:
        """Transmission loss 20 log10 |S21| (negative number = loss) for co-polar ``pol`` (TE or TM)."""
        return 20.0 * math.log10(abs(self.t_pol(pol)))

    def power_balance(self, pol: str) -> float:
        """|r|^2 + |t|^2 summed over both output polarisations for a unit ``pol`` input (<= 1; 1 if lossless)."""
        j = 0 if pol == "TE" else 1
        return float(np.sum(np.abs(self.r[:, j]) ** 2) + np.sum(np.abs(self.t[:, j]) ** 2))


def solve_cell(
    cell: OctagonPillarCell,
    freq_hz: float,
    theta_deg: float,
    phi_deg: float,
    n_harm: int,
    factorization: str = "normal",
) -> CellResponse:
    res = solve_stack(cell.layers(), cell.period_m, cell.period_m, freq_hz, theta_deg, phi_deg, n_harm, factorization)
    return CellResponse(res["r"], res["t"], theta_deg, phi_deg, freq_hz)
