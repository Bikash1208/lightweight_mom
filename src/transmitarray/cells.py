"""Per-cell transmission for an x-polarised horn: 2x2 Jones response at the cell's own (theta, phi),
and a direct (no table interpolation) solve of the octagon radius that gives a requested phase."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from src.dielectric_cell import OctagonPillarCell, solve_cell


def basis_matrix(theta_deg: float, phi_deg: float) -> np.ndarray:
    """Columns = transverse (x, y) components of the unit s and p E-vectors:  [[-sin f, cos t cos f],[cos f, cos t sin f]]."""
    t, f = math.radians(theta_deg), math.radians(phi_deg)
    return np.array([[-math.sin(f), math.cos(t) * math.cos(f)], [math.cos(f), math.cos(t) * math.sin(f)]])


def jones_xy(t_sp_matrix: np.ndarray, theta_deg: float, phi_deg: float) -> np.ndarray:
    """Transmission Jones matrix in the (x, y) transverse-E basis:  J = M T M^-1  (T is [out, in], s/p)."""
    m = basis_matrix(theta_deg, phi_deg)
    return m @ t_sp_matrix @ np.linalg.inv(m)


def wrap_deg(a):
    return (np.asarray(a) + 180.0) % 360.0 - 180.0


@dataclass(frozen=True)
class CellFamilyModel:
    """One lens-cell family at fixed stack parameters; ``solve`` returns the s/p transmission matrix."""

    period_m: float
    total_thickness_m: float
    plate_m: float
    eps_r: float
    tan_delta: float
    through_hole: bool
    freq_hz: float
    n_harm: int

    def cell(self, kind: str, radius_m: float) -> OctagonPillarCell:
        if kind == "hole" and self.through_hole:
            return OctagonPillarCell(
                self.period_m, radius_m, self.total_thickness_m, self.eps_r, self.tan_delta, 0.0, 22.5, "hole"
            )
        return OctagonPillarCell(
            self.period_m,
            radius_m,
            self.total_thickness_m - self.plate_m,
            self.eps_r,
            self.tan_delta,
            self.plate_m,
            22.5,
            kind,
        )

    def t_matrix(self, kind: str, radius_m: float, theta_deg: float, phi_deg: float) -> np.ndarray:
        return solve_cell(self.cell(kind, radius_m), self.freq_hz, theta_deg, phi_deg, self.n_harm).t

    def jones(self, kind, radius_m, theta_deg, phi_deg):
        return jones_xy(self.t_matrix(kind, radius_m, theta_deg, phi_deg), theta_deg, phi_deg)


def solve_radius_for_phase(
    model: CellFamilyModel, kind, target_deg, theta_deg, phi_deg, r0_m, r_lo_m, r_hi_m, tol_deg=0.05, max_iter=12
):
    """Secant solve of  arg J_xx(R) = target  on [r_lo, r_hi] starting from r0.

    Returns (R, residual_deg, J_xx, n_solves, J); residual is wrap(arg J_xx - target) at the returned R.
    The returned R is the best (smallest |residual|) iterate; the bounds are never exceeded.
    """

    def evaluate(r):
        j = model.jones(kind, r, theta_deg, phi_deg)
        return wrap_deg(math.degrees(np.angle(j[0, 0])) - target_deg), j

    n = 0
    r = float(np.clip(r0_m, r_lo_m, r_hi_m))
    e, j = evaluate(r)
    n += 1
    best = (abs(e), r, e, j)
    step = 2e-5  # finite-difference step [m]
    for _ in range(max_iter):
        if abs(e) < tol_deg:
            break
        r2 = r + step if r + step <= r_hi_m else r - step
        e2, j2 = evaluate(r2)
        n += 1
        if abs(e2) < best[0]:
            best = (abs(e2), r2, e2, j2)
        slope = wrap_deg(e2 - e) / (r2 - r)  # deg per metre
        if abs(slope) < 1e3:  # flat: cannot move usefully
            break
        r_new = float(np.clip(r - e / slope, max(r_lo_m, r - 3e-4), min(r_hi_m, r + 3e-4)))
        if abs(r_new - r) < 1e-9:
            break
        r = r_new
        e, j = evaluate(r)
        n += 1
        if abs(e) < best[0]:
            best = (abs(e), r, e, j)
    _a, r_b, e_b, j_b = best
    return r_b, float(e_b), complex(j_b[0, 0]), n, j_b
