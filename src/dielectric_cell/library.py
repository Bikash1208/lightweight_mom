"""28 GHz nylon octagon cell library (Phase 19 Part B): geometry grid, sweep, and 360-degree window selection.

Two cell families on the same pitch / total thickness (so the absolute phase offset is common):
  ``pillar``  nylon octagon in air, standing on a base plate;
  ``hole``    solid nylon with an octagonal air pocket (open on the feed side, closed by the plate).
A pillar of circumradius R and a hole of circumradius R are *different* cells; together they cover
nylon fill 0..1.  The minimum printable feature limits R to [R_min, R_max] in both families.
"""

from __future__ import annotations

import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from .cell import OctagonPillarCell, solve_cell

COS_FLAT = math.cos(math.radians(22.5))


@dataclass(frozen=True)
class LibraryParams:
    freq_hz: float
    period_m: float
    total_thickness_m: float
    base_plate_m: float
    eps_r: float
    tan_delta: float
    min_feature_m: float
    n_per_family: int
    thetas_deg: tuple
    phis_deg: tuple
    n_harm: int
    through_hole: bool = False  # hole-family pocket also cuts the base plate (continuous nylon lattice)

    @property
    def r_min_m(self) -> float:
        return self.min_feature_m / (2 * COS_FLAT)

    @property
    def r_max_m(self) -> float:
        return (self.period_m - self.min_feature_m) / (2 * COS_FLAT)


def cell_list(p: LibraryParams) -> list[tuple[str, float]]:
    """[(kind, circumradius_m)] : empty cell, pillars R_min..R_max, holes R_min..R_max, solid cell."""
    radii = np.linspace(p.r_min_m, p.r_max_m, p.n_per_family)
    return (
        [("pillar", 0.0)]
        + [("pillar", float(r)) for r in radii]
        + [("hole", float(r)) for r in radii]
        + [("hole", 0.0)]
    )


def make_cell(p: LibraryParams, kind: str, radius_m: float) -> OctagonPillarCell:
    if kind == "hole" and p.through_hole:
        return OctagonPillarCell(p.period_m, radius_m, p.total_thickness_m, p.eps_r, p.tan_delta, 0.0, 22.5, "hole")
    return OctagonPillarCell(
        p.period_m, radius_m, p.total_thickness_m - p.base_plate_m, p.eps_r, p.tan_delta, p.base_plate_m, 22.5, kind
    )


def _solve_one(args):
    p, kind, r = args
    c = make_cell(p, kind, r)
    rows = []
    for th in p.thetas_deg:
        for ph in p.phis_deg if th > 0 else (0.0,):
            resp = solve_cell(c, p.freq_hz, th, ph, p.n_harm)
            rows.append((kind, r, c.fill_fraction, th, ph, resp.t.copy(), resp.r.copy()))
    return rows


def sweep(p: LibraryParams, workers: int = 4):
    """Run the whole library; returns a list of rows (kind, R, fill, theta, phi, t(2x2), r(2x2))."""
    jobs = [(p, k, r) for k, r in cell_list(p)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        out = list(ex.map(_solve_one, jobs))
    return [row for chunk in out for row in chunk]


def normal_incidence_table(rows):
    """Per cell at theta=0: (kind, R, fill, phase_deg_wrapped, loss_db_of_|t_pp|)."""
    tab = []
    for kind, r, fill, th, _ph, t, _r in rows:
        if th == 0:
            tab.append((kind, r, fill, float(np.degrees(np.angle(t[1, 1]))), 20 * math.log10(abs(t[1, 1]))))
    return tab


def order_by_fill(tab):
    """Both families sorted by nylon fill, with the transmission phase unwrapped along the fill axis."""
    tab = sorted(tab, key=lambda x: x[2])
    ph = np.degrees(np.unwrap(np.radians([x[3] for x in tab])))
    return [(*x[:3], float(u), x[4]) for x, u in zip(tab, ph, strict=True)]


def best_window(ordered, span_deg: float = 360.0):
    """Contiguous (in fill) window of >= span_deg phase minimising the worst-case loss.

    Returns (i0, i1, worst_loss_db, span_deg).  Cells are ordered by fill; phase is monotone in fill.
    """
    best = None
    n = len(ordered)
    for i in range(n):
        for j in range(i, n):
            span = abs(ordered[j][3] - ordered[i][3])
            if span >= span_deg:
                worst = min(x[4] for x in ordered[i : j + 1])
                if best is None or worst > best[2]:
                    best = (i, j, worst, span)
                break
    return best
