"""1 x K strips of mixed pillar / hole cells (Phase 19 local-periodicity check).

Physical stack, z downward from the feed side:  upper layer (``height_m``) over a plate layer (``plate_m``).
  pillar cell : octagon (nylon) in air          | plate: solid nylon
  hole cell   : nylon square with octagonal pocket | plate: solid nylon (blind) or same octagon hole (through)
The plate layer is uniform nylon over the whole lens; a through-hole cell punches its octagon through it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .fourier import octagon_vertices
from .rcwa import Region, RegionLayer, UniformLayer, solve_stack


@dataclass(frozen=True)
class CellSpec:
    kind: str  # "pillar" | "hole"
    radius_m: float  # 0 = empty cell (pillar) / solid cell (hole)


def _rect(x0, x1, y0, y1):
    return ((x0, y0), (x1, y0), (x1, y1), (x0, y1))  # CCW


def grid_layers(grid, period_m, height_m, plate_m, eps_r, tan_delta, through_hole, rotation_deg=22.5):
    """Layers for a Kx x Ky periodic arrangement ``grid[ix][iy]`` of CellSpec (cell (0,0) at the lower-left).

    Same stack as :func:`strip_layers`; a hole cell gets an interface edge on every side whose (periodic) neighbour is
    not a hole cell.  ``strip_layers`` is the Ky = 1 case (top / bottom neighbours are the cell itself).
    """
    kx, ky = len(grid), len(grid[0])
    half = period_m / 2
    upper, plate = [], []
    for ix in range(kx):
        for iy in range(ky):
            sp = grid[ix][iy]
            x = (ix - (kx - 1) / 2) * period_m
            y = (iy - (ky - 1) / 2) * period_m
            octv = octagon_vertices(sp.radius_m, rotation_deg) + np.array([x, y]) if sp.radius_m > 0 else None
            octv_t = tuple(map(tuple, octv)) if octv is not None else None
            if sp.kind == "pillar":
                if octv_t is not None:
                    upper.append(Region(octv_t, eps_r, tan_delta, 1.0, 0.0))
                continue
            edges = []
            if grid[(ix - 1) % kx][iy].kind != "hole":
                edges.append(((x - half, y - half), (x - half, y + half)))
            if grid[(ix + 1) % kx][iy].kind != "hole":
                edges.append(((x + half, y - half), (x + half, y + half)))
            if grid[ix][(iy - 1) % ky].kind != "hole":
                edges.append(((x - half, y - half), (x + half, y - half)))
            if grid[ix][(iy + 1) % ky].kind != "hole":
                edges.append(((x - half, y + half), (x + half, y + half)))
            upper.append(
                Region(_rect(x - half, x + half, y - half, y + half), eps_r, tan_delta, 1.0, 0.0, tuple(edges))
            )
            if octv_t is not None:
                upper.append(Region(octv_t, 1.0, 0.0, eps_r, tan_delta))
                if through_hole:
                    plate.append(Region(octv_t, 1.0, 0.0, eps_r, tan_delta))
    layers = [RegionLayer(height_m, 1.0, 0.0, tuple(upper))]
    if plate_m > 0:
        layers.append(
            RegionLayer(plate_m, eps_r, tan_delta, tuple(plate)) if plate else UniformLayer(plate_m, eps_r, tan_delta)
        )
    return layers


def strip_layers(specs, period_m, height_m, plate_m, eps_r, tan_delta, through_hole, rotation_deg=22.5):
    """1 x K strip along x (see :func:`grid_layers`)."""
    return grid_layers(
        [[sp] for sp in specs], period_m, height_m, plate_m, eps_r, tan_delta, through_hole, rotation_deg
    )


def solve_grid(
    grid,
    period_m,
    height_m,
    plate_m,
    eps_r,
    tan_delta,
    through_hole,
    freq_hz,
    theta_deg,
    phi_deg,
    n_harm,
    orders=((0, 0),),
):
    kx, ky = len(grid), len(grid[0])
    layers = grid_layers(grid, period_m, height_m, plate_m, eps_r, tan_delta, through_hole)
    return solve_stack(
        layers, kx * period_m, ky * period_m, freq_hz, theta_deg, phi_deg, n_harm, "normal", orders=orders
    )


def solve_strip(
    specs,
    period_m,
    height_m,
    plate_m,
    eps_r,
    tan_delta,
    through_hole,
    freq_hz,
    theta_deg,
    phi_deg,
    n_harm,
    orders=((0, 0),),
):
    K = len(specs)
    layers = strip_layers(specs, period_m, height_m, plate_m, eps_r, tan_delta, through_hole)
    return solve_stack(layers, K * period_m, period_m, freq_hz, theta_deg, phi_deg, n_harm, "normal", orders=orders)


def local_periodicity_prediction(cell_fields, orders_m):
    """Order-m prediction of a strip from per-cell exit fields (uniform field over each cell's width).

    ``cell_fields`` : complex E_x (or any scalar) of each cell, in order along x.  For a field that is
    piecewise constant over K equal cells the order-m amplitude is
    (1/K) sum_c E_c sinc(m/K) exp(-j 2 pi m (c - (K-1)/2)/K)   [exp(+jwt) convention, order e^{-j...}].
    """
    K = len(cell_fields)
    out = {}
    for m in orders_m:
        acc = 0j
        for c, e in enumerate(cell_fields):
            acc += e * np.exp(1j * 2 * math.pi * m * (c - (K - 1) / 2) / K)
        out[m] = acc / K * (np.sinc(m / K) if m else 1.0)
    return out
