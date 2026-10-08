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


def strip_layers(specs, period_m, height_m, plate_m, eps_r, tan_delta, through_hole, rotation_deg=22.5):
    K = len(specs)
    half = period_m / 2
    xs = [(c - (K - 1) / 2) * period_m for c in range(K)]
    upper, plate = [], []
    for c, sp in enumerate(specs):
        x = xs[c]
        octv = octagon_vertices(sp.radius_m, rotation_deg) + np.array([x, 0.0]) if sp.radius_m > 0 else None
        octv_t = tuple(map(tuple, octv)) if octv is not None else None
        if sp.kind == "pillar":
            if octv_t is not None:
                upper.append(Region(octv_t, eps_r, tan_delta, 1.0, 0.0))
            continue
        # hole cell: nylon square (interfaces only towards non-hole neighbours) + air pocket
        edges = []
        if specs[(c - 1) % K].kind != "hole":
            edges.append(((x - half, -half), (x - half, half)))
        if specs[(c + 1) % K].kind != "hole":
            edges.append(((x + half, -half), (x + half, half)))
        upper.append(Region(_rect(x - half, x + half, -half, half), eps_r, tan_delta, 1.0, 0.0, tuple(edges)))
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
