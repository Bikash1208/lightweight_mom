"""PARALLEL TEST (model only): full Bloch sheet tensors of the static NOCUT cells and the island cell at kt = beta along any direction theta' relative to the slot axis.
Same MoM path as the cell scan (sub6_slotdisc_lib.node_task / sub6_island_cell.node_task: RooftopBasis, assemble_impedance_matrix, LumpedLoad, Y = J (I + G00 J)^-1), but kt = beta*k0*(cos t, sin t) in any direction.
MODEL_ONLY / NOT_HFSS_VALIDATED.  Cells are NOT rotated: rotating kt by theta' relative to the (fixed) slot axis is the same cell-frame tensor Y(theta') (lattice orientation differs, flagged).
"""
import math
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_slotdisc_lib as S
from src.basis_functions import RooftopBasis
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

MM = 1e-3
F = S.F
SUB = ("rt5880", 2.2, 0.0009, 3.175)
P = 11.6
FORDER = 128
GAP_CELLS = 5


class SlotDiscIsland:
    """Same as sub6_island_cell.SlotDiscIsland (duck-types Geometry) with the mesh size dx_mm as a parameter."""
    kind = "patch_circle_slot_island"
    polarity = "metal"
    rotation = 0.0

    def __init__(self, radius, n_l, n_w, n_i, cb, dx_mm, island=True):
        self.radius, self.n_l, self.n_w, self.n_i, self.cb, self.island, self.dx = radius, n_l, n_w, n_i, cb, island, dx_mm

    @property
    def dxm(self):
        return self.dx * MM

    @property
    def y_top_island(self):
        return (self.n_w - GAP_CELLS) * self.dxm

    @property
    def y_bot_island(self):
        return (-self.n_w + self.cb) * self.dxm

    @property
    def diode_y(self):
        return 0.5 * (self.y_top_island + self.n_w * self.dxm)

    @property
    def feature_area(self):
        return math.pi * self.radius**2

    def corners(self):
        r = self.radius
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])

    def in_feature(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        d = self.dxm
        inside = np.hypot(x, y) <= self.radius
        slot = (np.abs(x) <= self.n_l * d + 1e-9) & (np.abs(y) <= self.n_w * d + 1e-9)
        inside &= ~slot
        if self.island:
            isl = (np.abs(x) <= self.n_i * d + 1e-9) & (y >= self.y_bot_island - 1e-9) & (y <= self.y_top_island + 1e-9)
            sheet = (np.abs(x) <= S.DIODE_ALONG / 2) & (np.abs(y - self.diode_y) <= S.G_C / 2)
            inside = inside | (slot & (isl | sheet))
        return inside


def island_cell(cells, dx_mm, state, r_mm=5.435):
    n_l, n_w, n_i, cb = cells
    geo = SlotDiscIsland(r_mm * MM, n_l, n_w, n_i, cb, dx_mm)
    loads = ()
    if state != "NOCUT":
        ld = S.STATES[state][1]
        loads = (LumpedLoad(0.0, geo.diode_y, length=S.G_C, width=S.DIODE_ALONG, direction="y", **ld),)
    return UnitCell(RectangularLattice(P * MM, P * MM), geo, S.make_sub(SUB), Metal(), backing="pec", loads=loads)


def static_cell(row):
    r = row["row"]
    return S.make_cell(r["p"], r["r"], r["L"], r["W"], S.make_sub(tuple(r["sub"])), "NOCUT")


def tensor_task(a):
    """a = (kind, spec, mesh_n, theta_deg, kt_over_k0, states) -> dict(state -> 2x2 complex Y in the lattice (x,y) frame), seconds.
    kind 'static': spec = rows_RT5880 row dict;  kind 'island': spec = dict(cells=[n_l,n_w,n_i,cb])."""
    kind, spec, n, th, kt, states = a
    t0 = time.time()
    k0 = free_space_wavenumber(F)
    kx, ky = kt * k0 * math.cos(math.radians(th)), kt * k0 * math.sin(math.radians(th))
    mk = (lambda st: island_cell(spec["cells"], P / n, st)) if kind == "island" else (lambda st: static_cell(spec))
    c0 = mk(states[0])
    green = LayeredGreenFunction(c0.substrate, F, backing="pec")
    basis = RooftopBasis(c0, n, n)
    z0 = assemble_impedance_matrix(basis, green, kx, ky, FORDER)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for st in states:
        c = mk(st)
        z = z0 + load_matrix(basis, list(c.loads), green.omega) if c.loads else z0
        cur = solve(z, v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[st] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(kind=kind, theta=th, kt=kt, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))
