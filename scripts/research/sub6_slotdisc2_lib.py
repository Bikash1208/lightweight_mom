"""SUB6 slot-disc 2-GAP cell (MODEL_ONLY): SlotDiscNotch (+u notch, BAR50 bridge) PLUS a second, OPEN notch g2 through the opposite rim (-u slot end).
The disc is then two DC-isolated halves.  Nothing under src/ and nothing in sub6_slotdisc_lib.py is touched (this module wraps it).

OPTIONS: g2 = 0.12 and 0.36 mm, both at odd n = 97 (1 and 3 cell rows; 2 rows are not symmetric about y = 0 on an odd mesh).
RUN (user 2026-11-08, latest): g2 = 0.20 at n = 59 = ONE 0.197 mm cell (modelled width 0.197 mm); the 0.12/0.36 (n 97) and 0.10 (n 117) entries were not run.
Mesh resolution of g2 (cell-centre sampling of the metal mask, src/basis_functions.RooftopBasis): with ODD n the centre row y = 0 is in the gap, so the gap
is ONE cell wide (= dy = p/n) as long as dy > g2/2.  g2 = 0.10 mm -> n = 117 (dy 0.0991, 1 cell);  g2 = 0.20 mm -> n = 59 (dy 0.1966, 1 cell).
Even n would drop a gap narrower than dy altogether.  (Mesh 96 would give 0.12 mm for BOTH g2 values, so it cannot separate them.)
"""
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_slotdisc_lib as B
from src._validate import finite_real
from src.basis_functions import MeshSettings, RooftopBasis
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

MM, F = B.MM, B.F
state_summary, next_node, interp_y, make_sub = B.state_summary, B.next_node, B.interp_y, B.make_sub
MESH_FOR_G2 = {0.20: 59, 0.12: 97, 0.36: 97, 0.10: 117}  # odd n = 97: dy = 0.1196 mm; g2 0.12 -> 1 mesh row (0.1196 mm), g2 0.36 -> 3 mesh rows (0.3588 mm)
CLASS_STATES = {"OPEN2": ["NOCUT2"], "BRIDGE2": ["OFF_BAR50", "ON_BAR50"]}


@dataclass(frozen=True)
class SlotDiscNotch2(B.SlotDiscNotch):
    notch2_width: float = 0.0  # open notch through the -x rim, |y| <= g2/2, from x = -L/2 to the disc edge

    def __post_init__(self):
        B.SlotDiscNotch.__post_init__(self)
        object.__setattr__(self, "notch2_width", finite_real("notch2_width", self.notch2_width))

    def in_feature(self, x, y):
        inside = B.SlotDiscNotch.in_feature(self, x, y)
        if self.notch2_width > 0:
            x, y = np.asarray(x, float), np.asarray(y, float)
            inside &= ~((x <= -self.slot_length / 2) & (np.abs(y) <= self.notch2_width / 2))
        return inside


def make_cell(p_mm, r_mm, L_mm, W_mm, g2_mm, sub, state):
    ld = "open" if state == "NOCUT2" else B.STATES[state][1]
    geo = SlotDiscNotch2(r_mm * MM, L_mm * MM, W_mm * MM, B.G_C, bridge=isinstance(ld, dict), notch2_width=g2_mm * MM)
    loads = ()
    if isinstance(ld, dict):
        loads = (LumpedLoad(geo.diode_x, 0.0, length=B.G_C, width=B.DIODE_ALONG, direction="y", **ld),)
    return UnitCell(RectangularLattice(p_mm * MM, p_mm * MM), geo, sub, Metal(), backing="pec", loads=loads)


def node_task(a):
    spec, cls, direction, kt, n, forder, f = a
    sub, p, r, L, W, g2 = spec
    cell = make_cell(p, r, L, W, g2, make_sub(sub), CLASS_STATES[cls][0])
    k0 = free_space_wavenumber(f)
    kx, ky = (kt * k0, 0.0) if direction == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(cell.substrate, f, backing="pec")
    basis = RooftopBasis(cell, n, n)
    t0 = time.time()
    z0 = assemble_impedance_matrix(basis, green, kx, ky, forder)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for st in CLASS_STATES[cls]:
        c = make_cell(p, r, L, W, g2, make_sub(sub), st)
        z = z0 + load_matrix(basis, list(c.loads), green.omega) if c.loads else z0
        cur = solve(z, v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[st] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(cls=cls, direction=direction, kt=kt, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))


def gap_cells(p_mm, r_mm, L_mm, W_mm, g2_mm, n):
    """Number of empty mesh cells across the second gap on the line x = -(L/2 + r)/2 (must be >= 1), and the effective width [mm]."""
    cell = make_cell(p_mm, r_mm, L_mm, W_mm, g2_mm, B.make_sub(("rt5880", 2.2, 0.0009, 3.175)), "NOCUT2")
    dy = p_mm / n
    ys = -p_mm / 2 + (np.arange(n) + 0.5) * dy
    xm = -0.5 * (L_mm / 2 + r_mm)
    m = cell.is_metal(np.full(n, xm * MM), ys * MM)
    k = int(np.sum(~np.asarray(m)[np.abs(ys) < 1.5]))
    return k, k * dy
