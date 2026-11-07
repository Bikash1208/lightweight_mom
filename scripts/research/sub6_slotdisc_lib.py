"""SUB6 slot-disc cell scan (MODEL_ONLY): slotted circular patch (Phase 26 geometry) + rim notch bridged by a lumped BAR50 sheet.

Geometry choice: `SlottedCircularPatch` has no notch parameter and cannot be given one without changing it (frozen), so `SlotDiscNotch`
below DUCK-TYPES `geometry.Geometry` in the same way (a new class in this script; nothing under src/ is touched).  It reproduces
SlottedCircularPatch exactly when notch_width = 0 (tested in `selftest`).  Not PolyComposite: the circle stays a mask, the same
discretisation as the Phase 26 cells (so the NOCUT state is directly comparable with the Phase 26 table).
Axes: slot long axis along +x (psi = 0); the notch is at the +x end of the slot: |y| <= g_c/2, from the slot end x = L/2 to the disc edge.
The diode sheet (LumpedLoad, direction "y": current across the notch) is g_c across (y) x 0.35 mm along x, centred in the rim; its
region is metal in the pattern (so rooftops cross it), OFF_OPEN removes it (empty notch).
"""
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from src.errors import LwmomError

from src._validate import finite_real
from src.basis_functions import MeshSettings, RooftopBasis
from src.dispersion import tensor_sheet as ts
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal, Substrate
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

MM = 1e-3
F = 5.8e9
K0 = free_space_wavenumber(F)
G_C = 1.00 * MM  # BAR50-02V SC79 footprint: pads 0.35 mm, 1.35 mm c-c (user correction; was 0.65 -> g065)
DIODE_ALONG = 0.35 * MM
THETA_BEAM = math.radians(25.0)


@dataclass(frozen=True)
class SlotDiscNotch:
    radius: float
    slot_length: float
    slot_width: float
    notch_width: float = 0.0  # 0 = no notch (identical to SlottedCircularPatch at angle 0)
    bridge: bool = False  # True: the diode region (g_c x 0.35 mm, centred in the rim) is metal (a lumped load fills it)
    kind = "patch_circle_rect_slot_notch"
    polarity = "metal"
    rotation = 0.0

    def __post_init__(self):
        for n in ("radius", "slot_length", "slot_width", "notch_width"):
            object.__setattr__(self, n, finite_real(n, getattr(self, n)))

    @property
    def diode_x(self):
        return 0.5 * (self.slot_length / 2 + self.radius)

    @property
    def feature_area(self):
        return math.pi * self.radius**2 - self.slot_length * self.slot_width

    def corners(self):
        r = self.radius
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])

    def in_feature(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        inside = np.hypot(x, y) <= self.radius
        if self.slot_length > 0 and self.slot_width > 0:
            inside &= ~((np.abs(x) <= self.slot_length / 2) & (np.abs(y) <= self.slot_width / 2))
        if self.notch_width > 0:
            notch = (x >= self.slot_length / 2) & (np.abs(y) <= self.notch_width / 2)
            if self.bridge:
                notch &= ~((np.abs(x - self.diode_x) <= DIODE_ALONG / 2) & (np.abs(y) <= self.notch_width / 2))
            inside &= ~notch
        return inside


def _z_off_bar50(f):
    w = 2 * np.pi * f
    return 1j * w * 0.6e-9 + 1.0 / (1.0 / 5000.0 + 1j * w * 0.15e-12)


STATES = {
    # name: (notch?, load kwargs or None)
    "NOCUT": (False, None),
    "OFF_OPEN": (True, "open"),
    "OFF_BAR50": (True, dict(z_of_f=_z_off_bar50)),  # (0.15 pF || 5 kohm) in series with 0.6 nH
    "ON_BAR50": (True, dict(R=3.0, L=0.6e-9)),
    "ON_BAR50_MAX": (True, dict(R=4.5, L=0.6e-9)),
}


def make_cell(p_mm, r_mm, L_mm, W_mm, sub, state):
    """UnitCell for one state.  OFF_OPEN: empty notch (bridge False, no load)."""
    notch, ld = STATES[state]
    geo = SlotDiscNotch(r_mm * MM, L_mm * MM, W_mm * MM, G_C if notch else 0.0, bridge=isinstance(ld, dict))
    loads = ()
    if isinstance(ld, dict):
        loads = (LumpedLoad(geo.diode_x, 0.0, length=G_C, width=DIODE_ALONG, direction="y", **ld),)
    return UnitCell(RectangularLattice(p_mm * MM, p_mm * MM), geo, sub, Metal(), backing="pec", loads=loads)


def y_tensor(cell, mesh, kx0, ky0, f=F, basis=None):
    """Bloch sheet tensor with lumped loads (bloch_sheet_tensor + load matrix, as sidequest_cell_search.bloch_loaded)."""
    green = LayeredGreenFunction(cell.substrate, f, backing="pec")
    basis = basis or RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    z = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    if cell.loads:
        z = z + load_matrix(basis, list(cell.loads), green.omega)
    v = np.stack([assemble_excitation_vector(basis, e, kx0, ky0) for e in np.eye(2)], axis=1)
    cur = solve(z, v)
    j = np.stack([basis.current_coefficients(cur[:, k], kx0, ky0) for k in range(2)], axis=1)
    e_tot = np.eye(2) + np.asarray(green.spectral_dyad(kx0, ky0), dtype=complex) @ j
    return j @ np.linalg.inv(e_tot)


# ------------------------------------------------------------------------------------------------ nodes and self-consistent beta
CLASS_STATES = {"NOCUT": ["NOCUT"], "OPEN": ["OFF_OPEN"], "BRIDGE": ["OFF_BAR50", "ON_BAR50", "ON_BAR50_MAX"]}


def make_sub(sub):
    return Substrate(sub[0], sub[1], sub[2], sub[3] * MM)


def node_task(a):
    """(spec, cls, direction 'x'|'y', kt_over_k0_f) -> {state: 2x2 complex Y} for all states of the class (one Z assembly)."""
    spec, cls, direction, kt, n, forder, f = a
    sub, p, r, L, W = spec
    s0 = CLASS_STATES[cls][0]
    cell = make_cell(p, r, L, W, make_sub(sub), s0)
    mesh = MeshSettings(n, n, forder)
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
        c = make_cell(p, r, L, W, make_sub(sub), st)
        z = z0 + load_matrix(basis, list(c.loads), green.omega) if c.loads else z0
        cur = solve(z, v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[st] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(cls=cls, direction=direction, kt=kt, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))


def interp_y(nodes, st, kt):
    """Linear interpolation (extrapolation outside) of Y_par (x-dir node, [0,0]) and Y_perp (y-dir node, [1,1]) plus the off-diagonals, at kt."""
    ks = sorted(nodes)
    def at(direction):
        ys = [nodes[k][direction][st] for k in ks]
        if len(ks) == 1:
            return ys[0]
        i = int(np.clip(np.searchsorted(ks, kt) - 1, 0, len(ks) - 2))
        w = (kt - ks[i]) / (ks[i + 1] - ks[i])
        return (1 - w) * ys[i] + w * ys[i + 1]
    return at("x"), at("y")


def state_summary(nodes, st, sub, f=F, theta=THETA_BEAM, guess=1.05):
    """Self-consistent beta (mean sheet from Y_par(kt), Y_perp(kt) at kt = beta) and derived quantities for one state."""
    s = make_sub(sub)
    k0 = free_space_wavenumber(f)
    kt = guess
    for _ in range(40):
        yx, yy = interp_y(nodes, st, kt)
        ybar = 0.5 * (yx[0, 0].imag + yy[1, 1].imag)
        try:
            b = find_beta_from_sheet_admittance(1j * ybar, s, f).beta_rad_per_m / k0
        except LwmomError as e:  # no bound mode in (k0, 2.5 k0) for this state at this kt: flag it, do not abort the row
            nan = float("nan")
            return dict(beta_over_k0=nan, Y_par_mS=nan, Y_perp_mS=nan, Delta_mS=nan, mean_mS=nan, A_par_mS=nan, A_perp_mS=nan,
                        alpha_np_per_m=nan, extrapolated=False, nearest_node_dist=0.0, status="NO_ROOT", Im_Ybar_mS=ybar * 1e3, kt=kt)
        if abs(b - kt) < 1e-5:
            break
        kt = 0.5 * (kt + b)
    yx, yy = interp_y(nodes, st, b)
    ypar, yperp = yx[0, 0], yy[1, 1]
    delta = ypar - yperp
    mean = 0.5 * (ypar + yperp)
    alpha = leakage_alpha_first_order_tensor(b * k0, 1j * mean.imag, 0.5 * abs(delta.imag), theta, s, f)
    ks = sorted(nodes)
    return dict(beta_over_k0=float(b), Y_par_mS=ypar * 1e3, Y_perp_mS=yperp * 1e3, Delta_mS=delta * 1e3, mean_mS=mean * 1e3,
                A_par_mS=0.5 * (yx[0, 1] - yx[1, 0]) * 1e3, A_perp_mS=0.5 * (yy[0, 1] - yy[1, 0]) * 1e3,
                alpha_np_per_m=float(alpha), extrapolated=bool(b < ks[0] - 1e-9 or b > ks[-1] + 1e-9),
                nearest_node_dist=float(min(abs(b - k) for k in ks)))


def next_node(nodes, states, sub, tol=0.004, f=F):
    """kt of the next node (None if every state's beta is within tol of a node)."""
    bs = []
    for st in states:
        r = state_summary(nodes, st, sub, f)
        if r["nearest_node_dist"] > tol:
            bs.append(r["beta_over_k0"])
    return None if not bs else float(np.median(bs))
