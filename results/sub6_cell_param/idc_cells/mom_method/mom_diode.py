"""MoM (repo src/, read-only worktree) diode-cell tools: slot-disc with a centre strip and a gap carrying a LumpedLoad (true diode Z, not sheet values).
Board frame: slot long axis along x, strip across the slot along y at |x| < strip_w/2, gap |y| < gap/2, load direction 'y'. MODEL_ONLY."""
import math, sys, time
from dataclasses import dataclass
import numpy as np
sys.path.insert(0, r'C:\Users\Bikash\Desktop\lwmom_main')
from src.basis_functions import MeshSettings, RooftopBasis
from src.dispersion import tensor_sheet as ts
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.frequency import free_space_wavenumber, angular_frequency
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.solver import solve

F, P = 5.68e9, 11.6e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 3.175e-3)
K0 = free_space_wavenumber(F)
OMEGA = angular_frequency(F)


@dataclass(frozen=True)
class SlotStripDisc:
    """Disc (diameter) with a slot (length along x, width along y); optional centre strip (width along x, across the slot) with an optional gap (length along y)."""
    diameter: float
    slot_length: float
    slot_width: float
    strip_w: float = 0.0
    gap: float = 0.0
    kind = "patch_circle_rect_slot"

    @property
    def polarity(self): return "metal"
    @property
    def rotation(self): return 0.0
    @property
    def feature_area(self): return math.pi * self.diameter ** 2 / 4.0 - self.slot_length * self.slot_width
    def corners(self):
        r = self.diameter / 2.0
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])
    def in_feature(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        inside = np.hypot(x, y) <= self.diameter / 2.0
        slot = (np.abs(x) <= self.slot_length / 2) & (np.abs(y) <= self.slot_width / 2)
        if self.strip_w > 0:
            strip = (np.abs(x) <= self.strip_w / 2) & (np.abs(y) <= self.slot_width / 2)
            if self.gap > 0:
                strip &= np.abs(y) > self.gap / 2
            slot &= ~strip
        return inside & ~slot


def make_cell(geo, loads=()):
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec", loads=tuple(loads))


def y_tensor(cell, mesh, kx0, ky0, basis=None):
    """Bloch sheet tensor like tensor_sheet.bloch_sheet_tensor, plus the load matrix of cell.loads (the repo's function ignores loads)."""
    from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
    green = LayeredGreenFunction(cell.substrate, F, backing="pec")
    basis = basis or RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    z = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    if cell.loads:
        z = z + load_matrix(basis, cell.loads, OMEGA)
    v = np.stack([assemble_excitation_vector(basis, e, kx0, ky0) for e in np.eye(2)], axis=1)
    cur = solve(z, v)
    j = np.stack([basis.current_coefficients(cur[:, k], kx0, ky0) for k in range(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx0, ky0), dtype=complex)
    return j @ np.linalg.inv(np.eye(2) + g00 @ j)


def evaluate(cell, mesh_n, kt_over=1.10, theta_deg=30.0):
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    kt = kt_over * K0
    yp = y_tensor(cell, mesh, kt, 0.0, basis)[0, 0]       # wave along the slot (x), field along x
    yq = y_tensor(cell, mesh, 0.0, kt, basis)[1, 1]       # wave across the slot (y), field along y
    ybar, delta = 0.5 * (yp + yq), 0.5 * abs(yp - yq)
    beta = find_beta_from_sheet_admittance(ybar, SUB, F).beta_rad_per_m
    a = leakage_alpha_first_order_tensor(beta, ybar, delta, np.radians(theta_deg), SUB, F)
    return dict(y_par=yp, y_perp=yq, ybar=ybar, delta=delta, beta_over_k0=beta / K0, alpha=a)


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    D, L, W = 10.9e-3, 8.72e-3, 2.91e-3
    cases = []
    cases.append(('anchor (no strip)', make_cell(SlotStripDisc(D, L, W))))
    cases.append(('strip, continuous', make_cell(SlotStripDisc(D, L, W, 0.5e-3, 0.0))))
    g = SlotStripDisc(D, L, W, 0.5e-3, 0.3e-3)
    cases.append(('strip + gap, load R=1e-3 ohm (expect = continuous)', make_cell(g, (LumpedLoad(0, 0, 0.3e-3, 0.5e-3, 'y', R=1e-3),))))
    cases.append(('strip + gap, open_circuit load (expect = open gap)', make_cell(g, (LumpedLoad(0, 0, 0.3e-3, 0.5e-3, 'y', open_circuit=True),))))
    cases.append(('strip + gap, no load (open gap, reference)', make_cell(g)))
    for name, cell in cases:
        t0 = time.time()
        r = evaluate(cell, n)
        print('%-52s Ypar %.2f%+.2fj Yperp %.2f%+.2fj mS | Ybar %.2f delta %.2f beta/k0 %.3f alpha %.2f [%.0f s]' % (
            name, r['y_par'].real * 1e3, r['y_par'].imag * 1e3, r['y_perp'].real * 1e3, r['y_perp'].imag * 1e3, r['ybar'].imag * 1e3, r['delta'] * 1e3, r['beta_over_k0'], r['alpha'], time.time() - t0), flush=True)
