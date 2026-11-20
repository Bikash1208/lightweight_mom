"""Step 3: repo first-order tensor leakage alpha of the anchor slot-disc at 5.68 GHz (r 5.45, slot 8.72 x 2.91, P 11.6, RT5880 3.175 mm on PEC).
Operating wavenumber kt = beta_op * k0 along and across the slot -> Y_par, Y_perp -> Ybar, delta -> alpha (theta 30 deg). MODEL_ONLY."""
import sys, time
import numpy as np
sys.path.insert(0, r'C:\Users\Bikash\Desktop\lwmom_main')
from src.basis_functions import MeshSettings, RooftopBasis
from src.dispersion import tensor_sheet as ts
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.dispersion.slotted_cells import SlottedCircularPatch
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

F, P = 5.68e9, 11.6e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 3.175e-3)
K0 = free_space_wavenumber(F)
mesh_n = int(sys.argv[1]) if len(sys.argv) > 1 else 32
mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
cell = UnitCell(RectangularLattice(P, P), SlottedCircularPatch(10.9e-3, 8.72e-3, 2.91e-3, 0.0), SUB, Metal(), backing="pec")
basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
for kt_over in (1.05, 1.10, 1.16):
    t0 = time.time()
    kt = kt_over * K0
    y_par = ts.bloch_sheet_tensor(cell, F, mesh, kt, 0.0, basis=basis)
    y_perp = ts.bloch_sheet_tensor(cell, F, mesh, 0.0, kt, basis=basis)
    yp, yq = y_par[0, 0].imag, y_perp[1, 1].imag
    ybar, delta = 0.5 * (yp + yq), 0.5 * abs(yp - yq)
    beta = find_beta_from_sheet_admittance(1j * ybar, SUB, F).beta_rad_per_m
    a30 = leakage_alpha_first_order_tensor(beta, 1j * ybar, delta, np.radians(30.0), SUB, F)
    print('kt %.2f k0: Y_par %.2f mS Y_perp %.2f mS Ybar %.2f delta %.2f | beta_mean/k0 %.3f | alpha(30deg) %.2f Np/m -> converted over 150 mm %.2f; two-beam dither alpha ~ alpha/2 = %.2f -> %.2f  [%.0f s]' % (
        kt_over, yp * 1e3, yq * 1e3, ybar * 1e3, delta * 1e3, beta / K0, a30, 1 - np.exp(-2 * a30 * 0.15), a30 / 2, 1 - np.exp(-a30 * 0.15), time.time() - t0), flush=True)
