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
SUB = Substrate("rt5880", 2.2, 0.0009, 3.175e-3); K0 = free_space_wavenumber(F)
for (r, L, W) in ((5.45, 8.72, 2.91), (5.45, 8.72, 5.0), (5.65, 8.72, 2.91), (5.65, 9.3, 4.2)):
    for n in (32, 48, 64):
        t0 = time.time(); mesh = MeshSettings(n, n, n)
        cell = UnitCell(RectangularLattice(P, P), SlottedCircularPatch(2*r*1e-3, L*1e-3, W*1e-3, 0.0), SUB, Metal(), backing="pec")
        basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y); kt = 1.10*K0
        yp = ts.bloch_sheet_tensor(cell, F, mesh, kt, 0.0, basis=basis)[0, 0].imag; yq = ts.bloch_sheet_tensor(cell, F, mesh, 0.0, kt, basis=basis)[1, 1].imag
        yb, dl = 0.5*(yp+yq), 0.5*abs(yp-yq); beta = find_beta_from_sheet_admittance(1j*yb, SUB, F).beta_rad_per_m
        a = leakage_alpha_first_order_tensor(beta, 1j*yb, dl, np.radians(30.0), SUB, F)
        print('r %.2f L %.2f W %.2f mesh %d: Ybar %.2f delta %.2f beta/k0 %.3f alpha %.2f [%.0f s]' % (r, L, W, n, yb*1e3, dl*1e3, beta/K0, a, time.time()-t0), flush=True)
