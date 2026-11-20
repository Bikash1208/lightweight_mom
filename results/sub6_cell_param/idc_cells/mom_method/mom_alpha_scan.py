"""alpha map of the static slot-disc family at 5.68 GHz with the repo MoM + first-order tensor leakage (kt = 1.10 k0). MODEL_ONLY."""
import sys, time, json
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
out = []
for r in (5.45, 5.65):
    for L in (6.0, 7.0, 8.0, 8.72, 9.3):
        if L > 2 * (r - 1.0) + 1e-9:
            continue
        for W in (1.5, 2.0, 2.91, 3.5, 4.2, 5.0):
            chord = 2 * np.sqrt(r * r - (W / 2) ** 2)
            if L > chord - 0.3:
                continue
            try:
                cell = UnitCell(RectangularLattice(P, P), SlottedCircularPatch(2 * r * 1e-3, L * 1e-3, W * 1e-3, 0.0), SUB, Metal(), backing="pec")
                basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
                kt = 1.10 * K0
                yp = ts.bloch_sheet_tensor(cell, F, mesh, kt, 0.0, basis=basis)[0, 0].imag
                yq = ts.bloch_sheet_tensor(cell, F, mesh, 0.0, kt, basis=basis)[1, 1].imag
                ybar, delta = 0.5 * (yp + yq), 0.5 * abs(yp - yq)
                beta = find_beta_from_sheet_admittance(1j * ybar, SUB, F).beta_rad_per_m
                a = leakage_alpha_first_order_tensor(beta, 1j * ybar, delta, np.radians(30.0), SUB, F)
                out.append(dict(r=r, L=L, W=W, ybar_mS=ybar * 1e3, delta_mS=delta * 1e3, beta_over_k0=beta / K0, alpha=a))
                print('r %.2f L %.2f W %.2f: Ybar %.2f delta %.2f mS beta/k0 %.3f alpha %.2f Np/m (conv %.2f)' % (r, L, W, ybar * 1e3, delta * 1e3, beta / K0, a, 1 - np.exp(-2 * a * 0.15)), flush=True)
            except Exception as e:
                print('r %.2f L %.2f W %.2f FAIL %s' % (r, L, W, str(e)[:80]), flush=True)
json.dump(out, open(r'C:\Users\Bikash\Desktop\lwmom_labB\bar90_tune_work\mom_alpha_scan.json', 'w'), indent=1)
