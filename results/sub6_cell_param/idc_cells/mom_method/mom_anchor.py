"""Step 2 gate: repo MoM Bloch extraction of the anchor slot-disc (r 5.45, slot 8.72 x 2.91, P 11.6, RT5880 3.175 mm on PEC) vs the lidded HFSS rows.
HFSS anchor rows (beta/k0 at the eigenfrequency): across 85 deg: 5.553 GHz 1.099; across 95: 5.925 1.151; along 85: 5.743 1.062; along 95: 6.140 (and 6.287).
Repo convention: slot axis along x; phi = 0 along the slot, 90 across it."""
import sys, time
import numpy as np
sys.path.insert(0, r'C:\Users\Bikash\Desktop\lwmom_main')
from src.basis_functions import MeshSettings
from src.dispersion import tensor_sheet as ts
from src.dispersion.slotted_cells import SlottedCircularPatch
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

P = 11.6e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 3.175e-3)
mesh_n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
geo = SlottedCircularPatch(10.9e-3, 8.72e-3, 2.91e-3, 0.0)
cell = UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")
rows = (('across 85', 90.0, 5.553e9, 1.099), ('across 95', 90.0, 5.925e9, 1.151), ('along 85', 0.0, 5.743e9, 1.062))
for name, phi, f, hfss in rows:
    k0 = free_space_wavenumber(f)
    t0 = time.time()
    try:
        b, info = ts.bloch_beta(cell, f, mesh, np.radians(phi), n_grid=7, refine=2)
        print('%s: f %.3f GHz MoM beta/k0 %.4f vs HFSS %.3f  (%.1f %%)  [%.0f s]' % (name, f / 1e9, b / k0, hfss, 100 * (b / k0 / hfss - 1), time.time() - t0), flush=True)
    except Exception as e:
        print(name, 'FAIL', type(e).__name__, str(e)[:150], '[%.0f s]' % (time.time() - t0), flush=True)
