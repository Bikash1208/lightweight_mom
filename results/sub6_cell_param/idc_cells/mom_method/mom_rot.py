"""Rotation covariance of the cell tensor (repo phase26_slotted_family.rotation_check pattern): Y(cell rotated g, k along x) vs R(g) Y(cell 0, k at -g) R(g)^T. MODEL_ONLY."""
import sys, time, math
import numpy as np
import mom_diode as md
from mom_diode import *
from mom_family import Mask
from src.dispersion import tensor_sheet as ts
def cell_mask(r, L, W, strip, gap, ang):
    c, s = math.cos(ang), math.sin(ang)
    def fn(x, y):
        u = (x*c + y*s)*1e3; v = (-x*s + y*c)*1e3
        inside = np.hypot(x, y) <= r*1e-3
        slot = (abs(u) <= L/2) & (abs(v) <= W/2)
        if strip > 0:
            st = (abs(u) <= strip/2) & (abs(v) <= W/2)
            if gap > 0: st &= abs(v) > gap/2
            slot &= ~st
        return inside & ~slot
    return fn
n = int(sys.argv[1]) if len(sys.argv) > 1 else 64
mesh = MeshSettings(n, n, n); KT = 1.10 * K0
def cell(fn): return UnitCell(RectangularLattice(P, P), Mask(fn), md.SUB, Metal(), backing="pec")
for name, (L, W, strip, gap) in (('anchor W2.91', (8.72, 2.91, 0, 0)), ('open-gap strip W2.91', (8.72, 2.91, 0.5, 0.3)), ('open-gap strip W5.0', (8.72, 5.0, 0.5, 0.3))):
    c0 = cell(cell_mask(5.45, L, W, strip, gap, 0.0)); b0 = RooftopBasis(c0, mesh.cells_x, mesh.cells_y)
    for gdeg in (45.0,):
        t0 = time.time(); g = math.radians(gdeg)
        ca = cell(cell_mask(5.45, L, W, strip, gap, g)); ya = ts.bloch_sheet_tensor(ca, F, mesh, KT, 0.0)
        y0 = ts.bloch_sheet_tensor(c0, F, mesh, KT*math.cos(-g), KT*math.sin(-g), basis=b0)
        pred = ts.rotate_tensor(y0, g)
        err = np.linalg.norm(ya - pred) / np.linalg.norm(ya)
        print('%-22s rot %2.0f deg: |Y_rot - R Y0(k at -g) R^T| / |Y_rot| = %.3f   [Yxx rot %.2f vs pred %.2f mS; Yyy %.2f vs %.2f] [%.0f s]' % (name, gdeg, err, ya[0,0].imag*1e3, pred[0,0].imag*1e3, ya[1,1].imag*1e3, pred[1,1].imag*1e3, time.time()-t0), flush=True)
