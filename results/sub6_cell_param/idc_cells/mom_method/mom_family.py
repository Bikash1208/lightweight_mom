import sys, time, math
from dataclasses import dataclass
import numpy as np
import mom_diode as md
from mom_diode import *
@dataclass(frozen=True)
class Mask:
    fn: object
    name: str = 'x'
    kind = "patch_circle_rect_slot"
    @property
    def polarity(self): return "metal"
    @property
    def rotation(self): return 0.0
    @property
    def feature_area(self): return 1e-5
    def corners(self): return np.array([[5.5e-3,5.5e-3],[-5.5e-3,5.5e-3],[-5.5e-3,-5.5e-3],[5.5e-3,-5.5e-3]])
    def in_feature(self, x, y): return self.fn(np.asarray(x, float), np.asarray(y, float))
def disc_slot(r, L, W): return lambda x, y: (np.hypot(x, y) <= r*1e-3) & ~((abs(x) <= L*5e-4) & (abs(y) <= W*5e-4))
def sq_slot(s, L, W): return lambda x, y: ((abs(x) <= s*5e-4) & (abs(y) <= s*5e-4)) & ~((abs(x) <= L*5e-4) & (abs(y) <= W*5e-4))
def dogbone(r, L, W0, W1, Le): return lambda x, y: (np.hypot(x, y) <= r*1e-3) & ~(((abs(x) <= L*5e-4) & (abs(y) <= W0*5e-4)) | ((abs(x) <= L*5e-4) & (abs(x) >= (L/2-Le)*1e-3) & (abs(y) <= W1*5e-4)))
def run(name, fn, sub_t=3.175e-3, n=48):
    md.SUB = Substrate("rt5880", 2.2, 0.0009, sub_t)
    cell = UnitCell(RectangularLattice(P, P), Mask(fn), md.SUB, Metal(), backing="pec")
    t0 = time.time()
    try:
        r = md.evaluate(cell, n)
        print('%-44s Ybar %.2f delta %.2f beta/k0 %.3f alpha %.2f Np/m [%.0f s]' % (name, r['ybar'].imag*1e3, r['delta']*1e3, r['beta_over_k0'], r['alpha'], time.time()-t0), flush=True)
    except Exception as e:
        print(name, 'FAIL', str(e)[:80], flush=True)
