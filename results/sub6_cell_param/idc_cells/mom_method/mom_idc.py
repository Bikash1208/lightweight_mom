import sys, time, math
import numpy as np
import mom_diode as md
from mom_diode import *
from mom_family import Mask
def mn_mask(r, w, n, t, d, bridge, edge=True, ang=0.0):
    c, s = math.cos(ang), math.sin(ang)
    x0 = -(2*n*t)/2.0
    boxes = []
    for i in range(2*n+1):
        xi = x0 + i*t; boxes.append((xi - w/2, xi + w/2, -d/2 - w/2, d/2 + w/2))
    for i in range(2*n):
        yc = d/2 if i % 2 == 0 else -d/2
        boxes.append((x0 + i*t - w/2, x0 + i*t + t + w/2, yc - w/2, yc + w/2))
    if edge:
        xe = x0 + 2*n*t; rimx = math.sqrt(r*r - (d/2 + w/2)**2)
        xs = -(rimx - bridge) if bridge > 0 else -(r + 1); xr = (rimx - bridge) if bridge > 0 else (r + 1)
        boxes.append((xs, x0 + w/2, -d/2 - w/2, -d/2 + w/2))
        ylast = -(d/2 if (2*n - 1) % 2 == 0 else -d/2)
        boxes.append((xe - w/2, xr, ylast - w/2, ylast + w/2))
    def fn(x, y):
        u = (x*c + y*s)*1e3; v = (-x*s + y*c)*1e3
        inside = np.hypot(x, y) <= r*1e-3
        cut = np.zeros_like(inside)
        for (a, b, c0, d0) in boxes: cut |= (u >= a) & (u <= b) & (v >= c0) & (v <= d0)
        return inside & ~cut
    return fn
def run(name, fn, n):
    cell = UnitCell(RectangularLattice(P, P), Mask(fn), md.SUB, Metal(), backing="pec")
    t0 = time.time()
    try:
        r = md.evaluate(cell, n)
        print('%-40s mesh %d Ybar %.2f delta %.2f beta/k0 %.3f alpha %.3f [%.0f s]' % (name, n, r['ybar'].imag*1e3, r['delta']*1e3, r['beta_over_k0'], r['alpha'], time.time()-t0), flush=True)
    except Exception as e: print(name, n, 'FAIL', str(e)[:90], flush=True)
if __name__ == '__main__':
    for n in (48, 64, 80):
        run('D3p6 bridged B2.2 (r5.6)', mn_mask(5.6, 0.4, 2, 1.2, 3.6, 2.2), n)
        run('D3p6 OPEN twin (no bridge)', mn_mask(5.6, 0.4, 2, 1.2, 3.6, 0.0), n)
