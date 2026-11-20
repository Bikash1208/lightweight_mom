import sys, time
import numpy as np
import mom_diode as md
from mom_diode import *
D = 10.9e-3
open_gap = lambda: (LumpedLoad(0,0,0.3e-3,0.5e-3,'y',open_circuit=True),)
for n in (48, 64):
    for W in (2.91, 3.5, 4.2, 5.0):
        geo = SlotStripDisc(D, 8.72e-3, W*1e-3, 0.5e-3, 0.0)
        t0 = time.time()
        a = evaluate(make_cell(SlotStripDisc(D, 8.72e-3, W*1e-3)), n)
        o = evaluate(make_cell(geo, open_gap()), n)
        s = evaluate(make_cell(geo), n)
        print('mesh %d W %.2f: anchor alpha %.2f beta %.3f | open-gap (on proxy) alpha %.2f beta %.3f | continuous strip (off proxy) alpha %.3f -> leak ratio sqrt(a_s/a_anchor) %.2f, sqrt(a_s/a_open) %.2f [%.0f s]' % (
            n, W, a['alpha'], a['beta_over_k0'], o['alpha'], o['beta_over_k0'], s['alpha'], np.sqrt(s['alpha']/a['alpha']), np.sqrt(s['alpha']/o['alpha']), time.time()-t0), flush=True)
