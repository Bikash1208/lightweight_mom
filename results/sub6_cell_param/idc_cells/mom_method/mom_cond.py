import sys, time
import numpy as np
from mom_diode import *
D, L, W = 10.9e-3, 8.72e-3, 2.91e-3
strip = SlotStripDisc(D, L, W, 0.5e-3, 0.0)
cases = (('continuous strip', ()), ('BAR90 on 1.5 ohm + 0.2 nH', (LumpedLoad(0,0,0.3e-3,0.5e-3,'y',R=1.5,L=0.2e-9),)), ('0.8 nH sheet-misuse 1.6 + j29', (LumpedLoad(0,0,0.3e-3,0.5e-3,'y',R=1.6,L=29.2/OMEGA),)))
for n in (48, 64, 80):
    for name, loads in cases:
        t0 = time.time(); r = evaluate(make_cell(strip, loads), n)
        print('mesh %d %-30s Yperp %.2f%+.2fj Ypar %.2fj | Ybar %.2f delta %.2f beta/k0 %.3f alpha %.3f [%.0f s]' % (n, name, r['y_perp'].real*1e3, r['y_perp'].imag*1e3, r['y_par'].imag*1e3, r['ybar'].imag*1e3, r['delta']*1e3, r['beta_over_k0'], r['alpha'], time.time()-t0), flush=True)
