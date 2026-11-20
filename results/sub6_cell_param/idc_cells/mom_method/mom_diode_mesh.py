import sys, time
import numpy as np
from mom_diode import *
D, L, W = 10.9e-3, 8.72e-3, 2.91e-3
strip = SlotStripDisc(D, L, W, 0.5e-3, 0.0)
zoff = lambda f: 2j*np.pi*f*0.2e-9 + 1.0/(2j*np.pi*f*0.19e-12 + 1.0/(1.5e3*(2.5e9/f)))
for n in (48, 64, 80):
    for name, loads in (('open gap', (LumpedLoad(0,0,0.3e-3,0.5e-3,'y',open_circuit=True),)), ('BAR90 reverse (0.19 pF)', (LumpedLoad(0,0,0.3e-3,0.5e-3,'y',z_of_f=zoff),)), ('pure -j140', (LumpedLoad(0,0,0.3e-3,0.5e-3,'y',C=1/(2*np.pi*F*140)),))):
        t0 = time.time(); r = evaluate(make_cell(strip, loads), n)
        print('mesh %d %-26s Yperp %.2f%+.2fj Ypar %.2fj | Ybar %.2f delta %.2f beta/k0 %.3f alpha %.2f [%.0f s]' % (n, name, r['y_perp'].real*1e3, r['y_perp'].imag*1e3, r['y_par'].imag*1e3, r['ybar'].imag*1e3, r['delta']*1e3, r['beta_over_k0'], r['alpha'], time.time()-t0), flush=True)
