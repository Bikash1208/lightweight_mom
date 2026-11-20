import sys, time
import numpy as np
sys.argv = [sys.argv[0]] + sys.argv[1:]
from mom_diode import *
n = int(sys.argv[1]) if len(sys.argv) > 1 else 48
D, L, W = 10.9e-3, 8.72e-3, 2.91e-3
strip = SlotStripDisc(D, L, W, 0.5e-3, 0.0)   # metal strip continuous, the load region sits ON the metal
def ld(**kw): return (LumpedLoad(0, 0, 0.3e-3, 0.5e-3, 'y', **kw),)
cases = [('continuous strip, no load', make_cell(strip)),
         ('strip + load R=1e-3 (expect continuous)', make_cell(strip, ld(R=1e-3))),
         ('strip + open_circuit load (expect open gap: Yperp 10.47j)', make_cell(strip, ld(open_circuit=True))),
         ('strip + load 1e6 ohm (expect ~open gap)', make_cell(strip, ld(R=1e6))),
         ('strip + load 14.4 - j140 (BAR90 off, 0.19 pF||1.5k + 0.2nH)', make_cell(strip, (LumpedLoad(0, 0, 0.3e-3, 0.5e-3, 'y', z_of_f=lambda f: 2j*np.pi*f*0.2e-9 + 1.0/(2j*np.pi*f*0.19e-12 + 1.0/(1.5e3*(2.5e9/f))),),))),
         ('strip + load 1.5 + 0.2 nH (BAR90 on)', make_cell(strip, ld(R=1.5, L=0.2e-9)))]
for name, cell in cases:
    t0 = time.time(); r = evaluate(cell, n)
    print('%-62s Ypar %.2f%+.2fj Yperp %.2f%+.2fj mS | Ybar %.2f delta %.2f beta/k0 %.3f alpha %.2f [%.0f s]' % (name, r['y_par'].real*1e3, r['y_par'].imag*1e3, r['y_perp'].real*1e3, r['y_perp'].imag*1e3, r['ybar'].imag*1e3, r['delta']*1e3, r['beta_over_k0'], r['alpha'], time.time()-t0), flush=True)
