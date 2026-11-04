"""Recalibrate the mushroom-cell model on the HFSS eigenmode result (user run, PEC lid 30 mm, lid modes removed). MODEL_ONLY.
HFSS OFF (30..180 deg step 30): 3.86 / 7.66 / 10.98 / 12.54 / 13.17 / 13.36 GHz, beta/k0 ~ 1.06 at 10.54 GHz.
HFSS ON : lowest band flat at 4.9-5.1 GHz for 90-180 deg (f(180) = 4.92); second band 7.74 (60 deg) -> 11.14 (90 deg).
The first model (wire-medium slab, L_eff = 2.18 nH) put the ON resonance at 10.98 GHz: wrong by 2.2x.  New model 'PB' (plain slab + series branch):
    Y_down = Y_plain(beta) + 1 / (j w L_br + Z_diode(f))      (no wire-medium short; L_br = effective series inductance of via + pad/bridge path, CALIBRATED).
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
from scipy.optimize import brentq

import sidequest_mushroom_cell as mc

HFSS_OFF = {30: 3.86, 60: 7.66, 90: 10.98, 120: 12.54, 150: 13.17, 180: 13.36}
LID = 30e-3


class ModelPB(mc.Model):
    def __init__(self, table, l_br):
        super().__init__(table)
        self.l_br = l_br  # H

    def y_down(self, state, f, beta):
        w = 2 * np.pi * f
        zd = mc.dm.z_of_f("ON_10mA" if state == "ON" else "OFF_0V")(f)
        return self.y_layer(mc.SUB.eps_r, mc.H, f, beta) + 1.0 / (1j * w * self.l_br + zd)


def surface_band(m, state, ph, lid=LID, below_light=True):
    ll = mc.C0 * np.radians(ph) / mc.P / (2 * np.pi)
    fs = m.f_of_phase(state, ph, lid, fmin=0.3e9, fmax=20e9, n=2500)
    fs = [x for x in fs if x < ll - 1e6] if below_light else fs
    return fs


if __name__ == "__main__":
    import json

    tab = json.loads(Path("results/sidequest_hybrid/mushroom_theory_raw.json").read_text(encoding="utf-8"))["5.6"]["table"]
    # calibration: f_ON(180 deg) = 4.92 GHz  -> L_br
    def resid(lb):
        m = ModelPB(tab, lb)
        fs = surface_band(m, "ON", 180)
        return (min(fs) / 1e9 if fs else 30.0) - 4.92
    grid = np.geomspace(0.5e-9, 40e-9, 40)
    vals = [resid(x) for x in grid]
    print([(round(g * 1e9, 2), round(v, 2)) for g, v in zip(grid, vals)][::3])
    idx = [i for i in range(len(grid) - 1) if np.sign(vals[i]) != np.sign(vals[i + 1])]
    for i in idx:
        lb = brentq(resid, grid[i], grid[i + 1])
        m = ModelPB(tab, lb)
        print("L_br = %.2f nH" % (lb * 1e9))
        print(" ON  :", {ph: [round(x / 1e9, 2) for x in surface_band(m, "ON", ph)][-2:] for ph in (60, 90, 120, 180)})
        print(" OFF :", {ph: [round(x / 1e9, 2) for x in surface_band(m, "OFF", ph)][-1:] for ph in (30, 60, 90, 120, 150, 180)}, "HFSS", HFSS_OFF)
        print(" OFF beta/k0 @10.54:", [round(b / mc.free_space_wavenumber(10.54e9), 3) for b in m.beta_at("OFF", 10.54e9, None)][:3])
