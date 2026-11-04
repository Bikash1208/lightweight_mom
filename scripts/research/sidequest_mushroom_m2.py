"""Calibrated wire-medium model 'M2' (L_eff scaled by kappa so that the ON lowest band at 180 deg = 4.92 GHz, HFSS) and the alternative 'M1' (plain slab + series L_br = 6.0 nH).
M2 keeps the local-resonance character (flat band at f0, beta-independent via branch) that the HFSS ON bands show (flat 4.9-5.1 GHz for 90-180 deg); M1 reproduces the HFSS ON second band
(7.66 / 11.2 vs 7.74 / 11.14 GHz) but disperses (3.3 -> 4.9 GHz from 90 to 180 deg). The two are reported as the model-form spread."""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
from scipy.optimize import brentq

import sidequest_mushroom_calibrate as cal
import sidequest_mushroom_cell as mc

FREQS = [4, 5, 6, 8, 9.5, 10.54, 11.5, 12.5, 14, 16, 18]
LID = 30e-3


class ModelM2(mc.Model):
    def __init__(self, table, kappa):
        super().__init__(table)
        self.kappa = kappa

    def z_slab_tem(self, f):
        return self.kappa * super().z_slab_tem(f)


def lowest(m, state, ph):
    fs = cal.surface_band(m, state, ph)
    return min(fs) / 1e9 if fs else None


def sheet_job(item):
    name, (w, g) = item
    return name, mc.sheet_table(w * 1e-3, FREQS, 32, g * 1e-3, 0.8e-3)


if __name__ == "__main__":
    base = json.loads(Path("results/sidequest_hybrid/mushroom_theory_raw.json").read_text(encoding="utf-8"))["5.6"]["table"]
    base = mc.sheet_table(5.6e-3, FREQS, 32, 0.4e-3, 0.8e-3)  # same grid as the variants
    # calibrate kappa on f_ON(180) = 4.92 GHz
    def r(k):
        v = lowest(ModelM2(base, k), "ON", 180)
        return (v if v else 30.0) - 4.92
    ks = np.geomspace(1.0, 12.0, 30)
    vs = [r(k) for k in ks]
    idx = [i for i in range(len(ks) - 1) if np.sign(vs[i]) != np.sign(vs[i + 1])]
    kappa = brentq(r, ks[idx[0]], ks[idx[0] + 1])
    m = ModelM2(base, kappa)
    print("kappa", round(kappa, 3), "L_eff", round(kappa * 2.18, 2), "nH")
    print("M2 baseline ON lowest band:", {ph: round(lowest(m, "ON", ph), 2) for ph in (60, 90, 120, 150, 180)}, " HFSS: 90-180 flat 4.9-5.1")
    print("M2 baseline OFF:", {ph: [round(x / 1e9, 2) for x in cal.surface_band(m, "OFF", ph)][-1:] for ph in (30, 60, 90, 120, 150, 180)}, "beta", [round(b / mc.free_space_wavenumber(10.54e9), 3) for b in m.beta_at("OFF", 10.54e9, None) if b / mc.free_space_wavenumber(10.54e9) < 1.4][:1], " HFSS 3.86/7.66/10.98/12.54/13.17/13.36, 1.06")
    geoms = {"w3.2_g0.4": (3.2, 0.4), "w3.6_g0.4": (3.6, 0.4), "w4.0_g0.4": (4.0, 0.4), "w4.4_g0.4": (4.4, 0.4), "w3.6_g0.8": (3.6, 0.8), "w4.4_g0.8": (4.4, 0.8)}
    with ProcessPoolExecutor(4) as ex:
        tabs = dict(ex.map(sheet_job, geoms.items()))
    out = {"kappa": kappa, "L_eff_nH": kappa * 2.18, "geoms": {}}
    for name, tab in tabs.items():
        m2 = ModelM2(tab, kappa)
        m1 = cal.ModelPB(tab, 6.0e-9)
        row = {"b_10p54_mS": m2.b_sheet(10.54e9) * 1e3}
        for tag, mm in (("M2", m2), ("M1", m1)):
            row[tag] = {"ON_90": lowest(mm, "ON", 90), "ON_180": lowest(mm, "ON", 180),
                        "OFF_beta": [b / mc.free_space_wavenumber(10.54e9) for b in mm.beta_at("OFF", 10.54e9, None) if b / mc.free_space_wavenumber(10.54e9) < 1.4][:1]}
        out["geoms"][name] = row
        print(name, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items() if k == "b_10p54_mS"}, {t: {k: (round(x, 2) if isinstance(x, float) else x) for k, x in row[t].items()} for t in ("M2", "M1")}, flush=True)
    Path("results/sidequest_hybrid/mushroom_m2_scan.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
