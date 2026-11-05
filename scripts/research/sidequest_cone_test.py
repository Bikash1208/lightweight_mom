"""Cone / ring test: azimuthal median level of the MODEL far field at theta 12 / 24 deg relative to the peak, for the baseline board (RT5880, lam 0.28, mu 30): full, static-only, switched-only, and with the statics scaled (decoupled static depth).
HFSS (ffd): median at theta 12 = -4.6 / -4.2 dB, at theta 24 = -5.7 / -6.3 dB below the peak (state 1 / 2).   MODEL_ONLY.   -> results/sidequest_hybrid/k2_cone_test.json"""
import json
import math
import os
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
os.environ.setdefault("K4_A0", "0.844")
os.environ.setdefault("K4_BETA_K0", "1.0731")
import numpy as np
import sidequest_coded_k4 as K
import sidequest_k2_default_design as D

st = D.st


def pat(bd, amp, m, k, mask=None):
    j, _ = bd.unit(amp)
    mm = m[k] if mask is None else m[k] * mask[None, :]
    cur = (mm[0] * j["x"][0] + mm[1] * j["y"][0], mm[0] * j["x"][1] + mm[1] * j["y"][1])
    et, ep = st.far_field(bd.x, bd.y, *cur)
    return np.abs(et) ** 2 + np.abs(ep) ** 2


def ring(u, ref, th):
    i = int(np.argmin(np.abs(np.degrees(st.TG[:, 0]) - th)))
    row = u[i]
    return float(10 * np.log10(np.median(row) / ref + 1e-30)), float(10 * np.log10(row.max() / ref + 1e-30))


if __name__ == "__main__":
    stat, _p, _n = D.library()
    levels = sorted({round(v[0], 4) for v in stat.values()})
    D.setup(10.5e9, 0.15, D.MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    r = D.run_design(bd, D.options(levels, [(1.0, -1.0, 1.0)], [(0, 1)]), [0.0, np.radians(180)], 0.28, 2.0, mu=30.0)
    m0, amp, nd = r["m"], r["amp"], r["nd"]
    sw = nd > 0
    out = []
    for s in (1.0, 1.33, 2.0):
        m = m0.copy()
        m[:, :, ~sw] *= s
        for k in range(2):
            full = pat(bd, amp, m, k)
            ref = full.max()
            rec = {"static_scale": s, "state": k + 1}
            for nm, mask in (("full", None), ("static_only", (~sw).astype(float)), ("switched_only", sw.astype(float))):
                u = full if mask is None else pat(bd, amp, m, k, mask)
                rec[nm] = {f"theta{t}": ring(u, ref, t) for t in (12, 24)}
            out.append(rec)
            print("s=%.2f state %d  full median/max dB rel peak: th12 %s th24 %s | static-only th12 %s th24 %s | switched-only th12 %s th24 %s" % (s, k + 1, *[tuple(round(x, 1) for x in rec[n][t]) for n in ("full", "static_only", "switched_only") for t in ("theta12", "theta24")]), flush=True)
    Path("results/sidequest_hybrid/k2_cone_test.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline=chr(10))
