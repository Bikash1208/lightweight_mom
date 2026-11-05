"""Does 'static cells at full depth, switched cells at reduced contrast' reproduce the HFSS far field of board_K2_default?  Fixed baseline design (RT5880, lam 0.28, mu 30); the physics is perturbed: statics m x s, switched cells m x c (both states).
HFSS targets (ffd analysis): state 1: peak 9.47 dBi at (24,2), strongest sidelobe -0.7 dB (20,24), D(25,180) -10.2 dB below the peak; state 2: 8.27 dBi at (24,176), sidelobe -0.1 dB (20,154), D(25,0) -14.9 dB below the peak; common-mode power 97 %.  MODEL_ONLY.
    -> results/sidequest_hybrid/k2_decouple_test.json"""
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


def pattern(bd, amp, m, state):
    j, _ = bd.unit(amp)
    mm = m[state]
    cur = (mm[0] * j["x"][0] + mm[1] * j["y"][0], mm[0] * j["x"][1] + mm[1] * j["y"][1])
    et, ep = st.far_field(bd.x, bd.y, *cur)
    return np.abs(et) ** 2 + np.abs(ep) ** 2


def metrics(u, k):
    ref = u.max()
    ipk = np.unravel_index(np.argmax(u), u.shape)
    cosd = np.cos(st.TG) * np.cos(st.TG[ipk]) + np.sin(st.TG) * np.sin(st.TG[ipk]) * np.cos(st.PG - st.PG[ipk])
    dist = np.degrees(np.arccos(np.clip(cosd, -1, 1)))
    prad = float(np.sum(u * np.sin(st.TG) * math.radians(1) * math.radians(2)))
    hp = math.degrees(math.sqrt(16 / (4 * math.pi * ref / prad)))
    sl = np.where(dist > 1.2 * hp, u, 0)
    isl = np.unravel_index(np.argmax(sl), sl.shape)
    oth = 180 if k == 0 else 0
    ii = np.argmin(np.abs(np.degrees(st.TG[:, 0]) - 25))
    jj = np.argmin(np.abs(np.degrees(st.PG[0]) - oth))
    return {"peak_dir": [float(np.degrees(st.TG[ipk])), float(np.degrees(st.PG[ipk]))], "SLL_db": float(10 * np.log10(sl[isl] / ref)), "SLL_dir": [float(np.degrees(st.TG[isl])), float(np.degrees(st.PG[isl]))], "other_beam_rel_db": float(10 * np.log10(u[ii, jj] / ref)), "D_dbi": float(10 * np.log10(4 * math.pi * ref / prad))}


if __name__ == "__main__":
    stat, _p, _n = D.library()
    levels = sorted({round(v[0], 4) for v in stat.values()})
    D.setup(10.5e9, 0.15, D.MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    r = D.run_design(bd, D.options(levels, [(1.0, -1.0, 1.0)], [(0, 1)]), [0.0, np.radians(180)], 0.28, 2.0, mu=30.0)
    m0, amp, nd = r["m"], r["amp"], r["nd"]
    sw = nd > 0
    rows = []
    for s in (1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0):
        m = m0.copy()
        m[:, :, ~sw] *= s
        out = {"static_scale": s, "states": []}
        for k in range(2):
            u = pattern(bd, amp, m, k)
            out["states"].append(metrics(u, k))
        rows.append(out)
        print("static x%.2f:" % s, ["SLL %.1f other %.1f" % (o["SLL_db"], o["other_beam_rel_db"]) for o in out["states"]], ["pk (%.0f,%.0f) sl (%.0f,%.0f)" % (*o["peak_dir"], *o["SLL_dir"]) for o in out["states"]], flush=True)
    # static-only pattern lobe structure at s = 1: ring-like?
    u = pattern(bd, amp, m0 * (~sw)[None, None, :], 0)
    ref = pattern(bd, amp, m0, 0).max()
    ipk = np.unravel_index(np.argmax(u), u.shape)
    ring = {}
    for th in (12, 20, 26, 32):
        ii = np.argmin(np.abs(np.degrees(st.TG[:, 0]) - th))
        ring[th] = [float(10 * np.log10(u[ii].max() / ref)), float(10 * np.log10(u[ii].min() / ref + 1e-30))]
    print("static-only (state 1) rel. to full peak: max / min over phi at theta 12/20/26/32:", {k: [round(x, 1) for x in v] for k, v in ring.items()})
    Path("results/sidequest_hybrid/k2_decouple_test.json").write_text(json.dumps({"rows": rows, "static_only_ring": ring}, indent=1), encoding="utf-8", newline=chr(10))
