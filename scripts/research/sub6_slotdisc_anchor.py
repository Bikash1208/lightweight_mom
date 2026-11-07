"""Step 0 anchor: the Phase 26 10.54 GHz cell that went to HFSS (RT5880 1.575, p 6.4, r 3.0, slot 4.8 x 1.6, the cell used 869 times in
hex150_tensor_cp_single) and the same cell with all lengths scaled by 10.54/5.8 at 5.8 GHz on RT5880 3.175 mm.  NOCUT, no diode.
MODEL_ONLY.  mesh 64, F = 2n = 128.   python sub6_slotdisc_anchor.py"""
import json, os, sys, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
from multiprocessing import Pool
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import sub6_slotdisc_lib as S

N, FO = 64, 128
s = 10.54 / 5.8
CASES = {
    "orig_10p54GHz": dict(f=10.54e9, sub=("rt5880", 2.2, 0.0009, 1.575), p=6.4, r=3.0, L=4.8, W=1.6),
    "scaled_5p8GHz": dict(f=5.8e9, sub=("rt5880", 2.2, 0.0009, 3.175), p=6.4 * s, r=3.0 * s, L=4.8 * s, W=1.6 * s),
}

def jsonable(o):
    if isinstance(o, (complex, np.complexfloating)): return [float(o.real), float(o.imag)]
    if isinstance(o, dict): return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [jsonable(v) for v in o]
    if isinstance(o, np.generic): return o.item()
    return o

if __name__ == "__main__":
    t0 = time.time()
    res = {}
    nodes = {k: {} for k in CASES}
    with Pool(4) as pool:
        kts = {k: 1.05 for k in CASES}
        for rnd in range(4):
            tasks, meta = [], []
            for k, c in CASES.items():
                if kts[k] is None: continue
                for d in "xy":
                    tasks.append(((c["sub"], c["p"], c["r"], c["L"], c["W"]), "NOCUT", d, kts[k], N, FO, c["f"])); meta.append(k)
            if not tasks: break
            outs = pool.map(S.node_task, tasks)
            for k, o in zip(meta, outs):
                nodes[k].setdefault(o["kt"], {})[o["direction"]] = o["Y"]
                print(rnd, k, o["direction"], "kt %.4f" % o["kt"], "%.0f s" % o["seconds"], "unknowns", o["unknowns"], flush=True)
            for k, c in CASES.items():
                if kts[k] is None: continue
                kts[k] = S.next_node(nodes[k], ["NOCUT"], c["sub"], f=c["f"])
                print("next node", k, kts[k], flush=True)
    for k, c in CASES.items():
        r = S.state_summary(nodes[k], "NOCUT", c["sub"], f=c["f"])
        r30 = S.state_summary(nodes[k], "NOCUT", c["sub"], f=c["f"], theta=np.radians(30.0))
        r["alpha_np_per_m_theta30"] = r30["alpha_np_per_m"]
        res[k] = dict(case=c, result=r, nodes=sorted(nodes[k]))
        print(k, json.dumps(jsonable(r)), flush=True)
    res["seconds"] = time.time() - t0
    Path(__file__).resolve().parents[2].joinpath("results/sub6_slotdisc/anchor.json").write_text(json.dumps(jsonable(res), indent=1))
    a0, a1 = res["orig_10p54GHz"]["result"]["alpha_np_per_m"], res["scaled_5p8GHz"]["result"]["alpha_np_per_m"]
    print("alpha orig (25 deg) %.3f ; x5.8/10.54 = %.3f ; scaled-cell alpha %.3f ; ratio %.3f -> %s" % (a0, a0 * 5.8 / 10.54, a1, a1 / (a0 * 5.8 / 10.54),
          "ANCHOR_REPRODUCED" if 0.5 <= a1 / (a0 * 5.8 / 10.54) <= 2.0 else "ANCHOR_NOT_REPRODUCED"))
