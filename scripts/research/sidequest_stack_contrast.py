"""Board cell (K2 ls 2.8 rim disc, OFF = open slit, ON = 2.44 ohm R) on other stacks at 10.5 GHz, 6.4 lattice, mesh 64: Y OFF/ON -> X_s (lid) and contrast; leakage scaling A0 ~ (dX)^2.  MODEL_ONLY.  -> results/sidequest_hybrid/k2_stack_contrast.json"""
import json, math, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); warnings.simplefilter("ignore")
import numpy as np
import sidequest_gate_cell as gc
import sidequest_k2_default_lib as L
import sidequest_leaf_cell as L1
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Substrate
F = 10.5e9
STK = {"RT5880_1.575": (2.2, 0.0009, 1.575), "RO3006_1.27": (6.15, 0.002, 1.27), "RO3010_0.635": (11.2, 0.0022, 0.635), "RO3010_1.27": (11.2, 0.0022, 1.27)}
def job(a):
    warnings.simplefilter("ignore")
    n, mode = a; e, t, h = STK[n]
    gc.SUB = Substrate(n, e, t, h * 1e-3)
    y = sheet_tensor(L.cell(2.8, 0.8, mode), F, 64)
    return {"stack": n, "mode": mode, "Y_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3]}
if __name__ == "__main__":
    jobs = [(n, m) for n in STK for m in ("open", "ON")]
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(job, jobs))
    out = {}
    for n, (e, t, h) in STK.items():
        sub = Substrate(n, e, t, h * 1e-3); H = 3e8 / F / 4
        d = {}
        for r in res:
            if r["stack"] == n:
                try: d[r["mode"]] = {**r, **L1.xs_from_y(complex(r["Y_mS"][0] * 1e-3, r["Y_mS"][1] * 1e-3), sub, F, H)}
                except Exception as ex_: d[r["mode"]] = {**r, "err": str(ex_)}
        out[n] = d
        o, on = d["open"], d["ON"]
        print(n, "Y OFF %.2f%+.2fj ON %.2f%+.2fj mS" % (*o["Y_mS"], *on["Y_mS"]), "X_lid OFF %s ON %s" % (round(o.get("X_lid_ohm", float("nan")), 1), round(on.get("X_lid_ohm", float("nan")), 1)), "beta/k0 %s/%s" % (round(o.get("beta_over_k0_lid", float("nan")), 3), round(on.get("beta_over_k0_lid", float("nan")), 3)), flush=True)
    Path("results/sidequest_hybrid/k2_stack_contrast.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
