"""Re-scaled switched-cell scan: stack x disc radius r (ls = r - 0.2, rim 0.2), OFF = open slit, ON = 2.44 ohm R, lattice 6.4, mesh 64, 10.5 GHz; then X (lid) OFF/ON and direct A0.  MODEL_ONLY.
    python sidequest_family_scan.py [workers] -> results/sidequest_hybrid/k2_family_scan.json"""
import json, math, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); warnings.simplefilter("ignore")
import sidequest_gate_cell as gc
import sidequest_k2_default_lib as L
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Substrate
F = 10.5e9
STK = {"RT5880_1.575": (2.2, 0.0009, 1.575), "RT5880_3.175": (2.2, 0.0009, 3.175), "RO3006_1.27": (6.15, 0.002, 1.27)}
RS = (1.6, 2.0, 2.4, 2.7, 3.0, 3.1)
FINE = {"RT5880_1.575": (3.12, 3.15, 3.17), "RT5880_3.175": (3.1, 3.15, 3.17), "RO3006_1.27": (2.2, 2.3, 2.35, 2.45, 2.5, 2.6)}  # refined scan (second pass)


def job(a):
    warnings.simplefilter("ignore")
    n, r, mode = a; e, t, h = STK[n]
    gc.SUB = Substrate(n, e, t, h * 1e-3); L.R = r
    y = sheet_tensor(L.cell(round(r - 0.2, 3), min(0.8, r - 0.4), mode), F, 64)
    return {"stack": n, "r": r, "mode": mode, "Y_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3]}


if __name__ == "__main__":
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    fine = len(sys.argv) > 2 and sys.argv[2] == "fine"
    jobs = [(n, r, m) for n in STK for r in (FINE[n] if fine else RS) for m in ("open", "ON")]
    with ProcessPoolExecutor(nw) as ex:
        res = list(ex.map(job, jobs))
    Path("results/sidequest_hybrid/k2_family_scan_raw%s.json" % ("_fine" if fine else "")).write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    print(len(res), "done")
