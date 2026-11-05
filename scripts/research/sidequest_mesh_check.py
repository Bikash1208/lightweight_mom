"""Mesh convergence of the cell quantities vs disc radius (RO3006 1.27, eps_r 6.15): pair (open / ON) at mesh 64 (0.1 mm pixel) and 96 (0.0667 mm pixel), r in {2.19 .. 2.31}.  Radius steps of 0.03 mm are below the 0.1 mm pixel.  MODEL_ONLY.
    python sidequest_mesh_check.py -> results/sidequest_hybrid/k2_mesh_check.json"""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
F = 10.5e9
RS = (2.19, 2.22, 2.25, 2.28, 2.31)
H, TD, ER = 1.27, 0.002, 6.15


def job(a):
    warnings.simplefilter("ignore")
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.hologram.coded_cell_eval import sheet_tensor
    from src.materials import Substrate
    r, mode, mesh = a
    gc.SUB = Substrate("ro", ER, TD, H * 1e-3)
    L.R = r
    y = sheet_tensor(L.cell(round(r - 0.2, 3), 0.8, mode), F, mesh)
    return {"r": r, "mode": mode, "mesh": mesh, "Y_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3]}


if __name__ == "__main__":
    import sidequest_a0_direct as A
    jobs = [(r, m, n) for n in (64, 96) for r in RS for m in ("open", "ON")]
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(job, jobs))
    out = []
    for n in (64, 96):
        for r in RS:
            d = {x["mode"]: complex(*x["Y_mS"]) for x in res if x["r"] == r and x["mesh"] == n}
            mid = (d["open"].imag + d["ON"].imag) / 2
            dy = abs(d["open"] - d["ON"])
            a0, bk = A.a0_direct(mid * 1e-3, dy * 1e-3 * 3.13 / 2.63, ER, H * 1e-3, F, tand=TD)
            out.append({"mesh": n, "r": r, "Ybar_mS": mid, "dY_mS": dy, "A0_cal": float(a0), "beta_k0": float(bk)})
            print("mesh %d r %.2f Ybar %.2f dY %.2f A0cal %.2f beta/k0 %.3f" % (n, r, mid, dy, a0, bk), flush=True)
    Path("results/sidequest_hybrid/k2_mesh_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
