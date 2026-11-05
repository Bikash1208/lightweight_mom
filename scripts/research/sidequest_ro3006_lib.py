"""Final-cell library for the RO3006 1.27 mm board, disc r 2.22 (rim 0.2, slit 0.4), lattice 6.4, 10.5 GHz, MoM sheet tensor at mesh 96 (0.0667 mm pixel; normal incidence like the earlier library), analog statics + switched pairs.
Statics: solid disc, open slits of half-length 0.357 / 0.643 / 0.786 / 0.929 x (r - 0.2), strapped footprints (0.786, 1.0).  Switched pairs (OFF = open slit, ON = 2.44 ohm R on the 0.4 x 1.0 footprints at 0.8 mm): slit 1.0 / 0.643 / 0.786 / 0.929 x (r - 0.2).  MODEL_ONLY.
    python sidequest_ro3006_lib.py [mesh=96] [r=2.22] -> results/sidequest_hybrid/ro3006_r<r>_mesh<mesh>_lib.json"""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
F = 10.5e9
ER, TD, H = 6.15, 0.002, 1.27
MESH = int(sys.argv[1]) if len(sys.argv) > 1 else 96
R = float(sys.argv[2]) if len(sys.argv) > 2 else 2.22
FRS = (0.357, 0.643, 0.786, 0.929)


def job(a):
    warnings.simplefilter("ignore")
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.hologram.coded_cell_eval import sheet_tensor
    from src.materials import Substrate
    kind, fr = a
    gc.SUB = Substrate("ro3006", ER, TD, H * 1e-3)
    L.R = R
    lm = round(R - 0.2, 4)
    ls = round(lm * fr, 4)
    y = sheet_tensor(L.cell(ls, 0.8, kind), F, MESH)
    return {"kind": kind, "fr": fr, "ls": ls, "Y_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3], "Yyy_mS": [y[1, 1].real * 1e3, y[1, 1].imag * 1e3]}


if __name__ == "__main__":
    jobs = [("open", 1.0), ("ON", 1.0), ("none", 0.0)] + [("open", f) for f in FRS] + [("ON", f) for f in FRS[1:]] + [("strap", 0.786), ("strap", 1.0)]
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(job, jobs))
    Path(f"results/sidequest_hybrid/ro3006_r{R}_mesh{MESH}_lib.json").write_text(json.dumps({"r": R, "mesh": MESH, "stack": [ER, TD, H], "rows": res}, indent=1), encoding="utf-8", newline="\n")
    for r in res:
        print(r["kind"], r["fr"], "Y %.3f%+.3fj" % tuple(r["Y_mS"]), flush=True)
