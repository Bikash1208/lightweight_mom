"""Why is the MoM OFF eigenfrequency of the RO3006 r 2.2 cell 0.12 GHz below HFSS?  Bloch tensor at 90 deg (beta = pi/2 / 6.4 mm) at mesh 64 / 80 / 96, OFF (open slit) and ON (2.44 ohm), around the HFSS frequencies (OFF 10.0193, ON 10.3617 GHz).
Eigenfrequency = root of stack_b(beta, f) + Im Yxx(f) (transverse resonance with the 30 mm PEC lid, as the package prediction).  MODEL_ONLY.  -> results/sidequest_hybrid/k2_offcal.json"""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
import sidequest_cell_ro3006_package as P

HF = {"open": 10.0193, "ON": 10.3617}
FGRID = {"open": (9.80, 9.95, 10.0193, 10.10, 10.25), "ON": (10.15, 10.25, 10.3617, 10.45, 10.55)}


def job(a):
    mesh, mode, f = a
    import sidequest_cell_search as C
    import sidequest_gate_cell as gc
    import sidequest_k2_default_lib as L
    from src.materials import Substrate
    gc.SUB = Substrate("ro3006", P.ER, P.TAND, P.H * 1e-3)
    L.R = P.R
    y = C.bloch_loaded(L.cell(round(P.R - 0.2, 3), P.POS, mode), f * 1e9, mesh, P.BETA, 0.0)
    return {"mesh": mesh, "mode": mode, "f": f, "Y_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3]}


if __name__ == "__main__":
    from scipy.interpolate import CubicSpline
    from scipy.optimize import brentq
    meshes = [int(x) for x in sys.argv[1:]] or [64, 80, 96]
    jobs = [(m, mode, f) for m in meshes for mode in ("open", "ON") for f in FGRID[mode]]
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(job, jobs))
    out = []
    for m in meshes:
        for mode in ("open", "ON"):
            rows = sorted([x for x in res if x["mesh"] == m and x["mode"] == mode], key=lambda x: x["f"])
            fs = np.array([x["f"] for x in rows]) * 1e9
            b = np.array([x["Y_mS"][1] for x in rows]) * 1e-3
            cs = CubicSpline(fs, b)
            try:
                fr = brentq(lambda f, cs=cs: P.stack_b(P.BETA, f, P.ER, P.H) + float(cs(f)), fs[0], fs[-1], xtol=1e3) / 1e9
            except ValueError:
                fr = None
            bh = float(-P.stack_b(P.BETA, HF[mode] * 1e9, P.ER, P.H)) * 1e3  # cell susceptance HFSS requires [mS]
            bm = float(cs(HF[mode] * 1e9)) * 1e3
            out.append({"mesh": m, "mode": mode, "f_eigen_GHz": fr, "f_HFSS_GHz": HF[mode], "B_model_at_HFSS_f_mS": bm, "B_required_by_HFSS_mS": bh, "delta_B_mS": bh - bm})
            print("mesh %d %-4s eigen %s (HFSS %.4f)  B_model %.3f mS  B_HFSS-required %.3f mS  delta %.3f mS" % (m, mode, None if fr is None else round(fr, 4), HF[mode], bm, bh, bh - bm), flush=True)
    Path("results/sidequest_hybrid/k2_offcal.json").write_text(json.dumps({"rows": out, "raw": res}, indent=1), encoding="utf-8", newline="\n")
