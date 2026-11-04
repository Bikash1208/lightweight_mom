import json, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_gate_cell as gc
import numpy as np

MESH = 32
WD = 1.2e-3
YREF = 1j * 7.82e-3


def job(p):
    e, s, h, g, lt = p
    out = {"e_mm": e * 1e3, "s_mm": s * 1e3, "h_mm": h * 1e3, "g_mm": g * 1e3, "lt_mm": lt * 1e3}
    for st in ("OFF", "ON"):
        try:
            y, yx = gc.y_xx(gc.make_dogbone(e, s, h, g, st, lt=lt, wd=WD), MESH)
            out[st] = {"y_xx_ms": [yx.real * 1e3, yx.imag * 1e3], "y_yy_ms": [complex(y[1, 1]).real * 1e3, complex(y[1, 1]).imag * 1e3]}
            r = gc.strip_transmission(YREF, yx, 3)
            out[st]["strip3"] = r
        except Exception as exc:  # noqa: BLE001
            out[st] = {"error": f"{type(exc).__name__}: {exc}"[:120]}
    return out


if __name__ == "__main__":
    import warnings; warnings.simplefilter("ignore")
    grid = [(e * 1e-3, s * 1e-3, h * 1e-3, 0.4e-3, lt * 1e-3) for e in (1.2, 1.6, 2.0, 2.4) for s in (2.0, 2.8, 3.6, 4.4) for h in (5.6, 6.0) for lt in (6.2,)]
    t = time.time()
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, grid))
    Path("results/sidequest_hybrid/gate_screen5_mesh32.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    print(round(time.time() - t), "s")
    def key(r):
        v = r["ON"].get("strip3") if "error" not in r["ON"] else None
        o = r["OFF"].get("strip3") if "error" not in r["OFF"] else None
        return (v["S21_db"] if v else 0.0) + (3 * o["S21_db"] if o else -99)  # reward ON isolation, penalise OFF loss
    for r in sorted(res, key=key)[:14]:
        o, n = r["OFF"], r["ON"]
        f = lambda d: ("ERR" if "error" in d else f"y={d['y_xx_ms'][0]:.1f}{d['y_xx_ms'][1]:+.1f}j S21={d['strip3']['S21_db'] if d['strip3'] else None}")
        print(r["e_mm"], r["s_mm"], r["h_mm"], r["lt_mm"], "| OFF", f(o), "| ON", f(n))
