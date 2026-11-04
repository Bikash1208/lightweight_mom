"""Sensitivity of the two chosen switched cells to eps_r (+-0.3) and mesh. Quick job: Bloch Y_xx at a fixed kt (OFF/OFF and ON/OFF), alpha from the cell's own mean via the first-order formula.
RO3010 notch superpixel (g 0.2, window 0.4, footprint 0.4 x 0.4), 10.05 GHz;  RT5880 rim disc (diameter 6.0, rim 0.2, diode position 0.8), 10.54 GHz."""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_cell_search as C
from src.frequency import free_space_wavenumber
from src.materials import Substrate

CASES = {
    "ro3010": (dict(kind="notch", g=0.2, window=0.4, wd=0.4), 1.024 * 0 + 1.071),  # beta/k0 used for kt (first-pass value of the cell's own mean)
    "rt5880": (dict(diam=6.0, ls=2.8, pos=0.8), 1.079),
}


def job(a):
    warnings.simplefilter("ignore")
    stack, mesh, eps = a
    kw, bk = CASES[stack]
    kw = dict(kw)
    kind = kw.pop("kind", "disc")
    build = C.superpixel_notch if kind == "notch" else (lambda s_, sx, sy, sub=None, **k: C.cell(s_, sx, sy, sub=sub, **k))
    sb, f0 = C.STACKS[stack]
    if eps is not None:
        sb = Substrate(sb.name, eps, sb.tan_delta, sb.thickness)
    k0 = free_space_wavenumber(f0)
    ys = {}
    for sx in ("OFF", "ON"):
        y = C.bloch_loaded(build(stack, sx, "OFF", sub=sb, **kw), f0, mesh, bk * k0, 0.0)
        ys[sx] = complex(y[0, 0])
    d = ys["ON"] - ys["OFF"]
    iso = 0.5 * (ys["ON"] + ys["OFF"])
    try:
        b, ab, aa = C.alpha_pair(iso - d / 2, iso + d / 2, sb, f0)
    except Exception:  # noqa: BLE001
        b, ab, aa = None, None, None
    return {"stack": stack, "mesh": mesh, "eps": eps or C.STACKS[stack][0].eps_r, "Y_off_mS": [ys["OFF"].real * 1e3, ys["OFF"].imag * 1e3], "Y_on_mS": [ys["ON"].real * 1e3, ys["ON"].imag * 1e3],
            "dY_mS": [d.real * 1e3, d.imag * 1e3], "beta_k0": b and b / k0, "alpha_beam": ab, "alpha_abs": aa}


if __name__ == "__main__":
    jobs = [("ro3010", 64, None), ("ro3010", 64, 10.9), ("ro3010", 64, 11.5), ("rt5880", 64, None), ("rt5880", 64, 1.9), ("rt5880", 64, 2.5), ("rt5880", 32, None), ("ro3010", 48, None),
            ("ro3010", 96, None), ("rt5880", 96, None)]
    out = []
    with ProcessPoolExecutor(3) as ex:
        for r in ex.map(job, jobs):
            out.append(r)
            print(r, flush=True)
            Path("results/sidequest_hybrid/cell_search_sens.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
