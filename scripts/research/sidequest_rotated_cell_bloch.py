"""Rotated one-diode disc cell (RT5880, diameter 6.0, gap 0.4, footprint 1.2): loaded Bloch tensor at kt = beta (1.08 k0) along the diode axis, chi = 0 and 45 deg, mesh 64. Output rotated_cell_bloch.json."""
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_cell_search as C
import sidequest_coded_pol_cell as PC
from src.frequency import free_space_wavenumber

F0 = 10.54e9


def job(a):
    warnings.simplefilter("ignore")
    state, chi, mesh = a
    k0 = free_space_wavenumber(F0)
    c = PC.rot_cell(state, chi, diam=6.0e-3)
    t = np.radians(chi)
    y = C.bloch_loaded(c, F0, mesh, 1.08 * k0 * np.cos(t), 1.08 * k0 * np.sin(t))  # propagation along the diode axis
    u, v = np.array([np.cos(t), np.sin(t)]), np.array([-np.sin(t), np.cos(t)])
    return {"state": state, "chi": chi, "mesh": mesh, "Y_uu_mS": [complex(u @ y @ u).real * 1e3, complex(u @ y @ u).imag * 1e3], "Y_vv_mS": [complex(v @ y @ v).real * 1e3, complex(v @ y @ v).imag * 1e3]}


if __name__ == "__main__":
    jobs = [(s, chi, 64) for chi in (0, 45) for s in ("OFF", "ON")]
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, jobs))
    out = {"rows": res}
    for chi in (0, 45):
        o = [r for r in res if r["chi"] == chi]
        yo = complex(*[r for r in o if r["state"] == "OFF"][0]["Y_uu_mS"]) * 1e-3
        yn = complex(*[r for r in o if r["state"] == "ON"][0]["Y_uu_mS"]) * 1e-3
        yv = complex(*[r for r in o if r["state"] == "OFF"][0]["Y_vv_mS"]) * 1e-3
        iso = 0.5 * (0.5 * (yo + yn) + yv)
        bk, ab, aa = C.alpha_pair(iso - (yn - yo) / 2, iso + (yn - yo) / 2, C.STACKS["rt5880"][0], F0)
        out[f"chi{chi}"] = {"dY_mS": [(yn - yo).real * 1e3, (yn - yo).imag * 1e3], "Y_iso_mS": [iso.real * 1e3, iso.imag * 1e3], "beta_k0": bk, "alpha_beam": ab, "alpha_abs": aa}
        print(chi, out[f"chi{chi}"], flush=True)
    Path("results/sidequest_hybrid/rotated_cell_bloch.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
