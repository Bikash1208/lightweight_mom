"""Phase 20 R3 item 3: bandwidth over which the four-state ladder keeps every step within 90 +- 20 deg (nominal diode, theta = 0, search mesh, state
assignment fixed at 10 GHz), plus amplitude >= -1.5 dB. MODEL_ONLY.

    python scripts/research/phase20_r3_band.py <cp2|lp3> <madp|smp> <result.json> <out.json>
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_finalist import build, channels
from phase20_r3_opt import lattice_fit, mesh_for, solve_s11, wrap180

from src.ris.cells import madp_kit, smp1345_kit

FREQS = np.arange(8.5e9, 11.51e9, 0.25e9)


def main(family, kit_name, src, out):
    d = json.loads(Path(src).read_text())
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    if family == "cp2":
        cell = build("cp2", d["params_mm"], d["period_mm"])
    else:
        from phase20_r3_lp_report import build as lp_build

        cell = lp_build("lp3", d["x_mm"], d["names"], d["codes"])
    n = mesh_for(cell.period)
    with Pool(8) as pool:
        res = pool.map(solve_s11, [(cell, st, kit, f, 0.0, n) for f in FREQS for st in cell.states()])
    co = np.array([channels(s, family)[0] for s in res]).reshape(len(FREQS), 4)
    i0 = int(np.argmin(abs(FREQS - 10e9)))
    _, perm = lattice_fit(np.degrees(np.angle(co[i0])))
    order = np.argsort(perm)
    rows, ok = [], []
    for f, z in zip(FREQS, co, strict=True):
        p = np.degrees(np.angle(z))[order]
        steps = [wrap180(p[(k + 1) % 4] - p[k]) for k in range(4)]
        err = min(max(abs(s - 90) for s in steps), max(abs(s + 90) for s in steps))
        amin = float(20 * np.log10(abs(z).min()))
        ok.append(err <= 20 and amin >= -1.5)
        rows.append({"f_ghz": float(f / 1e9), "max_step_error_deg": float(err), "min_amp_db": amin})
    lo = hi = i0
    if ok[i0]:
        while lo > 0 and ok[lo - 1]:
            lo -= 1
        while hi < len(FREQS) - 1 and ok[hi + 1]:
            hi += 1
    band = float((FREQS[hi] - FREQS[lo]) / 1e9) if ok[i0] else 0.0
    Path(out).write_text(json.dumps({"rows": rows, "band_ghz_around_10": band, "mesh": n}, indent=1))
    print(json.dumps({"band_ghz": band, "rows": rows}))


if __name__ == "__main__":
    main(*sys.argv[1:5])
