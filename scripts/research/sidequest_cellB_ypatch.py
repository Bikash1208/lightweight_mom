"""Item 2 stage 1: patch-grid Bloch sheet admittance Y_uu(|k|/k0, psi) of the square-patch mushroom top (no via), several gaps / laminates / meshes. MODEL_ONLY.
    python sidequest_cellB_ypatch.py explore   -> results/sidequest_hybrid/cellB_ypatch_explore.json   (mesh 48)"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_cellB_core as K

F = 10.5e9
K0 = 2 * np.pi * F / K.C0
KR = [0.4226, 1.04, 1.10, 1.20, 1.30, 1.40, 1.50]  # 0.4226 = sin 25 deg (radiating harmonic), rest SW


def job(a):
    lam, h, gap, mesh, psi, kr, p = a
    er, tand = K.LAMINATES[lam]
    y = K.patch_Y((p, gap, er, tand, h, F, mesh, kr * K0, psi))
    c, s = np.cos(np.radians(psi)), np.sin(np.radians(psi))
    u = np.array([c, s])
    yuu = complex(u @ y @ u)
    yuv = complex(u @ y @ np.array([-s, c]))
    return {"lam": lam, "h": h, "gap": gap, "mesh": mesh, "psi": psi, "kr": kr, "p": p, "Yuu_S": [yuu.real, yuu.imag], "Yuv_S": [yuv.real, yuv.imag], "Y_full_S": [[[v.real, v.imag] for v in row] for row in y]}


def run(jobs, out, workers=24):
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(job, jobs))
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(res, indent=0), encoding="utf-8", newline="\n")
    print("done %d solves in %.0f s -> %s" % (len(res), time.time() - t0, out), flush=True)
    return res


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "explore":
        stacks = [("RO3006", 1.28), ("RO3003", 1.52), ("RO4003C", 1.524)]
        jobs = [(lam, h, g, 48, psi, kr, 6.4) for lam, h in stacks for g in (0.8, 1.6, 2.4, 3.2) for psi in (0.0, 45.0) for kr in KR]
        run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_explore.json")
    if what == "explore2":
        stacks = [("RO3006", 1.28), ("RO3003", 1.52), ("RO4003C", 1.524)]
        jobs = [(lam, h, g, 64, 0.0, kr, 6.4) for lam, h in stacks for g in (0.4, 0.8, 1.2, 1.6, 2.0) for kr in KR]
        run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_explore2.json")
    if what == "final":
        lam, h = "RO4003C", 1.524
        g0 = 0.8
        j32 = [(lam, h, g0, 32, psi, kr, 6.4) for psi in (0.0, 45.0, 90.0) for kr in KR]
        j64 = [(lam, h, g0, 64, psi, kr, 6.4) for psi in (0.0, 45.0, 90.0) for kr in KR]
        run(j32 + j64, ROOT / "results/sidequest_hybrid/cellB_ypatch_final_m32_m64.json", workers=24)
        gd = 6.4 / 48  # one pixel of the 96 mesh = 0.0667 mm per edge -> gap +-0.1333
        j96 = [(lam, h, g0, 96, psi, kr, 6.4) for psi in (0.0, 45.0) for kr in KR] + [(lam, h, g0 + s * gd, 96, 0.0, kr, 6.4) for s in (-1, 1) for kr in KR]
        run(j96, ROOT / "results/sidequest_hybrid/cellB_ypatch_final_m96.json", workers=14)
    if what == "er":
        jobs = [(lam, 1.524, 0.8, 64, 0.0, kr, 6.4) for lam in ("RO4003C_lo", "RO4003C_hi") for kr in KR]
        run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_er.json", workers=8)
    if what == "gap64":
        jobs = [("RO4003C", 1.524, g, 64, 0.0, kr, 6.4) for g in (0.6, 1.0) for kr in KR]
        run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_gap_m64.json", workers=8)
    if what == "m96b":
        jobs = [("RO4003C", 1.524, 0.8, 96, 0.0, kr, 6.4) for kr in KR]
        run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_final_m96.json", workers=6)
    if what == "p45":
        jobs = [("RO4003C", 1.524, g, 60, 0.0, kr, 4.5) for g in (0.3, 0.45, 0.6, 0.9) for kr in KR]
        run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_p45_m60.json", workers=8)
    if what == "d2":
        lam, h, g0 = "RO4003C", 1.524, 0.64
        gd = 6.4 / 60  # one px of the 60 mesh = 0.1067 mm per edge
        jobs = [(lam, h, g0, 40, psi, kr, 6.4) for psi in (0.0, 45.0) for kr in KR]
        jobs += [(lam, h, g0, 60, psi, kr, 6.4) for psi in (0.0, 45.0) for kr in KR]
        jobs += [(lam, h, g0 + s * 2 * gd, 60, 0.0, kr, 6.4) for s in (-1, 1) for kr in KR]
        jobs += [(l2, h, g0, 60, 0.0, kr, 6.4) for l2 in ("RO4003C_lo", "RO4003C_hi") for kr in KR]
        jobs += [(lam, h, g0, 80, psi, kr, 6.4) for psi in (0.0, 45.0) for kr in KR]
        run(jobs, ROOT / "results/sidequest_hybrid/cellB_ypatch_d2.json", workers=8)
