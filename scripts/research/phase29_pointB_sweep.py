"""Phase 29 design point B: parallel design sweep (BLAS 1 thread per worker, memory-guarded pool) over aperture size, ring radius, port-weight
exponent gamma and conversion strength q0; per-port metrics. Writes results/phase29/pointB_sweep.json."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import ctypes  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results/phase29/pointB_sweep.json"


def free_mem_gb():
    class MS(ctypes.Structure):
        _fields_ = [("l", ctypes.c_ulong), ("m", ctypes.c_ulong)] + [
            (n, ctypes.c_ulonglong) for n in ("tp", "ap", "tpf", "apf", "tv", "av", "ae")
        ]

    m = MS()
    m.l = ctypes.sizeof(MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ap / 2**30


def n_workers(per_worker_gb=1.5, reserve_gb=8.0, cap=None):
    return int(max(1, min(os.cpu_count() - 2 if cap is None else cap, (free_mem_gb() - reserve_gb) // per_worker_gb)))


_STATE = {}


def _init():
    from src.hologram import seven_hex_cluster as c
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    _STATE["table"] = SupercellTable.from_files(t.TABLE_GLOB)
    _STATE["green"] = c.green_function()


def port_metrics(F, theta, phi, beams):
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    rows = []
    for j, f in enumerate(F):
        g = sa.gain_db(*f)
        lb = c.lobes(g, theta, phi, within_db=3.0)
        pk = (lb[0][1], lb[0][2])
        hp = _hpbw(g, theta, phi)
        tg, pg = np.meshgrid(theta, phi, indexing="ij")
        ang = np.degrees(np.arccos(np.clip(np.sum(c._unit(tg, pg) * c._unit(*pk), axis=-1), -1, 1)))
        out = ang > 1.5 * hp
        sll = float(g[out].max() - g.max()) if out.any() else float("nan")
        rows.append(
            {
                "port": j + 1,
                "peak_gain_dbi": lb[0][0],
                "peak_theta_phi": [float(pk[0]), float(pk[1])],
                "pointing_error_deg": c.angle_between(pk, (np.degrees(beams[j][0]), np.degrees(beams[j][1]))),
                "hpbw_deg": hp,
                "sll_db": sll,
            }
        )
    return rows


def _hpbw(g_db, theta, phi):
    tg, _ = np.meshgrid(theta, phi, indexing="ij")
    w = np.sin(np.radians(np.maximum(tg, 0.25))) * np.radians(theta[1] - theta[0]) * np.radians(phi[1] - phi[0])
    return float(np.degrees(2 * np.sqrt(np.sum(w[g_db >= g_db.max() - 3.0]) / np.pi)))


def run_one(args):
    from src.hologram import shared_aperture_multiport as sa

    r, ring, gamma, w0, q0 = args
    t0 = time.time()
    try:
        d = sa.design_shared(_STATE["table"], r, ring, w0, q0, gamma=gamma)
        th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
        F = sa.port_fields(d, _STATE["green"], th, ph)
        rows = port_metrics(F, th, ph, d["beams"])
        scan = sa.scan_summary(sa.scan_table(F, th, ph))
        return {
            "r_mm": r,
            "ring_mm": ring,
            "gamma": gamma,
            "w0_mm": w0,
            "q0": q0,
            "n_cells": int(d["x_m"].size),
            "slots": np.bincount(d["slot"], minlength=4).tolist(),
            "eta_conv": [float(e) for e in d["eta_conv"]],
            "fit_error": d["fit_error_history"][-1],
            "ports": rows,
            "scan": scan,
            "seconds": time.time() - t0,
        }
    except Exception as e:  # a failed candidate is reported, never fatal
        return {"r_mm": r, "ring_mm": ring, "gamma": gamma, "w0_mm": w0, "q0": q0, "error": repr(e)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", type=float, nargs="+", default=[40, 50, 60, 80])
    ap.add_argument("--ring", type=float, nargs="+", default=[16, 22, 28])
    ap.add_argument("--gamma", type=float, nargs="+", default=[0, 1, 4])
    ap.add_argument("--w0", type=float, nargs="+", default=[1000.0])
    ap.add_argument("--q0", type=float, nargs="+", default=[20000.0])
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    jobs = [j for j in itertools.product(a.r, a.ring, a.gamma, a.w0, a.q0) if j[1] + 6.5 <= j[0] * 0.866 + 4]
    nw = n_workers()
    print(f"{len(jobs)} candidates, {nw} workers, free mem {free_mem_gb():.0f} GB", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(nw, initializer=_init) as pool:
        res = []
        for r in pool.imap_unordered(run_one, jobs):
            res.append(r)
            print(f"[{time.time() - t0:4.0f}s] {len(res)}/{len(jobs)}", flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(
        json.dumps(
            {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "INTER_PORT_COUPLING_NOT_MODELLED"], "results": res},
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
