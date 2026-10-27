"""Phase 29 B3 sanity check: reproduce the TREND of Frey et al. (IEEE OJAP 2023) Fig. 15-17 with the repo's cell model: N = 2..8 feeds equidistant on a circle
(adjacent feeds 1.5 lambda0 apart, eq. 23-24), hologram offsets 0, pencil beams (tilt 5 deg, azimuth n 360/N: the paper gives +-5 deg elevation only for N = 2),
aperture = hexagon of the same area as their 48 x 48 mm (12.3 lambda0)^2 square (R = 7.62 lambda0 = 217 mm), strong slots, no defocus.
`frey` = angle-dependent circular-mean superposition, `aia` = pixel-wise averaging of the unit tensors (their state of the art). Compared: mean peak RHCP gain,
mean RHCP SLL (outside 1.5 HPBW), mean XPD at the peak. Their curves are plots (values not readable from the text): only the trends are compared:
gain falls with N for both methods, the angle-dependent one >= ~0.7 dB above AiA, SLL deteriorates with N, polarisation purity higher for the angle-dependent one.
Writes results/phase29/frey_trend.json. MODEL_ONLY."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/research"))
LAM = 299792458.0 / 10.54e9 * 1e3
R_MM = 7.62 * LAM
_S = {}


def _init():
    from src.hologram import seven_hex_cluster as c
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    _S["table"], _S["green"] = SupercellTable.from_files(t.TABLE_GLOB), c.green_function()


def one(args):
    n, mode = args
    from phase29_pointB_sweep import port_metrics

    from src.hologram import shared_aperture_frey as fr
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t

    ring = 1.5 * LAM / (2 * np.sin(np.pi / n))
    ports = []
    for k in range(n):
        az = np.radians(360.0 * k / n)
        ports.append((ring * np.cos(az), ring * np.sin(az), 5.0, 360.0 * k / n))
    d = fr.design_frey(_S["table"], R_MM, ports, [0.0] * n, 1000.0, 1e5, mode=mode, n_iter=2, n_psi=240)
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 2.0)
    F = sa.port_fields(d, _S["green"], th, ph)
    pm = port_metrics(F, th, ph, d["beams"])
    xpd = [t.pol_metrics(f[0], f[1], th, ph, tuple(p["peak_theta_phi"]))["xpd_db"] for f, p in zip(F, pm, strict=True)]
    return {
        "n": n,
        "mode": mode,
        "ring_mm": ring,
        "cells": int(d["x_m"].size),
        "mean_peak_gain_dbi": float(np.mean([p["peak_gain_dbi"] for p in pm])),
        "mean_sll_db": float(np.mean([p["sll_db"] for p in pm])),
        "mean_xpd_db": float(np.mean(xpd)),
        "mean_hpbw_deg": float(np.mean([p["hpbw_deg"] for p in pm])),
        "eta_conv_mean": float(np.mean(d["eta_conv"])),
    }


def main():
    from phase29_pointB_sweep import n_workers

    jobs = [(n, m) for n in range(2, 9) for m in ("frey", "aia")]
    with mp.get_context("spawn").Pool(min(len(jobs), n_workers(per_worker_gb=1.5)), initializer=_init) as pool:
        res = pool.map(one, jobs)
    by = {}
    for r in res:
        by.setdefault(r["n"], {})[r["mode"]] = r
    trend = [
        {
            "n": n,
            "gain_frey_dbi": v["frey"]["mean_peak_gain_dbi"],
            "gain_aia_dbi": v["aia"]["mean_peak_gain_dbi"],
            "gain_diff_db": v["frey"]["mean_peak_gain_dbi"] - v["aia"]["mean_peak_gain_dbi"],
            "sll_frey_db": v["frey"]["mean_sll_db"],
            "sll_aia_db": v["aia"]["mean_sll_db"],
            "xpd_frey_db": v["frey"]["mean_xpd_db"],
            "xpd_aia_db": v["aia"]["mean_xpd_db"],
        }
        for n, v in sorted(by.items())
    ]
    g = [t_["gain_frey_dbi"] for t_ in trend]
    out = {
        "labels": ["MODEL_ONLY", "TREND_ONLY_NOT_NUMERIC_COMPARISON"],
        "aperture_R_mm": R_MM,
        "trend": trend,
        "rows": res,
        "checks": {
            "gain_falls_with_n_frey": bool(g[-1] < g[0]),
            "gain_monotone_non_increasing_frey": bool(all(b <= a + 0.05 for a, b in zip(g, g[1:], strict=False))),
            "frey_minus_aia_min_db": min(t_["gain_diff_db"] for t_ in trend),
            "frey_minus_aia_max_db": max(t_["gain_diff_db"] for t_ in trend),
            "sll_deteriorates_with_n_frey": bool(trend[-1]["sll_frey_db"] > trend[0]["sll_frey_db"]),
        },
    }
    (ROOT / "results/phase29/frey_trend.json").write_text(json.dumps(out, indent=1))
    for t_ in trend:
        print({k: round(v, 2) for k, v in t_.items()})
    print(out["checks"])


if __name__ == "__main__":
    main()
