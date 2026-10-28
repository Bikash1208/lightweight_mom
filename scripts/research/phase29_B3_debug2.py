"""Phase 29 B3 debugging, cause of the rung-2 SLL gap (MODEL_ONLY): (T1) Frey's own N = 2 analytic case (their Section III: feeds at +-0.75 lambda0 = 1.5 lambda0 apart, beams +-5 deg in elevation, offset of feed 2 swept,
published: angle-dependent 28.3 dBi / SLL 18.5 dB / polarisation purity ~25 dB vs averaging 27.7 dBi / ~15 dB / 16 dB) with both syntheses; (T2) N = 4 tilt scan; (T3) N = 2 with the SAME beam on both feeds (only the feed
positions conflict); (T4) N = 1. Aperture: hexagon of the area of their 12.3 lambda0 square, Phase 26 design rule per feed. Writes results/phase29/B3debug2_summary.json."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/research"))
LAM = 299792458.0 / 10.54e9 * 1e3
R = 7.62 * LAM
ATTR = "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01UaAn25MEtRcvrHeagn5k3z"


def task(spec):
    import phase29_B3_debug as dbg

    from src.hologram import shared_aperture_frey as fr
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram.tensor_leaky_cp import march_rays

    k0 = 2 * np.pi * 10.54e9 / 299792458.0
    mf = lambda j, f, b: march_rays(dbg._S["table"], sa.hexagon_m(R), f, k0, b[0], b[1], 0.9, taper_c=0.5, n_psi=240)  # noqa: E731
    try:
        d = fr.design_frey(
            dbg._S["table"],
            R,
            spec["ports"],
            spec["off"],
            1000.0,
            1.0,
            mode=spec["mode"],
            n_iter=2,
            n_psi=240,
            march_fn=mf,
        )
        th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.5)
        rows = [
            dbg.field_metrics(
                J, d["x_m"], d["y_m"], (spec["ports"][j][2], spec["ports"][j][3]), d["eta_conv"][j], th, ph
            )
            for j, J in enumerate(d["J"])
        ]
        return {
            "name": spec["name"],
            "orient_err_rad": d["mean_abs_orientation_error_rad_history"][-1],
            "gain": [round(r["gain_launch_bare_dbi"], 1) for r in rows],
            "sll_15": [round(r["sll_15deg_cone_db"], 1) for r in rows],
            "sll_1p5": [round(r["sll_1p5hpbw_db"], 1) for r in rows],
            "xpd": [round(r["xpd_peak_db"], 1) for r in rows],
            "ar": [round(r["ar_peak_db"], 1) for r in rows],
        }
    except Exception as exc:  # recorded
        return {"name": spec["name"], "error": repr(exc)[:200]}


def main():
    import phase29_B3_debug as dbg
    from phase29_pointB_sweep import n_workers

    h = 0.75 * LAM
    two = [(h, 0.0, 5.0, 0.0), (-h, 0.0, 5.0, 180.0)]
    same = [(h, 0.0, 5.0, 0.0), (-h, 0.0, 5.0, 0.0)]
    specs = []
    for mode in ("frey", "aia"):
        for psi in range(0, 360, 30):
            specs.append({"name": f"T1_{mode}_psi2={psi}", "mode": mode, "ports": two, "off": [0.0, float(psi)]})
        specs.append({"name": f"T3_same_beam_{mode}", "mode": mode, "ports": same, "off": [0.0, 0.0]})
        for tilt in (0.0, 1.25, 2.5, 5.0):
            ring = 1.35 * LAM / (2 * np.sin(np.pi / 4))
            p4 = [
                (ring * np.cos(np.radians(az + 180)), ring * np.sin(np.radians(az + 180)), tilt, az)
                for az in (0.0, 180.0, 90.0, 270.0)
            ]
            specs.append({"name": f"T2_N4_{mode}_tilt={tilt}", "mode": mode, "ports": p4, "off": [0.0] * 4})
    specs.append({"name": "T4_N1", "mode": "frey", "ports": [(0.0, 0.0, 0.0, 0.0)], "off": [0.0]})
    with mp.get_context("spawn").Pool(
        min(len(specs), n_workers(per_worker_gb=2.0, reserve_gb=8.0, cap=30)), initializer=dbg._init
    ) as pool:
        res = pool.map(task, specs)
    S = {r["name"]: r for r in res}
    (ROOT / "results/phase29/B3debug2_summary.json").write_text(json.dumps(S, indent=1))
    for a in (
        ["git", "add", "-A"],
        [
            "git",
            "commit",
            "-q",
            "-m",
            f"Phase 29 B3 debugging: Frey N=2 analytic case, tilt scan, same-beam control (MODEL_ONLY)\n\n{ATTR}",
        ],
        ["git", "push", "hub", "phase29-seven-hex"],
    ):
        subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
    (ROOT / "results/phase29/B3debug2_DONE").write_text("done\n")


if __name__ == "__main__":
    main()
