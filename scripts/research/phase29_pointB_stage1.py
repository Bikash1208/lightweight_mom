"""Phase 29 design point B, stage 1 report data (MODEL_ONLY): best shared-aperture design, per-port patterns, coupling estimate, scan test with
maximum-ratio weights (all directions) and with optimised weights (parallel), and the 7-small-tile comparator. Writes results/phase29/pointB_stage1.json."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
S = {}


def _init(fields, th, ph):
    S["F"], S["th"], S["ph"] = fields, th, ph


def _opt(d):
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    F, th, ph = S["F"], S["th"], S["ph"]
    w = sa.optimise_weights(F, th, ph, d)
    lb = c.lobes(sa.gain_db(*sa.combine(F, w)), th, ph)
    return {
        "theta_deg": d[0],
        "phi_deg": d[1],
        "peak_gain_dbi": lb[0][0],
        "pointing_error_deg": c.angle_between((lb[0][1], lb[0][2]), d),
    }


def coupling_db(design, calib_bias_db):
    from src.hologram import tensor_cp_design as t

    lam_sw = 299792458.0 / t.F / 1.06
    w_eff = t.LAUNCH * lam_sw / (2 * np.pi)
    out = {}
    feeds = design["feeds_m"]
    for a in range(7):
        for b in range(a + 1, 7):
            d = float(np.hypot(feeds[a][0] - feeds[b][0], feeds[a][1] - feeds[b][1]))
            m = design["march"][a]
            psi = np.arctan2(feeds[b][1] - feeds[a][1], feeds[b][0] - feeds[a][0])
            ip = int(np.argmin(np.abs(((m["psi"] - psi + np.pi) % (2 * np.pi)) - np.pi)))
            ir = min(int(round(d / m["drho"])), m["rho"].size - 1)
            raw = t.LAUNCH / (2 * np.pi) * m["w"][ip, ir] / d * w_eff
            out[f"{a + 1}-{b + 1}"] = {
                "distance_mm": d * 1e3,
                "raw_db": float(10 * np.log10(raw)),
                "calibrated_db": float(10 * np.log10(raw) - calib_bias_db),
            }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", type=float, default=90.0)
    ap.add_argument("--ring", type=float, default=32.0)
    ap.add_argument("--gamma", type=float, default=8.0)
    ap.add_argument("--tiles", type=float, nargs="*", default=[20.0, 30.0, 40.0])
    a = ap.parse_args()
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    table = SupercellTable.from_files(t.TABLE_GLOB)
    green = c.green_function()
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    d = sa.design_shared(table, a.r, a.ring, 1000.0, 20000.0, gamma=a.gamma)
    F = sa.port_fields(d, green, th, ph)
    sys_path = str(ROOT / "scripts/research")
    import sys

    sys.path.insert(0, sys_path)
    from phase29_pointB_sweep import port_metrics

    bias = json.loads((ROOT / "results/phase26/fourhex_scope.json").read_text())["coupling"]["one_bit_board"][
        "bias_db_raw_minus_mean_hfss"
    ]
    rows = sa.scan_table(F, th, ph)
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_PORT_COUPLING_NOT_MODELLED"],
        "shared": {
            "r_mm": a.r,
            "ring_mm": a.ring,
            "gamma": a.gamma,
            "n_cells": int(d["x_m"].size),
            "slots": np.bincount(d["slot"], minlength=4).tolist(),
            "eta_conv": d["eta_conv"],
            "ports": port_metrics(F, th, ph, d["beams"]),
            "coupling_db": coupling_db(d, bias),
            "coupling_state": "ESTIMATE_BIAS_FROM_ONE_BOARD_NOT_CALIBRATED_AT_THIS_SPACING",
            "scan_mrc": sa.scan_summary(rows),
            "scan_mrc_rows": rows,
        },
    }
    dirs = sa.scan_directions()
    th2, ph2 = np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 2.0)
    F2 = sa.port_fields(d, green, th2, ph2)
    with mp.get_context("spawn").Pool(min(40, os.cpu_count() - 2), initializer=_init, initargs=(F2, th2, ph2)) as pool:
        opt = pool.map(_opt, dirs)
    gmax = max(r["peak_gain_dbi"] for r in opt)
    inner = [r for r in opt if r["theta_deg"] <= 30]
    res["shared"]["scan_optimised_weights"] = {
        "inner_le_30": {
            "max_pointing_error_deg": max(r["pointing_error_deg"] for r in inner),
            "median_pointing_error_deg": float(np.median([r["pointing_error_deg"] for r in inner])),
            "fraction_within_3deg": float(np.mean([r["pointing_error_deg"] <= 3 for r in inner])),
            "max_gain_dip_db": gmax - min(r["peak_gain_dbi"] for r in inner),
        },
        "note": "penalised maximisation, 3 active ports, 2 deg x 2 deg grid",
    }
    tiles = {}
    for R in a.tiles:
        cen = c.tile_centres_mm(R)
        own = c.owned_cells_mm(cen, R)
        Ft = []
        for j, (x, y, az) in enumerate(cen):
            Ft.append(
                c.tile_field(
                    c.design_tile(table, (x, y), R, *((0.0, 0.0) if j == 0 else (30.0, az)), cells=own[j]),
                    green,
                    th,
                    ph,
                )
            )
        rt = sa.scan_table(Ft, th, ph)
        tiles[f"R{R:g}"] = {
            "summary": sa.scan_summary(rt),
            "median_pointing_inner_deg": float(
                np.median([r["pointing_error_deg"] for r in rt if r["theta_deg"] <= 30])
            ),
            "ports": port_metrics(
                Ft, th, ph, [(0.0, 0.0)] + [(np.radians(30.0), np.radians(60.0 * k)) for k in range(6)]
            ),
        }
    res["tiles_comparator"] = tiles
    (ROOT / "results/phase29/pointB_stage1.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
