"""Phase 29 design point B3 (MODEL_ONLY): 7-port shared 150 mm hexagon, angle-dependent hologram superposition after Frey et al. (IEEE OJAP 2023), GA over feed
positions / hologram offsets / defocus / conversion strength, scan tests for two cones. Parallel (BLAS 1 thread per worker, memory-guarded pool).

    python scripts/research/phase29_pointB3.py ga    --cone a|b [--pop 64 --gens 25]   # writes results/phase29/pointB3_<cone>_ga.json
    python scripts/research/phase29_pointB3.py final --cone a|b                       # writes results/phase29/pointB3_<cone>.json + weights CSV

Cone a: +-20 deg, beams ~13 deg HPBW (native 150 mm), outer tilt 13 deg.  Cone b: +-25 deg, beams ~17 deg (broadened by the objective defocus), outer tilt 17 deg.
Genes (23): ring radius (1.35-1.5 lambda0), ring rotation, 6 radial jitters, 7 hologram offsets, 7 defocus phases, log10 q0 (conversion strength).
Fitness (dB-like): mean port peak gain - 0.5 max(0, worst pointing error in the cone - 2 deg) - 1.0 max(0, gain dip - 2.5 dB) - 0.15 mean_ports max(0, |HPBW - target| - 2),
maximum-ratio weights on <= 3 ports, coarse grid 2 x 3 deg."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/research"))
LAM_MM = 299792458.0 / 10.54e9 * 1e3
R_MM, W0 = 150.0, 1000.0
CONES = {
    "a": {"tilt": 13.0, "hpbw": 13.0, "cone": 20.0, "defocus_max": 7.0, "ceiling_dbi": 23.7},
    "b": {"tilt": 17.0, "hpbw": 17.0, "cone": 25.0, "defocus_max": 14.0, "ceiling_dbi": 21.7},
}
LO0 = np.array([1.35 * LAM_MM, -30.0] + [-4.0] * 6 + [0.0] * 7 + [0.0] * 7 + [1.5])
HI0 = np.array([1.50 * LAM_MM, 30.0] + [4.0] * 6 + [360.0] * 7 + [3.0] * 7 + [4.0])
_S = {}


def bounds(cone):
    hi = HI0.copy()
    hi[15:22] = CONES[cone]["defocus_max"]
    return LO0, hi


def inner_stats(rows):
    """Pointing / gain-dip statistics over ALL rows (the cone): dip = best peak gain of the table minus this beam's."""
    gmax = max(r["peak_gain_dbi"] for r in rows)
    return {
        "max_pointing_error_deg": max(r["pointing_error_deg"] for r in rows),
        "max_gain_dip_db": max(gmax - r["peak_gain_dbi"] for r in rows),
        "min_gain_dbi": min(r["peak_gain_dbi"] for r in rows),
        "gain_ref_dbi": gmax,
    }


def _init():
    from src.hologram import seven_hex_cluster as c
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    _S["table"], _S["green"] = SupercellTable.from_files(t.TABLE_GLOB), c.green_function()


def design(g, cone, mode="frey", n_psi=240, n_iter=2):
    from src.hologram import shared_aperture_frey as fr

    ports = fr.feed_ports(g[0], g[1], g[2:8], tilt_deg=CONES[cone]["tilt"])
    return fr.design_frey(
        _S["table"],
        R_MM,
        ports,
        g[8:15],
        W0,
        10.0 ** g[22],
        mode=mode,
        n_iter=n_iter,
        n_psi=n_psi,
        defocus_rad=g[15:22],
    )


def evaluate(g, cone, grid=(2.0, 3.0), n_psi=180, full=False):
    from phase29_pointB_sweep import port_metrics

    from src.hologram import shared_aperture_multiport as sa

    cfg = CONES[cone]
    d = design(g, cone, n_psi=n_psi)
    th, ph = np.arange(0.0, 90.0 + 1e-9, grid[0]), np.arange(0.0, 360.0, grid[1])
    F = sa.port_fields(d, _S["green"], th, ph)
    pm = port_metrics(F, th, ph, d["beams"])
    i = inner_stats(sa.scan_table(F, th, ph, sa.scan_directions(cfg["cone"], cfg["cone"] / 5, 30.0)))
    hp = np.array([p["hpbw_deg"] for p in pm])
    score = (
        np.mean([p["peak_gain_dbi"] for p in pm])
        - 0.5 * max(0.0, i["max_pointing_error_deg"] - 2.0)
        - 1.0 * max(0.0, i["max_gain_dip_db"] - 2.5)
        - 0.15 * np.mean(np.maximum(0.0, np.abs(hp - cfg["hpbw"]) - 2.0))
    )
    out = {
        "score": float(score),
        "pointing_max": i["max_pointing_error_deg"],
        "dip": i["max_gain_dip_db"],
        "mean_gain": float(np.mean([p["peak_gain_dbi"] for p in pm])),
        "hpbw": hp.tolist(),
    }
    if full:
        out.update({"design": d, "fields": F, "ports_metrics": pm})
    return out


def _fit(args):
    g, cone = args
    try:
        return evaluate(np.asarray(g), cone)["score"]
    except Exception:  # an infeasible individual is simply bad
        return -99.0


def run_ga(cone, pop_n, gens, minutes, seed=1, warm=False):
    from phase29_pointB_sweep import n_workers

    lo, hi = bounds(cone)
    rng = np.random.default_rng(seed)
    pop = lo + rng.random((pop_n, lo.size)) * (hi - lo)
    pop[0] = np.concatenate([[40.0, 0.0], np.zeros(6), np.zeros(7), np.full(7, 3.0 if cone == "a" else 7.0), [3.0]])
    if warm:  # second design pass: start from the first pass's best genes (and mutated copies of them)
        g0 = np.array(json.loads((ROOT / f"results/phase29/pointB3_{cone}_ga.json").read_text())["best_genes"])
        pop[1] = g0
        for k in range(2, 12):
            pop[k] = np.clip(g0 + rng.normal(size=lo.size) * 0.05 * (hi - lo), lo, hi)
    nw = n_workers(per_worker_gb=1.2, reserve_gb=8.0)
    print(f"cone {cone}: pop {pop_n}, gens {gens}, {nw} workers", flush=True)
    t0, hist = time.time(), []
    with mp.get_context("spawn").Pool(nw, initializer=_init) as pool:
        score = np.array(pool.map(_fit, [(p, cone) for p in pop]))
        for gen in range(gens):
            order = np.argsort(-score)
            hist.append(
                {"gen": gen, "best": float(score[order[0]]), "mean": float(np.mean(score)), "t_s": time.time() - t0}
            )
            print(
                f"gen {gen}: best {score[order[0]]:.2f} mean {np.mean(score):.2f} [{time.time() - t0:.0f}s]", flush=True
            )
            if (time.time() - t0) / 60 > minutes:
                break
            new = [pop[order[0]], pop[order[1]]]
            while len(new) < pop_n:
                i, j = (min(rng.integers(0, pop_n, 3), key=lambda k: -score[k]) for _ in range(2))
                child = pop[i] + (rng.random(lo.size) * 1.6 - 0.3) * (pop[j] - pop[i])
                child = child + (rng.random(lo.size) < 0.12) * rng.normal(size=lo.size) * 0.1 * (hi - lo)
                new.append(np.clip(child, lo, hi))
            pop = np.array(new)
            score = np.concatenate([score[order[:2]], np.array(pool.map(_fit, [(p, cone) for p in pop[2:]]))])
    best = pop[int(np.argmax(score))]
    (ROOT / f"results/phase29/pointB3_{cone}_ga.json").write_text(
        json.dumps(
            {
                "cone": cone,
                "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"],
                "best_genes": best.tolist(),
                "best_score": float(score.max()),
                "history": hist,
                "gene_layout": "ring_mm, rot_deg, 6 radial jitters mm, 7 offsets deg, 7 defocus rad, log10 q0",
            },
            indent=1,
        )
    )


def crossovers(F, th, ph, pm):
    """Level (dB below the own peak) of each single-port pattern at the spherical midpoint of the peak directions of neighbouring ports
    (centre-outer x6, outer-outer x6); returns lists and means."""
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    def mid(a, b):
        v = c._unit(*a) + c._unit(*b)
        v = v / np.linalg.norm(v)
        return float(np.degrees(np.arccos(np.clip(v[2], -1, 1)))), float(np.degrees(np.arctan2(v[1], v[0])) % 360.0)

    G = [sa.gain_db(*f) for f in F]
    pk = [tuple(p["peak_theta_phi"]) for p in pm]
    out = {}
    for name, pairs in (
        ("centre_outer", [(0, k) for k in range(1, 7)]),
        ("outer_outer", [(k, k % 6 + 1) for k in range(1, 7)]),
    ):
        lv = []
        for a, b in pairs:
            m = mid(pk[a], pk[b])
            it, ip = int(round(m[0] / (th[1] - th[0]))), int(round(m[1] / (ph[1] - ph[0]))) % ph.size
            lv += [float(G[a][it, ip] - G[a].max()), float(G[b][it, ip] - G[b].max())]
        out[name] = {"mean_db": float(np.mean(lv)), "min_db": float(np.min(lv)), "max_db": float(np.max(lv))}
    return out


def sum_diff(F, th, ph):
    from phase29_pointB_sweep import _hpbw

    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    az = np.radians([60.0 * k for k in range(6)])
    cases = {
        "sum": np.array([1.0] * 7),
        "delta_x": np.concatenate([[0.0], np.cos(az)]),
        "delta_y": np.concatenate([[0.0], np.sin(az)]),
    }
    res = {}
    gs = {}
    for name, w in cases.items():
        w = w.astype(complex) / np.linalg.norm(w)
        gs[name] = sa.gain_db(*sa.combine(F, w))
    pk = np.unravel_index(np.argmax(gs["sum"]), gs["sum"].shape)
    lb = c.lobes(gs["sum"], th, ph)
    tg, pg = np.meshgrid(th, ph, indexing="ij")
    hp = _hpbw(gs["sum"], th, ph)
    ang = np.degrees(np.arccos(np.clip(np.sum(c._unit(tg, pg) * c._unit(th[pk[0]], ph[pk[1]]), axis=-1), -1, 1)))
    res["sum"] = {
        "peak_gain_dbi": float(gs["sum"].max()),
        "peak_theta_phi_deg": [float(lb[0][1]), float(lb[0][2])],
        "hpbw_deg": hp,
        "sll_db": float(gs["sum"][ang > 1.5 * hp].max() - gs["sum"].max()),
    }
    for name in ("delta_x", "delta_y"):
        g = gs[name]
        res[name] = {
            "peak_gain_dbi": float(g.max()),
            "null_depth_at_sum_peak_db": float(g[pk] - g.max()),
            "gain_at_broadside_rel_sum_peak_db": float(g[0, 0] - gs["sum"].max()),
        }
    return res


def final(cone):
    import phase29_compact7_tune as ct
    import phase29_pointB_stage1 as st1
    from phase29_pointB_sweep import n_workers

    from src.hologram import shared_aperture_multiport as sa

    cfg = CONES[cone]
    gj = json.loads((ROOT / f"results/phase29/pointB3_{cone}_ga.json").read_text())
    g = np.array(gj["best_genes"])
    _init()
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    th2, ph2 = np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 2.0)
    o = evaluate(g, cone, grid=(1.0, 1.0), n_psi=240, full=True)
    d, F, pm = o["design"], o["fields"], o["ports_metrics"]
    step = 2.0 if cone == "a" else 2.5
    dirs = sa.scan_directions(cfg["cone"], step, 15.0)
    mrc_rows = sa.scan_table(F, th, ph, dirs)
    F2 = sa.port_fields(d, _S["green"], th2, ph2)
    with mp.get_context("spawn").Pool(
        min(40, n_workers(per_worker_gb=1.0)), initializer=ct._init, initargs=(F, th, ph, F2, th2, ph2, 1.5)
    ) as pool:
        tuned = pool.map(ct.one, dirs)
    gmax = max(r["peak_gain_dbi"] for r in tuned)
    err = np.array([r["pointing_error_deg"] for r in tuned])
    dip = gmax - np.array([r["peak_gain_dbi"] for r in tuned])
    edge = [r for r in tuned if abs(r["theta_deg"] - dirs[-1][0]) < 1e-9 and abs((r["phi_deg"] - 30.0) % 60.0) < 1e-9]
    bias = json.loads((ROOT / "results/phase26/fourhex_scope.json").read_text())["coupling"]["one_bit_board"][
        "bias_db_raw_minus_mean_hfss"
    ]
    cp = st1.coupling_db(d, bias)
    gains = [p["peak_gain_dbi"] for p in pm]
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_PORT_COUPLING_NOT_MODELLED"],
        "cone": cone,
        "config": cfg,
        "genes": g.tolist(),
        "ga_score": gj["best_score"],
        "ring_mm": float(g[0]),
        "ring_over_lambda0": float(g[0] / LAM_MM),
        "n_cells": int(d["x_m"].size),
        "slots": np.bincount(d["slot"], minlength=4).tolist(),
        "eta_conv": d["eta_conv"],
        "ports": pm,
        "crossover": crossovers(F, th, ph, pm),
        "sum_diff": sum_diff(F, th, ph),
        "scan_mrc": {
            **inner_stats(mrc_rows),
            "median_pointing_error_deg": float(np.median([r["pointing_error_deg"] for r in mrc_rows])),
            "fraction_within_2deg": float(np.mean([r["pointing_error_deg"] <= 2.0 for r in mrc_rows])),
        },
        "scan_tuned": {
            "n": len(tuned),
            "gain_ref_dbi": gmax,
            "max_pointing_error_deg": float(err.max()),
            "median_pointing_error_deg": float(np.median(err)),
            "fraction_within_2deg": float(np.mean(err <= 2.0)),
            "fraction_within_3deg": float(np.mean(err <= 3.0)),
            "max_gain_dip_db": float(dip.max()),
            "min_gain_dbi": float(min(r["peak_gain_dbi"] for r in tuned)),
            "pass_pointing_2deg": bool(err.max() <= 2.0),
            "pass_dip_3db": bool(dip.max() <= 3.0),
            "cone_edge_between_outer_beams": [
                {
                    "theta_deg": r["theta_deg"],
                    "phi_deg": r["phi_deg"],
                    "pointing_error_deg": r["pointing_error_deg"],
                    "peak_gain_dbi": r["peak_gain_dbi"],
                }
                for r in edge
            ],
            "method_counts": {m: sum(r["method"] == m for r in tuned) for m in {r["method"] for r in tuned}},
        },
        "coupling_db": cp,
        "coupling_state": "ESTIMATE_BIAS_FROM_ONE_BOARD_NOT_CALIBRATED_AT_THIS_SPACING",
        "ceiling_gap_db": {
            "ceiling_dbi_as_given": cfg["ceiling_dbi"],
            "best_port_dbi": max(gains),
            "mean_port_dbi": float(np.mean(gains)),
            "gap_best_port_db": cfg["ceiling_dbi"] - max(gains),
            "gap_mean_port_db": cfg["ceiling_dbi"] - float(np.mean(gains)),
            "gap_min_scan_db": cfg["ceiling_dbi"] - float(min(r["peak_gain_dbi"] for r in tuned)),
        },
        "tuned_rows": tuned,
        "cells": [
            {
                "x_mm": float(x * 1e3),
                "y_mm": float(y * 1e3),
                "slot_L_mm": float(a),
                "slot_W_mm": float(b),
                "slot_angle_deg": float(np.degrees(c_)),
            }
            for x, y, a, b, c_ in zip(d["x_m"], d["y_m"], d["slot_L_mm"], d["slot_W_mm"], d["axis_rad"], strict=True)
        ],
        "feeds_mm": [[p[0], p[1]] for p in d["ports"]],
    }
    (ROOT / f"results/phase29/pointB3_{cone}.json").write_text(json.dumps(res, indent=1, default=float))
    with open(ROOT / f"results/phase29/pointB3_{cone}_weights.csv", "w", newline="", encoding="utf-8") as f:
        w_ = csv.writer(f)
        w_.writerow(
            ["theta_deg", "phi_deg"]
            + [f"P{k}_{q}" for k in range(1, 8) for q in ("power_W", "phase_deg")]
            + ["method", "pred_peak_gain_dbi", "pred_pointing_error_deg"]
        )
        for r in tuned:
            w = np.array([complex(a, b) for a, b in r["weights"]])
            row = []
            for x in w:
                row += [f"{abs(x) ** 2:.5f}", f"{np.degrees(np.angle(x)) if abs(x) > 1e-9 else 0.0:.2f}"]
            w_.writerow(
                [f"{r['theta_deg']:g}", f"{r['phi_deg']:g}"]
                + row
                + [r["method"], f"{r['peak_gain_dbi']:.2f}", f"{r['pointing_error_deg']:.2f}"]
            )
    t = res["scan_tuned"]
    print(
        cone,
        "ports gain",
        [round(x, 1) for x in gains],
        "hpbw",
        [round(p["hpbw_deg"]) for p in pm],
        "sll",
        [round(p["sll_db"], 1) for p in pm],
    )
    print(
        "  tuned: max point %.1f med %.1f f2 %.2f dip %.2f minG %.1f | crossover %s | sum %s"
        % (
            t["max_pointing_error_deg"],
            t["median_pointing_error_deg"],
            t["fraction_within_2deg"],
            t["max_gain_dip_db"],
            t["min_gain_dbi"],
            {k: round(v["mean_db"], 1) for k, v in res["crossover"].items()},
            {k: (round(v, 1) if isinstance(v, float) else v) for k, v in res["sum_diff"]["sum"].items()},
        )
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["ga", "final", "eval"])
    ap.add_argument("--cone", choices=["a", "b"], required=True)
    ap.add_argument("--pop", type=int, default=64)
    ap.add_argument("--gens", type=int, default=25)
    ap.add_argument("--max-minutes", type=float, default=35.0)
    ap.add_argument("--warm", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.cmd == "ga":
        run_ga(a.cone, a.pop, a.gens, a.max_minutes, a.seed, a.warm)
    elif a.cmd == "final":
        final(a.cone)
    else:
        _init()
        t0 = time.time()
        g = np.concatenate([[40.0, 0.0], np.zeros(6), np.zeros(7), np.full(7, 3.0 if a.cone == "a" else 7.0), [3.0]])
        r = evaluate(g, a.cone)
        print(time.time() - t0, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
