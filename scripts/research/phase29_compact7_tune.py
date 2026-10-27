"""Phase 29 compact 7-tile layout: per-direction weight tuning for R30 / R40 tiles (parallel, BLAS 1 per worker). For every direction of the table
several candidate weight sets are evaluated on the 1 deg x 1 deg RHCP pattern (MRC with 2 / 3 active ports, constrained optimisations with different
margins; never more than 3 active ports); chosen: the highest-gain candidate pointing within `--prefer` deg, else the smallest pointing error. Writes results/phase29/compact7_tune.json
(summary per tile size + the chosen weight table of the best size) and a compact `compact7_fields.npz` (port fields of the best size)."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
S = {}


def _init(F, th, ph, F2, th2, ph2, prefer):
    S.update(F=F, th=th, ph=ph, F2=F2, th2=th2, ph2=ph2, prefer=prefer)


def evaluate(F, th, ph, w, d):
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    g = sa.gain_db(*sa.combine(F, w))
    lb = c.lobes(g, th, ph, within_db=3.0)
    return lb[0][0], c.angle_between((lb[0][1], lb[0][2]), d), (lb[0][1], lb[0][2]), len(lb)


def one(d):
    from src.hologram import shared_aperture_multiport as sa

    F, th, ph, F2, th2, ph2 = (S[k] for k in ("F", "th", "ph", "F2", "th2", "ph2"))
    cands = {"mrc2": sa.mrc_weights(F, th, ph, d, 2), "mrc3": sa.mrc_weights(F, th, ph, d, 3)}
    for name, kw in (
        ("opt3_m63", dict(n_active=3)),
        ("opt3_m40", dict(n_active=3, margin=0.40)),
        ("opt2_m63", dict(n_active=2)),
        ("opt3_p2", dict(n_active=3, point_deg=2.0, inner_deg=15.0, margin=0.5)),
    ):
        cands[name] = sa.optimise_weights(F2, th2, ph2, d, **kw)
    ev = {k: evaluate(F, th, ph, w, d) for k, w in cands.items()}
    ok = [k for k, v in ev.items() if v[1] <= S["prefer"]]
    best = max(ok, key=lambda k: ev[k][0]) if ok else min(ev, key=lambda k: ev[k][1])
    g, err, pk, nl = ev[best]
    return {
        "theta_deg": d[0],
        "phi_deg": d[1],
        "method": best,
        "weights": [[float(x.real), float(x.imag)] for x in cands[best]],
        "peak_gain_dbi": float(g),
        "pointing_error_deg": float(err),
        "peak_theta_phi_deg": [float(pk[0]), float(pk[1])],
        "n_lobes_3db": nl,
        "mrc3_pointing_error_deg": float(ev["mrc3"][1]),
    }


def summary(rows):
    gmax = max(r["peak_gain_dbi"] for r in rows)
    out = {"gain_ref_dbi": gmax}
    for name, sel in (("inner_le_30", lambda r: r["theta_deg"] <= 30), ("band_30_45", lambda r: r["theta_deg"] > 30)):
        rs = [r for r in rows if sel(r)]
        e = np.array([r["pointing_error_deg"] for r in rs])
        out[name] = {
            "n": len(rs),
            "max_pointing_error_deg": float(e.max()),
            "median_pointing_error_deg": float(np.median(e)),
            "fraction_within_3deg": float(np.mean(e <= 3)),
            "fraction_within_8deg": float(np.mean(e <= 8)),
            "max_gain_dip_db": float(gmax - min(r["peak_gain_dbi"] for r in rs)),
            "min_gain_dbi": float(min(r["peak_gain_dbi"] for r in rs)),
        }
    i = out["inner_le_30"]
    out["pass_relaxed_inner"] = bool(i["max_pointing_error_deg"] <= 8.0 and i["max_gain_dip_db"] <= 3.0)
    return out


def tile_fields(table, green, R, th, ph, tilt=30.0):
    from src.hologram import seven_hex_cluster as c

    cen = c.tile_centres_mm(R)
    own = c.owned_cells_mm(cen, R)
    tiles = [
        c.design_tile(table, (x, y), R, *((0.0, 0.0) if j == 0 else (tilt, az)), cells=own[j])
        for j, (x, y, az) in enumerate(cen)
    ]
    return tiles, [c.tile_field(tl, green, th, ph) for tl in tiles]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=float, nargs="+", default=[30.0, 40.0])
    ap.add_argument("--prefer", type=float, default=3.0)
    a = ap.parse_args()
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    table, green = SupercellTable.from_files(t.TABLE_GLOB), c.green_function()
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    th2, ph2 = np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 2.0)
    dirs = sa.scan_directions()
    res, best = (
        {
            "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_TILE_COUPLING_NOT_MODELLED"],
            "criterion": "pointing <= 8 deg and gain dip <= 3 dB for theta <= 30; 30-45 reported as is",
            "sizes": {},
        },
        None,
    )
    t0 = time.time()
    for R in a.sizes:
        tiles, F = tile_fields(table, green, R, th, ph)
        _, F2 = tile_fields(table, green, R, th2, ph2)
        with mp.get_context("spawn").Pool(
            min(40, os.cpu_count() - 2), initializer=_init, initargs=(F, th, ph, F2, th2, ph2, a.prefer)
        ) as pool:
            rows = pool.map(one, dirs)
        base = sa.scan_summary(sa.scan_table(F, th, ph, dirs))
        res["sizes"][f"R{R:g}"] = {
            "tuned": summary(rows),
            "mrc3_untuned": base,
            "method_counts": {m: sum(r["method"] == m for r in rows) for m in {r["method"] for r in rows}},
            "rows": rows,
        }
        print(
            f"[{time.time() - t0:4.0f}s] R{R:g} tuned: "
            + json.dumps({k: summary(rows)[k] for k in ("inner_le_30", "pass_relaxed_inner")}, default=float),
            flush=True,
        )
        sm = summary(rows)
        key = (
            sm["pass_relaxed_inner"],
            -round(sm["inner_le_30"]["max_pointing_error_deg"], 1),
            sm["inner_le_30"]["min_gain_dbi"],
        )  # worst case to 0.1 deg, ties by gain
        if best is None or key > best[0]:
            best = (key, R, tiles, F)
    res["chosen_R_mm"] = best[1]
    (ROOT / "results/phase29/compact7_tune.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(
        ROOT / "results/phase29/compact7_fields.npz",
        R=best[1],
        theta=th,
        phi=ph,
        **{f"e{j}_{k}": f[k] for j, f in enumerate(best[3], 1) for k in (0, 1)},
    )


if __name__ == "__main__":
    main()
