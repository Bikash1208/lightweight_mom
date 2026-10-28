"""Phase 29 B3 diagnostic at Frey's scale (MODEL_ONLY): shared hexagon R = 300 mm (circumradius, project convention; = 10.5 lambda0), same angle-dependent synthesis + GA.
Question: is the high sidelobe level of the 150 mm B3 designs due to the aperture size or to the repo's cell model / synthesis?

    python scripts/research/phase29_B3_scale300.py all      # unattended: GA 7 feeds, GA 4 feeds, N = 1 baseline, evaluation, summary, docs, commit, push

7 feeds: centre + ring (1.35-1.5 lambda0), cone +-11 deg, outer tilt = HPBW target = 7.3 deg (tilt ~ one beamwidth, crossover ~ -3 dB; cone = tilt + HPBW/2).
4 feeds: ring (adjacent spacing 1.35-1.5 lambda0), pencil beams (no defocus), tilt 7.3 deg at azimuths 0, 90, 180, 270 (the comparison with Frey: 23 dBi, SLL 15-20 dB).
N = 1: the centre feed alone (broadside), same cells model: the control for superposition error.
GA fitness now includes an SLL term (the 150 mm GA never optimised sidelobes): 7 feeds: mean gain - 0.5 max(0, pointing - 2) - 1.0 max(0, dip - 2.5) - 0.15 mean max(0, |HPBW - 7.3| - 1.5)
- 0.3 mean max(0, SLL + 15); 4 feeds: mean gain - 0.3 mean max(0, SLL + 20). DoF ceiling (user): D_avg <= 4 pi N / Omega, Omega = 2 pi (1 - cos theta_max)."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import datetime as dt  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/research"))
RES = ROOT / "results/phase29"
LAM = 299792458.0 / 10.54e9 * 1e3
R_MM, TILT, HPBW_T, CONE = 300.0, 7.3, 7.3, 11.0
ATTR = "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01UaAn25MEtRcvrHeagn5k3z"
LO7 = np.array([1.35 * LAM, -30.0] + [-4.0] * 6 + [0.0] * 7 + [0.0] * 7 + [1.0])
HI7 = np.array([1.50 * LAM, 30.0] + [4.0] * 6 + [360.0] * 7 + [14.0] * 7 + [4.0])
LO4 = np.array([1.35 * LAM, -30.0] + [0.0] * 4 + [1.0])
HI4 = np.array([1.50 * LAM, 30.0] + [360.0] * 4 + [4.0])


def log(msg):
    with open(RES / "B3s300_progress.txt", "a", encoding="utf-8") as f:
        f.write(f"{dt.datetime.now():%m-%d %H:%M} {msg}\n")


def decode(g, n):
    from src.hologram import shared_aperture_frey as fr

    if n == 7:
        return fr.feed_ports(g[0], g[1], g[2:8], tilt_deg=TILT), g[8:15], g[15:22], 10.0 ** g[22]
    if n == 4:
        ring = g[0] / (2 * np.sin(np.pi / 4))
        ports = [
            (
                ring * np.cos(np.radians(90 * k + 180 + g[1])),
                ring * np.sin(np.radians(90 * k + 180 + g[1])),
                TILT,
                90.0 * k,
            )
            for k in range(4)
        ]
        return ports, g[2:6], np.zeros(4), 10.0 ** g[6]
    return [(0.0, 0.0, 0.0, 0.0)], np.zeros(1), np.zeros(1), 10.0 ** g[0]


def design(g, n, n_psi=240):
    import phase29_pointB3 as b3

    from src.hologram import shared_aperture_frey as fr

    ports, off, dfo, q0 = decode(g, n)
    return fr.design_frey(b3._S["table"], R_MM, ports, off, 1000.0, q0, n_iter=2, n_psi=n_psi, defocus_rad=dfo)


def evaluate(g, n, grid=(2.0, 3.0), n_psi=180, full=False):
    import phase29_pointB3 as b3
    from phase29_pointB_sweep import port_metrics

    from src.hologram import shared_aperture_multiport as sa

    d = design(g, n, n_psi)
    th, ph = np.arange(0.0, 90.0 + 1e-9, grid[0]), np.arange(0.0, 360.0, grid[1])
    F = sa.port_fields(d, b3._S["green"], th, ph)
    pm = port_metrics(F, th, ph, d["beams"])
    gain = float(np.mean([p["peak_gain_dbi"] for p in pm]))
    sll = np.array([p["sll_db"] for p in pm])
    hp = np.array([p["hpbw_deg"] for p in pm])
    out = {"mean_gain": gain, "mean_sll": float(sll.mean()), "hpbw": hp.tolist()}
    if n == 7:
        i = b3.inner_stats(sa.scan_table(F, th, ph, sa.scan_directions(CONE, CONE / 5, 30.0)))
        out.update(pointing_max=i["max_pointing_error_deg"], dip=i["max_gain_dip_db"])
        out["score"] = float(
            gain
            - 0.5 * max(0.0, i["max_pointing_error_deg"] - 2.0)
            - 1.0 * max(0.0, i["max_gain_dip_db"] - 2.5)
            - 0.15 * np.mean(np.maximum(0.0, np.abs(hp - HPBW_T) - 1.5))
            - 0.3 * np.mean(np.maximum(0.0, sll + 15.0))
        )
    else:
        out["score"] = float(gain - 0.3 * np.mean(np.maximum(0.0, sll + 20.0)))
    if full:
        out.update({"design": d, "fields": F, "ports_metrics": pm})
    return out


def _fit(a):
    g, n = a
    try:
        return evaluate(np.asarray(g), n, grid=(1.0, 1.5) if n == 4 else (2.0, 3.0))["score"]
    except Exception:  # an infeasible individual is simply bad
        return -99.0


def run_ga(n, lo, hi, init, pop_n, gens, minutes, seed=1):
    import phase29_pointB3 as b3
    from phase29_pointB_sweep import n_workers

    rng = np.random.default_rng(seed)
    pop = lo + rng.random((pop_n, lo.size)) * (hi - lo)
    pop[0] = init
    nw = n_workers(per_worker_gb=1.5, reserve_gb=8.0)
    t0, hist = time.time(), []
    log(f"GA n={n}: pop {pop_n}, gens {gens}, {nw} workers")
    with mp.get_context("spawn").Pool(nw, initializer=b3._init) as pool:
        score = np.array(pool.map(_fit, [(p, n) for p in pop]))
        for gen in range(gens):
            order = np.argsort(-score)
            hist.append({"gen": gen, "best": float(score[order[0]]), "t_s": time.time() - t0})
            log(f"GA n={n} gen {gen}: best {score[order[0]]:.2f} mean {np.mean(score):.2f} [{time.time() - t0:.0f}s]")
            if (time.time() - t0) / 60 > minutes:
                break
            new = [pop[order[0]], pop[order[1]]]
            while len(new) < pop_n:
                i, j = (min(rng.integers(0, pop_n, 3), key=lambda k: -score[k]) for _ in range(2))
                child = pop[i] + (rng.random(lo.size) * 1.6 - 0.3) * (pop[j] - pop[i])
                child = child + (rng.random(lo.size) < 0.12) * rng.normal(size=lo.size) * 0.1 * (hi - lo)
                new.append(np.clip(child, lo, hi))
            pop = np.array(new)
            score = np.concatenate([score[order[:2]], np.array(pool.map(_fit, [(p, n) for p in pop[2:]]))])
    best = pop[int(np.argmax(score))]
    (RES / f"B3s300_ga{n}.json").write_text(
        json.dumps({"n": n, "best_genes": best.tolist(), "best_score": float(score.max()), "history": hist})
    )
    return best


def sidelobe_diag(f, th, ph, pm):
    """Where the power and the worst sidelobe of one port pattern are: fraction of the radiated power inside 1.5 HPBW of the peak, LHCP fraction of the total,
    worst RHCP sidelobe outside that cone (level re peak, angular distance, theta/phi), worst LHCP level re the RHCP peak."""
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa
    from src.pattern.vector_aperture_pattern import copol_component

    et, ep = f
    tg, pg = np.meshgrid(th, ph, indexing="ij")
    w = np.sin(np.radians(np.maximum(tg, 0.25)))
    inten = np.abs(et) ** 2 + np.abs(ep) ** 2
    gr = sa.gain_db(et, ep)
    rh, lh = np.abs(sa.rhcp_field(et, ep)) ** 2, np.abs(copol_component(et, ep, "lhcp")) ** 2
    pk = tuple(pm["peak_theta_phi"])
    ang = np.degrees(np.arccos(np.clip(np.sum(c._unit(tg, pg) * c._unit(*pk), axis=-1), -1, 1)))
    out = ang > 1.5 * pm["hpbw_deg"]
    k = np.unravel_index(np.argmax(np.where(out, gr, -999)), gr.shape)
    return {
        "power_fraction_in_1p5hpbw": float(np.sum(inten * w * ~out) / np.sum(inten * w)),
        "lhcp_power_fraction": float(np.sum(lh * w) / np.sum((lh + rh) * w)),
        "worst_rhcp_sidelobe_db": float(gr[k] - gr.max()),
        "worst_sidelobe_deg_from_peak": float(ang[k]),
        "worst_sidelobe_theta_phi": [float(tg[k]), float(pg[k])],
        "worst_lhcp_re_rhcp_peak_db": float(10 * np.log10(max(lh.max(), 1e-30) / rh.max())),
    }


def evaluate_full(n, g):
    import phase29_compact7_tune as ct
    import phase29_pointB3 as b3
    from phase29_pointB_sweep import n_workers

    from src.hologram import shared_aperture_multiport as sa

    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    o = evaluate(g, n, grid=(1.0, 1.0), n_psi=240, full=True)
    d, F, pm = o["design"], o["fields"], o["ports_metrics"]
    res = {
        "n_feeds": n,
        "n_cells": int(d["x_m"].size),
        "slots": np.bincount(d["slot"], minlength=4).tolist(),
        "eta_conv": d["eta_conv"],
        "genes": list(map(float, g)),
        "ports": pm,
        "diag": [sidelobe_diag(f, th, ph, p) for f, p in zip(F, pm, strict=True)],
        "mean_gain_dbi": o["mean_gain"],
        "mean_sll_db": o["mean_sll"],
    }
    if n == 7:
        res["crossover"], res["sum_diff"] = b3.crossovers(F, th, ph, pm), b3.sum_diff(F, th, ph)
        dirs = sa.scan_directions(CONE, 1.0, 15.0)
        th2, ph2 = np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 2.0)
        F2 = sa.port_fields(d, b3._S["green"], th2, ph2)
        with mp.get_context("spawn").Pool(
            min(40, n_workers(per_worker_gb=1.5)), initializer=ct._init, initargs=(F, th, ph, F2, th2, ph2, 1.5)
        ) as pool:
            tuned = pool.map(ct.one, dirs)
        gmax = max(r["peak_gain_dbi"] for r in tuned)
        err = np.array([r["pointing_error_deg"] for r in tuned])
        dip = gmax - np.array([r["peak_gain_dbi"] for r in tuned])
        edge = [
            r for r in tuned if abs(r["theta_deg"] - dirs[-1][0]) < 1e-9 and abs((r["phi_deg"] - 30.0) % 60.0) < 1e-9
        ]
        res["scan_tuned"] = {
            "n": len(tuned),
            "max_pointing_error_deg": float(err.max()),
            "median_pointing_error_deg": float(np.median(err)),
            "fraction_within_2deg": float(np.mean(err <= 2)),
            "max_gain_dip_db": float(dip.max()),
            "min_gain_dbi": float(min(r["peak_gain_dbi"] for r in tuned)),
            "pass_pointing_2deg": bool(err.max() <= 2),
            "pass_dip_3db": bool(dip.max() <= 3),
            "cone_edge_between_outer_beams": [
                [r["phi_deg"], r["pointing_error_deg"], r["peak_gain_dbi"]] for r in edge
            ],
        }
        om = 2 * np.pi * (1 - np.cos(np.radians(CONE)))
        ceil = 10 * np.log10(4 * np.pi * 7 / om)
        res["ceiling_gap_db"] = {
            "ceiling_dbi_user_formula": float(ceil),
            "gap_mean_port_db": float(ceil - o["mean_gain"]),
            "gap_min_scan_db": float(ceil - res["scan_tuned"]["min_gain_dbi"]),
        }
    return res


def summary_line(r):
    d = r["diag"]
    return {
        "gain": [round(p["peak_gain_dbi"], 1) for p in r["ports"]],
        "hpbw": [round(p["hpbw_deg"], 1) for p in r["ports"]],
        "sll": [round(p["sll_db"], 1) for p in r["ports"]],
        "worst_sidelobe_deg_from_peak": [round(x["worst_sidelobe_deg_from_peak"]) for x in d],
        "main_lobe_power_fraction": [round(x["power_fraction_in_1p5hpbw"], 2) for x in d],
        "lhcp_power_fraction": [round(x["lhcp_power_fraction"], 2) for x in d],
    }


def main():
    import phase29_pointB3 as b3

    S = {"started": f"{dt.datetime.now():%Y-%m-%d %H:%M}", "R_mm": R_MM, "cone_deg": CONE, "tilt_deg": TILT}
    b3._init()
    log("start")
    t0 = time.time()
    g7 = run_ga(
        7, LO7, HI7, np.concatenate([[40.0, 0.0], np.zeros(6), np.zeros(7), np.full(7, 5.0), [3.0]]), 92, 40, 70
    )
    g4 = run_ga(4, LO4, HI4, np.array([40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0]), 46, 30, 25)
    log(f"GAs done [{(time.time() - t0) / 60:.0f} min]")
    g1 = np.array([g4[6]])
    for n, g in ((7, g7), (4, g4), (1, g1)):
        try:
            r = evaluate_full(n, g)
            (RES / f"B3s300_n{n}.json").write_text(json.dumps(r, indent=1, default=float))
            S[f"n{n}"] = summary_line(r) | {
                k: (r[k] if not isinstance(r[k], float) else round(r[k], 2))
                for k in ("n_cells", "mean_gain_dbi", "mean_sll_db")
            }
            if n == 7:
                S["n7"].update(
                    {
                        "crossover": {k: round(v["mean_db"], 1) for k, v in r["crossover"].items()},
                        "sum": {
                            k: (round(v, 1) if isinstance(v, float) else [round(x) for x in v])
                            for k, v in r["sum_diff"]["sum"].items()
                        },
                        "delta_null_db": {
                            k: round(r["sum_diff"][k]["null_depth_at_sum_peak_db"], 1) for k in ("delta_x", "delta_y")
                        },
                        "scan": {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r["scan_tuned"].items()},
                        "ceiling": {k: round(v, 1) for k, v in r["ceiling_gap_db"].items()},
                    }
                )
            log(f"final n={n} done")
        except Exception as exc:  # recorded, never fatal
            S[f"n{n}_error"] = repr(exc)[:300]
            log(f"final n={n} FAILED {exc!r}"[:300])
    (RES / "B3s300_summary.json").write_text(json.dumps(S, indent=1))
    # docs, commit, push (no fast suite: only new files; the new tests are run by the caller)
    p = ROOT / "docs/STATUS.md"
    s = p.read_text(encoding="latin-1")
    i = s.index("## Open next steps")
    row = (
        "\n## Phase 29 B3 diagnostic at Frey's scale (R 300 mm, 2026-10-28) - MODEL_ONLY, no HFSS package\n\n"
        "Is the high SLL of the 150 mm B3 designs due to aperture size or to the cell model / synthesis? Numbers: `results/phase29/B3s300_summary.json` (7 feeds, 4 feeds, single-feed baseline).\n\n"
    )
    p.write_text(s[:i] + row + s[i:], encoding="latin-1", newline="")
    for a in (
        ["git", "add", "-A"],
        [
            "git",
            "commit",
            "-q",
            "-m",
            f"Phase 29 B3 diagnostic at R 300 mm (Frey scale): 7-feed, 4-feed and single-feed runs (MODEL_ONLY)\n\n{ATTR}",
        ],
        ["git", "push", "hub", "phase29-seven-hex"],
    ):
        r = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
        log(f"{a[1]} rc={r.returncode} {r.stderr.strip()[-100:]}")
    (RES / "B3s300_DONE").write_text("done\n")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "all":
        main()
