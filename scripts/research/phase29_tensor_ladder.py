"""Phase 29: re-synthesis on the angle-dependent tensor model (MODEL_ONLY). Unattended driver, parallel (BLAS 1 per worker, memory-guarded):

    python scripts/research/phase29_tensor_ladder.py all

V   validation of the new forward model (item 2): centre-fed single feed vs Phase 26 (info), translation covariance (hexagon + feed shifted by whole lattice pitches) and 180 deg rotation covariance
    (beam azimuth + 180): gain within 0.5 dB, SLL within 2 dB of the unshifted run (pass), plus an off-centre feed in the fixed hexagon (info).
R1  4 feeds at R = 300 mm (Frey scale), pencil beams tilt 7.3 deg at azimuth 0/90/180/270, GA over spacing 1.35-1.5 lambda0, rotation, 4 offsets; fitness mean gain - 0.3 mean max(0, SLL + 20); target SLL <= -15 dB;
    single-feed control at the same size.
R2  7 feeds at R = 150 mm: GA pass 1 at outer tilt 13 deg (defocus allowed), the achieved outer HPBW then sets the tilt of pass 2 (adjacent beams one HPBW apart: crossover ~ -3 dB); +-20 deg scan test with tuned
    weights (<= 3 ports), also reported inside the covered cone (tilt + HPBW / 2). Fitness: mean gain - 0.5 max(0, pointing - 2) - 1.0 max(0, dip - 2.5) - 0.15 mean max(0, |HPBW - target| - 2) - 0.3 mean max(0, SLL + 15).
Also writes the new-model predictions of the two-feed HFSS package (`phase29_twofeed_export.py --new-model`)."""

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
ATTR = "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01UaAn25MEtRcvrHeagn5k3z"
_S = {}


def log(msg):
    with open(RES / "tensor_progress.txt", "a", encoding="utf-8") as f:
        f.write(f"{dt.datetime.now():%m-%d %H:%M} {msg}\n")


def init():
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_tensor as st
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    _S["table"] = SupercellTable.from_files(t.TABLE_GLOB)
    _S["green"] = c.green_function()
    _S["lut"] = st.AngleLUT(RES / "angle_lut_v1.npz", _S["table"])


# ---------------------------------------------------------------- R1 / R2 designs
def ports4(g, tilt=7.3):
    ring = g[0] / (2 * np.sin(np.pi / 4))
    return [
        (ring * np.cos(np.radians(90 * k + 180 + g[1])), ring * np.sin(np.radians(90 * k + 180 + g[1])), tilt, 90.0 * k)
        for k in range(4)
    ], g[2:6]


def ports7(g, tilt):
    from src.hologram import shared_aperture_frey as fr

    return fr.feed_ports(g[0], g[1], g[2:8], tilt_deg=tilt), g[8:15], g[15:22]


def design_for(g, kind, tilt, r_mm, n_psi=180):
    from src.hologram import shared_aperture_tensor as st

    if kind == "n4":
        p, off = ports4(g, tilt)
        return st.design_tensor(_S["table"], _S["lut"], r_mm, p, off, n_psi=n_psi)
    if kind == "n1":
        return st.design_tensor(_S["table"], _S["lut"], r_mm, [(0.0, 0.0, 0.0, 0.0)], [0.0], n_psi=n_psi)
    p, off, dfo = ports7(g, tilt)
    return st.design_tensor(_S["table"], _S["lut"], r_mm, p, off, defocus_rad=dfo, n_psi=n_psi)


def evaluate(g, kind, tilt, r_mm, hpbw_t=None, cone=20.0, grid=(2.0, 3.0), n_psi=180, full=False):
    import phase29_pointB3 as b3
    from phase29_pointB_sweep import port_metrics

    from src.hologram import shared_aperture_multiport as sa

    d = design_for(g, kind, tilt, r_mm, n_psi)
    th, ph = np.arange(0.0, 90.0 + 1e-9, grid[0]), np.arange(0.0, 360.0, grid[1])
    F = sa.port_fields(d, _S["green"], th, ph)
    pm = port_metrics(F, th, ph, d["beams"])
    gain = float(np.mean([p["peak_gain_dbi"] for p in pm]))
    sll = np.array([p["sll_db"] for p in pm])
    hp = np.array([p["hpbw_deg"] for p in pm])
    out = {"mean_gain": gain, "mean_sll": float(sll.mean()), "hpbw": hp.tolist()}
    if kind == "n4":
        out["score"] = float(gain - 0.3 * np.mean(np.maximum(0.0, sll + 20.0)))
    elif kind == "n7":
        i = b3.inner_stats(sa.scan_table(F, th, ph, sa.scan_directions(cone, cone / 5, 30.0)))
        out.update(pointing_max=i["max_pointing_error_deg"], dip=i["max_gain_dip_db"])
        out["score"] = float(
            gain
            - 0.5 * max(0.0, i["max_pointing_error_deg"] - 2.0)
            - 1.0 * max(0.0, i["max_gain_dip_db"] - 2.5)
            - 0.15 * np.mean(np.maximum(0.0, np.abs(hp - hpbw_t) - 2.0))
            - 0.3 * np.mean(np.maximum(0.0, sll + 15.0))
        )
    else:
        out["score"] = gain
    if full:
        out.update({"design": d, "fields": F, "ports_metrics": pm})
    return out


def _fit(a):
    g, kind, tilt, r_mm, hpbw_t, cone, grid = a
    try:
        return evaluate(np.asarray(g), kind, tilt, r_mm, hpbw_t, cone, grid)["score"]
    except Exception:  # an infeasible individual is simply bad
        return -99.0


def run_ga(kind, lo, hi, init_pop, tilt, r_mm, hpbw_t, cone, grid, pop_n, gens, minutes, tag, seed=1):
    from phase29_pointB_sweep import n_workers

    rng = np.random.default_rng(seed)
    pop = lo + rng.random((pop_n, lo.size)) * (hi - lo)
    for k, g0 in enumerate(init_pop):
        pop[k] = g0
    nw = n_workers(per_worker_gb=1.5, reserve_gb=8.0)
    t0, hist = time.time(), []
    log(f"GA {tag}: pop {pop_n}, gens {gens}, {nw} workers")
    with mp.get_context("spawn").Pool(nw, initializer=init) as pool:
        mk = lambda p: (p, kind, tilt, r_mm, hpbw_t, cone, grid)  # noqa: E731
        score = np.array(pool.map(_fit, [mk(p) for p in pop]))
        for gen in range(gens):
            order = np.argsort(-score)
            hist.append({"gen": gen, "best": float(score[order[0]]), "t_s": time.time() - t0})
            log(f"GA {tag} gen {gen}: best {score[order[0]]:.2f} mean {np.mean(score):.2f} [{time.time() - t0:.0f}s]")
            if (time.time() - t0) / 60 > minutes:
                break
            new = [pop[order[0]], pop[order[1]]]
            while len(new) < pop_n:
                i, j = (min(rng.integers(0, pop_n, 3), key=lambda k: -score[k]) for _ in range(2))
                ch = (
                    pop[i]
                    + (rng.random(lo.size) * 1.6 - 0.3) * (pop[j] - pop[i])
                    + (rng.random(lo.size) < 0.12) * rng.normal(size=lo.size) * 0.1 * (hi - lo)
                )
                new.append(np.clip(ch, lo, hi))
            pop = np.array(new)
            score = np.concatenate([score[order[:2]], np.array(pool.map(_fit, [mk(p) for p in pop[2:]]))])
    best = pop[int(np.argmax(score))]
    (RES / f"tensor_ga_{tag}.json").write_text(
        json.dumps({"best_genes": best.tolist(), "best_score": float(score.max()), "history": hist})
    )
    return best


# ---------------------------------------------------------------- validation
def validation():
    import phase29_B3_debug as dbg

    from src.hologram import shared_aperture_tensor as st

    table, lut = _S["table"], _S["lut"]
    dbg._S.update(table=table, green=_S["green"])
    out = {}

    def run(ports, centre=(0.0, 0.0)):
        d = st.design_tensor(table, lut, 150.0, ports, [0.0] * len(ports), centre_mm=centre, n_psi=360)
        m = dbg.field_metrics(d["J"][0], d["x_m"], d["y_m"], (ports[0][2], ports[0][3]), d["eta_conv"][0])
        return {
            k: round(m[k], 2)
            for k in (
                "D_ap_dbi",
                "gain_launch_bare_dbi",
                "sll_15deg_cone_db",
                "sll_1p5hpbw_db",
                "hpbw_deg",
                "ar_peak_db",
                "eta_conv",
            )
        }

    out["V0_centre_fed_new_model"] = run([(0.0, 0.0, 30.0, 225.0)])
    out["V1_translated_by_(12.8,6.4)mm"] = run([(12.8, 6.4, 30.0, 225.0)], (12.8, 6.4))
    out["V2_rotated_180deg"] = run([(0.0, 0.0, 30.0, 45.0)])
    out["V3_off_centre_feed_(25.6,0)_in_fixed_hexagon_info"] = run([(25.6, 0.0, 30.0, 225.0)])
    a = out["V0_centre_fed_new_model"]
    ok = True
    for k in ("V1_translated_by_(12.8,6.4)mm", "V2_rotated_180deg"):
        b = out[k]
        b["delta_gain_db"] = round(b["gain_launch_bare_dbi"] - a["gain_launch_bare_dbi"], 2)
        b["delta_sll_db"] = round(b["sll_15deg_cone_db"] - a["sll_15deg_cone_db"], 2)
        ok &= abs(b["delta_gain_db"]) <= 0.5 and abs(b["delta_sll_db"]) <= 2.0
    out["old_model_reference_phase26"] = {"D_ap_dbi": 28.48, "gain_launch_bare_dbi": 23.51, "sll_15deg_cone_db": -22.77}
    out["covariance_pass"] = bool(ok)
    return out


def main():
    import phase29_B3_debug as dbg
    import phase29_pointB3 as b3
    from phase29_pointB_sweep import n_workers  # noqa: F401

    S = {"started": f"{dt.datetime.now():%Y-%m-%d %H:%M}"}
    log("start")
    init()
    b3._init()
    dbg._init()
    lut_meta = json.loads((RES / "angle_lut_v1.json").read_text())
    S["lut"] = {
        "generated": lut_meta["generated"],
        "mesh": lut_meta["mesh"],
        "symmetry": lut_meta["symmetry_checks"],
        "nan_beta": lut_meta["nan_beta_count"],
    }
    S["V"] = validation()
    (RES / "tensor_summary.json").write_text(json.dumps(S, indent=1, default=float))
    log(f"validation covariance_pass={S['V']['covariance_pass']}")
    # R1: 4 feeds at 300 mm
    lo4, hi4 = np.array([1.35 * LAM, -30.0] + [0.0] * 4), np.array([1.50 * LAM, 30.0] + [360.0] * 4)
    g4 = run_ga(
        "n4",
        lo4,
        hi4,
        [np.array([1.35 * LAM, 0.0, 0.0, 0.0, 0.0, 0.0])],
        7.3,
        300.0,
        None,
        20.0,
        (1.0, 1.5),
        46,
        14,
        40,
        "n4_300",
    )
    res4 = {}
    for kind, g in (("n4", g4), ("n1", None)):
        o = evaluate(g if g is not None else np.zeros(1), kind, 7.3, 300.0, grid=(1.0, 1.0), n_psi=240, full=True)
        res4[kind] = {
            "gain": [round(p["peak_gain_dbi"], 1) for p in o["ports_metrics"]],
            "sll_1p5hpbw": [round(p["sll_db"], 1) for p in o["ports_metrics"]],
            "hpbw": [round(p["hpbw_deg"], 1) for p in o["ports_metrics"]],
            "eta_conv": [round(e, 2) for e in o["design"]["eta_conv"]],
            "mean_gain": round(o["mean_gain"], 1),
            "mean_sll": round(o["mean_sll"], 1),
        }
    res4["target_sll_le_-15"] = bool(res4["n4"]["mean_sll"] <= -15.0)
    S["R1_four_feed_300mm"] = res4
    (RES / "tensor_summary.json").write_text(json.dumps(S, indent=1, default=float))
    log(f"R1 mean SLL {res4['n4']['mean_sll']} (control {res4['n1']['mean_sll']})")
    # R2: 7 feeds at 150 mm, tilt from the achieved HPBW
    lo7 = np.array([1.35 * LAM, -30.0] + [-4.0] * 6 + [0.0] * 7 + [0.0] * 7)
    hi7 = np.array([1.50 * LAM, 30.0] + [4.0] * 6 + [360.0] * 7 + [7.0] * 7)
    seed7 = np.concatenate([[40.0, 0.0], np.zeros(6), np.zeros(7), np.full(7, 3.0)])
    tilt = 13.0
    gbest, tilts = None, []
    for ps in (1, 2):
        gbest = run_ga(
            "n7",
            lo7,
            hi7,
            [seed7] if gbest is None else [gbest, seed7],
            tilt,
            150.0,
            tilt,
            20.0,
            (2.0, 3.0),
            92,
            30,
            35 if ps == 1 else 25,
            f"n7_pass{ps}",
            seed=ps,
        )
        o = evaluate(gbest, "n7", tilt, 150.0, tilt, 20.0, grid=(1.0, 1.0), n_psi=240, full=True)
        hp = float(np.mean([p["hpbw_deg"] for p in o["ports_metrics"][1:]]))
        tilts.append({"pass": ps, "tilt_deg": tilt, "achieved_outer_hpbw_deg": hp})
        log(f"R2 pass {ps}: tilt {tilt:.1f} achieved outer HPBW {hp:.1f}")
        if ps == 1:
            tilt = round(hp, 1)
    S["R2_tilt_iteration"] = tilts
    # final R2 evaluation (tuned weights, <= 3 ports)
    import phase29_compact7_tune as ct

    from src.hologram import shared_aperture_multiport as sa

    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    th2, ph2 = np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 2.0)
    d, F, pm = o["design"], o["fields"], o["ports_metrics"]
    dirs = sa.scan_directions(20.0, 2.0, 15.0)
    F2 = sa.port_fields(d, _S["green"], th2, ph2)
    with mp.get_context("spawn").Pool(
        min(40, n_workers(per_worker_gb=1.5)), initializer=ct._init, initargs=(F, th, ph, F2, th2, ph2, 1.5)
    ) as pool:
        tuned = pool.map(ct.one, dirs)
    gmax = max(r["peak_gain_dbi"] for r in tuned)
    err = np.array([r["pointing_error_deg"] for r in tuned])
    dip = gmax - np.array([r["peak_gain_dbi"] for r in tuned])
    covered = tilt + hp / 2
    inside = np.array([r["theta_deg"] <= covered + 1e-9 for r in tuned])
    S["R2_seven_feed_150mm"] = {
        "tilt_deg": tilt,
        "outer_hpbw_deg": hp,
        "gain": [round(p["peak_gain_dbi"], 1) for p in pm],
        "hpbw": [round(p["hpbw_deg"], 1) for p in pm],
        "sll_1p5hpbw": [round(p["sll_db"], 1) for p in pm],
        "crossover_mean_db": {k: round(v["mean_db"], 1) for k, v in b3.crossovers(F, th, ph, pm).items()},
        "sum_diff": {
            k: {
                kk: (round(vv, 1) if isinstance(vv, float) else vv)
                for kk, vv in v.items()
                if kk in ("peak_gain_dbi", "sll_db", "null_depth_at_sum_peak_db")
            }
            for k, v in b3.sum_diff(F, th, ph).items()
        },
        "scan_pm20": {
            "n": len(tuned),
            "max_pointing_error_deg": round(float(err.max()), 1),
            "median_pointing_error_deg": round(float(np.median(err)), 1),
            "fraction_within_2deg": round(float(np.mean(err <= 2)), 2),
            "max_gain_dip_db": round(float(dip.max()), 1),
            "min_gain_dbi": round(float(min(r["peak_gain_dbi"] for r in tuned)), 1),
        },
        "covered_cone_deg": round(covered, 1),
        "scan_inside_covered_cone": {
            "n": int(inside.sum()),
            "max_pointing_error_deg": round(float(err[inside].max()), 1),
            "max_gain_dip_db": round(float(dip[inside].max()), 1),
            "fraction_within_2deg": round(float(np.mean(err[inside] <= 2)), 2),
        },
        "pass_pm20": bool(err.max() <= 2.0 and dip.max() <= 3.0),
        "pass_covered_cone": bool(err[inside].max() <= 2.0 and dip[inside].max() <= 3.0),
    }
    (RES / "tensor_summary.json").write_text(json.dumps(S, indent=1, default=float))
    # two-feed package: new-model predictions
    r = subprocess.run(
        [sys.executable, "scripts/research/phase29_twofeed_export.py", "--new-model"],
        cwd=ROOT,
        env=dict(os.environ, PYTHONPATH=str(ROOT)),
        capture_output=True,
        text=True,
    )
    log(f"twofeed new-model predictions rc={r.returncode} {r.stderr.strip()[-200:]}")
    p = ROOT / "docs/STATUS.md"
    s = p.read_text(encoding="latin-1")
    i = s.index("## Open next steps")
    row = (
        "\n## Phase 29 angle-dependent tensor model + re-synthesis (2026-10-28) - MODEL_ONLY\n\n"
        "LUT `results/phase29/angle_lut_v1.npz` (slot x orientation x incidence, Bloch sheet tensor from periodic MoM), forward model `src/hologram/shared_aperture_tensor.py`, validation, 4 feeds at 300 mm and 7 feeds at 150 mm re-synthesised; "
        "numbers in `results/phase29/tensor_summary.json`. Two-feed HFSS check package: `hfss_exports/phase29_twofeed/` (old and new model predictions).\n\n"
    )
    p.write_text(s[:i] + row + s[i:], encoding="latin-1", newline="")
    for a in (
        ["git", "add", "-A"],
        [
            "git",
            "commit",
            "-q",
            "-m",
            f"Phase 29: angle-dependent tensor cell LUT + forward model + Frey-style re-synthesis (4 feeds 300 mm, 7 feeds 150 mm), new-model predictions for the two-feed package (MODEL_ONLY)\n\n{ATTR}",
        ],
        ["git", "push", "hub", "phase29-seven-hex"],
    ):
        rr = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
        log(f"{a[1]} rc={rr.returncode}")
    (RES / "tensor_DONE").write_text("done\n")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "all":
        main()
