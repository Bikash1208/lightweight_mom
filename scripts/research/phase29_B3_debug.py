"""Phase 29 B3 debugging ladder (MODEL_ONLY): is the B3 code path correct? Stops at the first failing rung.

    python scripts/research/phase29_B3_debug.py        # writes results/phase29/B3debug_summary.json (+ progress lines in B3debug_progress.txt)

Rung 1 (150 mm, ONE feed through the B3 code path, beam theta 30 / phi 225): must reproduce Phase 26 hex150_tensor_cp_single (D_ap 28.48 dBi, gain 23.5 dBi with the bare launch fraction 0.592,
SLL -22.8 dB, `hfss_exports/hex150_tensor_cp_single/predictions.json`). Three experiments isolate the cause of a difference:
  (A) the committed Phase 26 cells run through the B3 CHAIN (analysis march, incident vector, tensor current, far field) -> chain check;
  (B) the B3 DESIGN rule (`march_design`, global conversion density q0) for several q0 -> what B3 actually does;
  (C) the Phase 26 design rule (`march_rays`, eta_design 0.9, taper 0.5) inside the B3 path (`design_frey(march_fn=...)`) -> design-rule check.
Rung 2 (only if rung 1 passes): Frey replication scaled to 10.54 GHz: 4 feeds, hexagon of the same area as their 12.3 lambda0 square (R = 7.62 lambda0), adjacent feed spacing 1.35 lambda0 (the user's
value; the paper says 1.5), pairwise-tilted pencil beams (tilt 2.5 deg: feeds 1/2 at azimuth 0/180, 3/4 at 90/270 -> FOV ~10 deg with ~5 deg beams; Table 1 of the paper is not readable from the text),
Phase 26 design rule per feed, Frey superposition, small GA over the 4 hologram offsets (SLL penalised). Target: ~23 dBi per feed, SLL 15-20 dB below the peak.
Also reported for both rungs: slot snapping (discrete levels), continuous-kappa variant, GA fitness terms, pins in the forward model."""

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
LAUNCH_BARE = 0.592
ATTR = "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01UaAn25MEtRcvrHeagn5k3z"
_S = {}


def log(msg):
    with open(RES / "B3debug_progress.txt", "a", encoding="utf-8") as f:
        f.write(f"{dt.datetime.now():%m-%d %H:%M} {msg}\n")


def _init():
    from src.hologram import seven_hex_cluster as c
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    _S["table"], _S["green"] = SupercellTable.from_files(t.TABLE_GLOB), c.green_function()


def field_metrics(J, x, y, beam, eta_conv, theta=None, phi=None):
    """Phase 26 `metrics` (D_ap, SLL outside a 15 deg cone, AR, XPD) of the far field of the currents `J` (n, 2) at the cells, plus the SLL outside 1.5 HPBW and the model gain with the bare
    launch fraction (D_ap + 10 log(0.592 eta_conv): the figure the Phase 26 prediction quotes)."""
    from phase29_pointB_sweep import _hpbw

    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t

    theta = np.arange(0.0, 90.0 + 1e-9, 1.0) if theta is None else theta
    phi = np.arange(0.0, 360.0, 2.0) if phi is None else phi
    et, ep = t.pattern({"jx": J[:, 0], "jy": J[:, 1]}, x, y, _S["green"], theta, phi)
    m = t.metrics(et, ep, theta, phi, beam=tuple(beam))
    g = sa.gain_db(et, ep)
    hp = _hpbw(g, theta, phi)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    k = np.unravel_index(np.argmax(g), g.shape)
    ang = np.degrees(np.arccos(np.clip(np.sum(c._unit(tg, pg) * c._unit(tg[k], pg[k]), axis=-1), -1, 1)))
    d_ap = m["directivity_ap_rhcp_dbi"]
    return {
        "D_ap_dbi": float(d_ap),
        "gain_launch_bare_dbi": float(d_ap + 10 * np.log10(LAUNCH_BARE * eta_conv)),
        "eta_conv": float(eta_conv),
        "sll_15deg_cone_db": float(m["rhcp_sll_db"]),
        "sll_1p5hpbw_db": float(g[ang > 1.5 * hp].max() - g.max()),
        "hpbw_deg": float(hp),
        "ar_peak_db": float(m["ar_at_peak_db"]),
        "peak_theta_phi": [float(m["peak_rhcp_theta_deg"]), float(m["peak_rhcp_phi_deg"])],
    }


def slot_levels(design, table):
    s = np.bincount(design["slot"], minlength=len(table.slots))
    return {
        "slot_LW_mm": [list(x) for x in table.slots],
        "cells_per_slot": s.tolist(),
        "levels_used": int(np.sum(s > 0)),
    }


def design_currents(d, kappa_continuous=False):
    """Per-port currents of a `design_frey` result; `kappa_continuous`: the design strength `kappa_des` (not snapped to the 4 table slots) with the same axes and incident waves."""
    from src.hologram import shared_aperture_multiport as sa

    if not kappa_continuous:
        return d["J"]
    out = []
    for j in range(len(d["J"])):
        s = sa.incident(d["march"][j], d["feeds_m"][j], d["x_m"], d["y_m"])[0]
        out.append(sa.tensor_current(s, d["kappa_des"][j], d["axis_rad"]))
    return out


def rung1(table):
    import json as _j

    from src.hologram import shared_aperture_frey as fr
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import march_rays

    ref_p = _j.loads((ROOT / "hfss_exports/hex150_tensor_cp_single/predictions.json").read_text())
    ap = ref_p["aperture_pattern"]
    ref = {
        "D_ap_dbi": ap["directivity_ap_rhcp_dbi"],
        "gain_launch_bare_dbi": ref_p["realized_gain_rhcp_dbi"]["launch_bare_0.592"],
        "sll_15deg_cone_db": ap["rhcp_sll_db"],
        "ar_peak_db": ap["ar_at_peak_db"],
    }
    k0, th0, ph0 = t.K0, np.radians(30.0), np.radians(225.0)
    out = {"reference_phase26": ref}
    # (A) Phase 26 cells through the B3 chain
    cells = t.single_hexagon_cells(table)
    x, y = np.array([c["x_mm"] for c in cells]) * 1e-3, np.array([c["y_mm"] for c in cells]) * 1e-3
    lw = {tuple(v): i for i, v in enumerate(table.slots)}
    slot = np.array([lw[(round(c["slot_L_mm"], 1), round(c["slot_W_mm"], 1))] for c in cells])
    axis = np.radians([c["slot_angle_deg"] for c in cells])
    poly, feed = sa.hexagon_m(150.0), (0.0, 0.0)
    ana = sa.march_analysis(table, poly, feed, k0, th0, ph0, np.c_[x, y], slot, 720)
    s, alpha, _, _, _ = sa.incident(ana, feed, x, y)
    # the slot angle of Phase 26 is the long axis from +x; the B3 axis `a` is the same quantity (chi = 2a - psi)
    J = sa.tensor_current(s, np.sqrt(2.0 * alpha), axis)
    out["A_phase26_cells_through_B3_chain"] = field_metrics(J, x, y, (30.0, 225.0), ana["eta_conv"])
    log("rung1 A done")
    # (B) B3 design rule, several q0
    ports = [(0.0, 0.0, 30.0, 225.0)]
    out["B_B3_design_rule"] = {}
    best = None
    for q0 in (30.0, 100.0, 300.0, 1000.0, 1e4, 1e5):
        d = fr.design_frey(table, 150.0, ports, [0.0], 1000.0, q0, n_iter=2, n_psi=360)
        m = field_metrics(d["J"][0], d["x_m"], d["y_m"], (30.0, 225.0), d["eta_conv"][0])
        m.update(slot_levels(d, table))
        out["B_B3_design_rule"][f"q0={q0:g}"] = m
        if best is None or m["sll_15deg_cone_db"] < best[1]["sll_15deg_cone_db"]:
            best = (q0, m, d)
    out["B_best_q0"] = best[0]
    mc = field_metrics(
        design_currents(best[2], True)[0], best[2]["x_m"], best[2]["y_m"], (30.0, 225.0), best[2]["eta_conv"][0]
    )
    out["B_best_q0_continuous_kappa"] = mc
    log("rung1 B done")
    # (C) Phase 26 design rule inside the B3 path
    mf = lambda j, f, b: march_rays(table, sa.hexagon_m(150.0), f, k0, b[0], b[1], 0.9, taper_c=0.5, n_psi=720)  # noqa: E731
    d = fr.design_frey(table, 150.0, ports, [0.0], 1000.0, 1.0, n_iter=2, n_psi=720, march_fn=mf)
    out["C_phase26_rule_in_B3_path"] = field_metrics(
        d["J"][0], d["x_m"], d["y_m"], (30.0, 225.0), d["eta_conv"][0]
    ) | slot_levels(d, table)
    out["C_continuous_kappa"] = field_metrics(
        design_currents(d, True)[0], d["x_m"], d["y_m"], (30.0, 225.0), d["eta_conv"][0]
    )
    log("rung1 C done")
    a, b_, c_ = (
        out["A_phase26_cells_through_B3_chain"],
        out["B_B3_design_rule"][f"q0={best[0]:g}"],
        out["C_phase26_rule_in_B3_path"],
    )
    out["verdict"] = {
        "chain_reproduces_phase26 (A: |dD_ap| <= 1 dB and SLL within 3 dB)": bool(
            abs(a["D_ap_dbi"] - ref["D_ap_dbi"]) <= 1.0
            and abs(a["sll_15deg_cone_db"] - ref["sll_15deg_cone_db"]) <= 3.0
        ),
        "B3_design_rule_reproduces_phase26 (B best q0)": bool(
            abs(b_["D_ap_dbi"] - ref["D_ap_dbi"]) <= 1.5 and b_["sll_15deg_cone_db"] <= ref["sll_15deg_cone_db"] + 4.0
        ),
        "phase26_rule_in_B3_path_reproduces (C)": bool(
            abs(c_["D_ap_dbi"] - ref["D_ap_dbi"]) <= 1.5 and c_["sll_15deg_cone_db"] <= ref["sll_15deg_cone_db"] + 4.0
        ),
    }
    out["rung1_pass"] = bool(out["verdict"]["B3_design_rule_reproduces_phase26 (B best q0)"])
    return out


def _fit2(a):
    g, spec = a
    try:
        return rung2_eval(np.asarray(g), spec)["score"]
    except Exception:  # an infeasible individual is simply bad
        return -99.0


def rung2_ports(ring_mm, rot, tilt):
    return [
        (ring_mm * np.cos(np.radians(az + 180 + rot)), ring_mm * np.sin(np.radians(az + 180 + rot)), tilt, az)
        for az in (0.0, 180.0, 90.0, 270.0)
    ]


def rung2_eval(g, spec, grid=(1.0, 1.5), full=False, n_psi=240):
    from src.hologram import shared_aperture_frey as fr
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram.tensor_leaky_cp import march_rays

    k0 = 2 * np.pi * 10.54e9 / 299792458.0
    R = spec["R"]
    ports = rung2_ports(spec["ring"], g[0], spec["tilt"])
    mf = lambda j, f, b: march_rays(_S["table"], sa.hexagon_m(R), f, k0, b[0], b[1], 0.9, taper_c=0.5, n_psi=n_psi)  # noqa: E731
    d = fr.design_frey(_S["table"], R, ports, g[1:5], 1000.0, 1.0, n_iter=2, n_psi=n_psi, march_fn=mf)
    th, ph = np.arange(0.0, 90.0 + 1e-9, grid[0]), np.arange(0.0, 360.0, grid[1])
    rows = [
        field_metrics(J, d["x_m"], d["y_m"], (ports[j][2], ports[j][3]), d["eta_conv"][j], th, ph)
        for j, J in enumerate(d["J"])
    ]
    gain = np.array([r["gain_launch_bare_dbi"] for r in rows])
    sll = np.array([r["sll_1p5hpbw_db"] for r in rows])
    out = {
        "score": float(gain.mean() - 0.3 * np.mean(np.maximum(0.0, sll + 20.0))),
        "gain": gain.tolist(),
        "sll": sll.tolist(),
    }
    if full:
        out.update({"design": d, "rows": rows})
    return out


def rung2(table):
    import phase29_pointB3 as b3  # noqa: F401
    from phase29_pointB_sweep import n_workers

    from src.hologram import shared_aperture_frey as fr
    from src.hologram import tensor_cp_design as t

    R = 7.62 * LAM
    spec = {"R": R, "ring": 1.35 * LAM / (2 * np.sin(np.pi / 4)), "tilt": 2.5}
    out = {"aperture_R_mm": R, "ring_mm": spec["ring"], "spacing_lambda0": 1.35, "tilt_deg": 2.5}
    rng = np.random.default_rng(1)
    lo, hi = np.array([-30.0, 0, 0, 0, 0]), np.array([30.0, 360, 360, 360, 360])
    pop = lo + rng.random((24, 5)) * (hi - lo)
    pop[0] = 0.0
    nw = n_workers(per_worker_gb=2.0, reserve_gb=8.0, cap=24)
    t0 = time.time()
    with mp.get_context("spawn").Pool(nw, initializer=_init) as pool:
        score = np.array(pool.map(_fit2, [(p, spec) for p in pop]))
        for gen in range(10):
            order = np.argsort(-score)
            log(f"rung2 GA gen {gen}: best {score[order[0]]:.2f} [{time.time() - t0:.0f}s]")
            new = [pop[order[0]], pop[order[1]]]
            while len(new) < len(pop):
                i, j = (min(rng.integers(0, len(pop), 3), key=lambda k: -score[k]) for _ in range(2))
                ch = (
                    pop[i]
                    + (rng.random(5) * 1.6 - 0.3) * (pop[j] - pop[i])
                    + (rng.random(5) < 0.15) * rng.normal(size=5) * 0.1 * (hi - lo)
                )
                new.append(np.clip(ch, lo, hi))
            pop = np.array(new)
            score = np.concatenate([score[order[:2]], np.array(pool.map(_fit2, [(p, spec) for p in pop[2:]]))])
    gbest = pop[int(np.argmax(score))]
    _init()
    r0 = rung2_eval(np.zeros(5), spec, full=True)
    r1 = rung2_eval(gbest, spec, full=True)
    out["offsets0"] = {
        "gain_dbi": [round(x, 1) for x in r0["gain"]],
        "sll_1p5hpbw_db": [round(x, 1) for x in r0["sll"]],
        "sll_15deg_cone_db": [round(r["sll_15deg_cone_db"], 1) for r in r0["rows"]],
    }
    out["ga_offsets"] = {
        "genes": gbest.tolist(),
        "gain_dbi": [round(x, 1) for x in r1["gain"]],
        "sll_1p5hpbw_db": [round(x, 1) for x in r1["sll"]],
        "sll_15deg_cone_db": [round(r["sll_15deg_cone_db"], 1) for r in r1["rows"]],
        "hpbw_deg": [round(r["hpbw_deg"], 1) for r in r1["rows"]],
        "ar_peak_db": [round(r["ar_peak_db"], 1) for r in r1["rows"]],
        "peak_theta_phi": [r["peak_theta_phi"] for r in r1["rows"]],
    }
    d = r1["design"]
    out["slot_levels"] = slot_levels(d, table)
    cont = design_currents(d, True)
    from phase29_pointB_sweep import _hpbw  # noqa: F401

    out["continuous_kappa"] = [
        field_metrics(
            cont[j],
            d["x_m"],
            d["y_m"],
            (d["ports"][j][2], d["ports"][j][3]),
            d["eta_conv"][j],
            np.arange(0.0, 90.1, 1.0),
            np.arange(0.0, 360.0, 1.5),
        )
        for j in range(4)
    ]
    out["continuous_kappa"] = [
        {k: round(v, 1) for k, v in r.items() if k in ("gain_launch_bare_dbi", "sll_1p5hpbw_db", "sll_15deg_cone_db")}
        for r in out["continuous_kappa"]
    ]
    # control: the same aperture, single feed
    from src.hologram.tensor_leaky_cp import march_rays

    k0 = t.K0
    from src.hologram import shared_aperture_multiport as sa

    mf = lambda j, f, b: march_rays(table, sa.hexagon_m(R), f, k0, b[0], b[1], 0.9, taper_c=0.5, n_psi=240)  # noqa: E731
    d1 = fr.design_frey(table, R, [(0.0, 0.0, 0.0, 0.0)], [0.0], 1000.0, 1.0, n_iter=2, n_psi=240, march_fn=mf)
    m1 = field_metrics(
        d1["J"][0],
        d1["x_m"],
        d1["y_m"],
        (0.0, 0.0),
        d1["eta_conv"][0],
        np.arange(0.0, 90.1, 1.0),
        np.arange(0.0, 360.0, 1.5),
    )
    out["single_feed_control"] = {
        k: round(v, 1)
        for k, v in m1.items()
        if k in ("gain_launch_bare_dbi", "sll_1p5hpbw_db", "sll_15deg_cone_db", "hpbw_deg")
    }
    g_ = np.array(out["ga_offsets"]["gain_dbi"])
    s_ = np.array(out["ga_offsets"]["sll_1p5hpbw_db"])
    out["rung2_pass"] = bool(np.all(g_ >= 21.0) and np.all(s_ <= -15.0))
    return out


def main():
    S = {"started": f"{dt.datetime.now():%Y-%m-%d %H:%M}"}
    _init()
    log("start")
    S["questions"] = {
        "slot_snapping": "YES: the slot (L, W) is one of the 4 table slots {(4.4,1.2), (4.6,1.6), (4.8,1.6), (4.8,2.0)} chosen per cell (as in Phase 26); the axis is continuous; alpha/beta come from the table at the nearest of 5 K columns. Continuous-kappa variants below.",
        "ga_fitness_150mm_B3": "mean port gain - 0.5 max(0, worst pointing - 2 deg) - 1.0 max(0, gain dip - 2.5 dB) - 0.15 mean max(0, |HPBW - target| - 2): SLL is NOT in the fitness (scan_300 run: + 0.3 mean max(0, SLL + 15))",
        "other_pins_in_forward_model": "NO: the model is the aperture current of the slotted cells only; the other monopole pins / coax holes appear only as a 6.5 mm cell keep-out; no pin scattering, no pin direct radiation (the Phase 26 launch fraction 0.592/0.62 is a constant).",
    }
    r1 = rung1(_S["table"])
    S["rung1"] = r1
    (RES / "B3debug_summary.json").write_text(json.dumps(S, indent=1, default=float))
    log(f"rung1 pass={r1['rung1_pass']}")
    if r1["rung1_pass"]:
        S["rung2"] = rung2(_S["table"])
        log(f"rung2 pass={S['rung2']['rung2_pass']}")
        S["rung3"] = (
            "not run by this script (needs pins in the forward model, which the repo does not have); run only if rungs 1-2 pass and the user agrees"
            if S["rung2"]["rung2_pass"]
            else "not reached"
        )
    else:
        S["rung2"] = "not reached (rung 1 failed)"
        S["rung3"] = "not reached"
    (RES / "B3debug_summary.json").write_text(json.dumps(S, indent=1, default=float))
    for a in (
        ["git", "add", "-A"],
        [
            "git",
            "commit",
            "-q",
            "-m",
            f"Phase 29 B3 debugging ladder: rung 1 (single feed vs Phase 26), rung 2 (Frey replication) (MODEL_ONLY)\n\n{ATTR}",
        ],
        ["git", "push", "hub", "phase29-seven-hex"],
    ):
        r = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
        log(f"{a[1]} rc={r.returncode}")
    (RES / "B3debug_DONE").write_text("done\n")


if __name__ == "__main__":
    main()
