"""Phase 29: weak-modulation synthesis test (Frey: every sub-hologram a small first-order perturbation around a common mean impedance), MODEL_ONLY. Unattended, parallel (BLAS 1 / worker).

    python scripts/research/phase29_weakmod.py all       # results/phase29/weakmod_summary.json

Modulation index of the model cell: M = |dY| / |Ybar| (orientation-dependent part of the Bloch sheet tensor over its orientation-averaged mean; `AngleLUT.M_slot`); a sub-hologram of fixed index M has target
strength `M * k_unit * g(r)` (g: the Phase 26 taper shape, max 1). The GA does NOT tune M (only spacing, ring rotation, hologram offsets). Two cell realisations: `table` = the 4 real slots (min-error assignment), `virtual` =
the tensor shape of slot 2 with a continuous strength factor (not realisable: it separates the synthesis physics from the table's resolution). Cases: 4 feeds at R = 300 mm (GA, target SLL <= -15 dB) with a single-feed control;
the two-feed 150 mm case of the HFSS package (psi2 scan, 1.35 lambda0, beams +-3.65 deg); reference = the unmodified design (no fixed M). SLL = highest RHCP lobe outside 1.5 HPBW of the peak."""

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
V_SLOT = 2
ATTR = "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01UaAn25MEtRcvrHeagn5k3z"
MS = (0.08, 0.13, 0.20)


def log(msg):
    with open(RES / "weakmod_progress.txt", "a", encoding="utf-8") as f:
        f.write(f"{dt.datetime.now():%m-%d %H:%M} {msg}\n")


def design(kind, g, cfg, n_psi=180):
    import phase29_tensor_ladder as L

    from src.hologram import shared_aperture_tensor as st

    kw = {"weak_M": cfg["M"], "virtual_slot": V_SLOT if cfg["mode"] == "virtual" else None}
    if kind == "n4":
        ports, off = L.ports4(g, 7.3)
        return st.design_tensor(L._S["table"], L._S["lut"], 300.0, ports, off, n_psi=n_psi, **kw)
    if kind == "n1_300":
        return st.design_tensor(L._S["table"], L._S["lut"], 300.0, [(0.0, 0.0, 0.0, 0.0)], [0.0], n_psi=n_psi, **kw)
    if kind == "n1_150":
        return st.design_tensor(L._S["table"], L._S["lut"], 150.0, [(0.0, 0.0, 0.0, 0.0)], [0.0], n_psi=n_psi, **kw)
    s, tilt = 1.35 * LAM, 7.3 / 2  # two-feed 150 mm (HFSS package geometry)
    return st.design_tensor(
        L._S["table"],
        L._S["lut"],
        150.0,
        [(s / 2, 0.0, tilt, 0.0), (-s / 2, 0.0, tilt, 180.0)],
        [0.0, float(g[0])],
        n_psi=n_psi,
        **kw,
    )


def evaluate(kind, g, cfg, grid=(1.0, 1.5), n_psi=180, full=False):
    import phase29_tensor_ladder as L
    from phase29_pointB_sweep import port_metrics

    from src.hologram import shared_aperture_multiport as sa

    d = design(kind, g, cfg, n_psi)
    th, ph = np.arange(0.0, 90.0 + 1e-9, grid[0]), np.arange(0.0, 360.0, grid[1])
    F = sa.port_fields(d, L._S["green"], th, ph)
    pm = port_metrics(F, th, ph, d["beams"])
    gain = float(np.mean([p["peak_gain_dbi"] for p in pm]))
    sll = float(np.mean([p["sll_db"] for p in pm]))
    out = {"score": gain - 0.3 * max(0.0, sll + 20.0), "mean_gain": gain, "mean_sll": sll}
    if full:
        out.update(
            {
                "gain": [round(p["peak_gain_dbi"], 1) for p in pm],
                "sll": [round(p["sll_db"], 1) for p in pm],
                "hpbw": [round(p["hpbw_deg"], 1) for p in pm],
                "eta_conv": [round(e, 3) for e in d["eta_conv"]],
                "index_rms": [round(x, 3) for x in d["index_rms"]],
                "slot_usage": np.bincount(d["slot"], minlength=4).tolist(),
                "n_cells": int(d["x_m"].size),
            }
        )
    return out


def _fit(a):
    kind, g, cfg = a
    try:
        return evaluate(kind, np.asarray(g), cfg)["score"]
    except Exception:  # an infeasible individual is simply bad
        return -99.0


def ga4(cfg, pop_n=32, gens=8):
    import phase29_tensor_ladder as L
    from phase29_pointB_sweep import n_workers

    rng = np.random.default_rng(1)
    lo, hi = np.array([1.35 * LAM, -30.0] + [0.0] * 4), np.array([1.50 * LAM, 30.0] + [360.0] * 4)
    pop = lo + rng.random((pop_n, 6)) * (hi - lo)
    pop[0] = [1.35 * LAM, 0.0, 0.0, 0.0, 0.0, 0.0]
    with mp.get_context("spawn").Pool(n_workers(per_worker_gb=1.5, reserve_gb=8.0), initializer=L.init) as pool:
        score = np.array(pool.map(_fit, [("n4", p, cfg) for p in pop]))
        for _ in range(gens):
            order = np.argsort(-score)
            new = [pop[order[0]], pop[order[1]]]
            while len(new) < pop_n:
                i, j = (min(rng.integers(0, pop_n, 3), key=lambda k: -score[k]) for _ in range(2))
                ch = (
                    pop[i]
                    + (rng.random(6) * 1.6 - 0.3) * (pop[j] - pop[i])
                    + (rng.random(6) < 0.15) * rng.normal(size=6) * 0.1 * (hi - lo)
                )
                new.append(np.clip(ch, lo, hi))
            pop = np.array(new)
            score = np.concatenate([score[order[:2]], np.array(pool.map(_fit, [("n4", p, cfg) for p in pop[2:]]))])
    return pop[int(np.argmax(score))]


def scan2(cfg):
    import phase29_tensor_ladder as L
    from phase29_pointB_sweep import n_workers

    psis = list(range(0, 360, 30))
    with mp.get_context("spawn").Pool(min(12, n_workers(per_worker_gb=1.5)), initializer=L.init) as pool:
        sc = pool.map(_fit, [("n2", [p], cfg) for p in psis])
    return float(psis[int(np.argmax(sc))])


def main():
    import phase29_tensor_ladder as L

    S = {"started": f"{dt.datetime.now():%Y-%m-%d %H:%M}"}
    L.init()
    lut, table = L._S["lut"], L._S["table"]
    kappa = lut.kappa_old
    S["table_resolution"] = {
        "M_slot": [round(float(x), 3) for x in lut.M_slot],
        "kappa_slot": [round(float(x), 2) for x in kappa],
        "alpha_slot_np_per_m": [round(float(np.nanmean(table.alpha[i])), 2) for i in range(4)],
        "k_unit": round(lut.k_unit, 1),
        "kappa_target_for_M": {str(m): round(m * lut.k_unit, 2) for m in MS},
        "alpha_target_np_per_m_for_M": {str(m): round(0.5 * (m * lut.k_unit) ** 2, 3) for m in MS},
        "N4_total_modulation_sqrtN_M_and_N_M": {str(m): [round(2 * m, 2), round(4 * m, 2)] for m in MS},
        "table_M_max": round(float(lut.M_slot.max()), 3),
    }
    log("start")
    cfgs = [{"name": "reference", "mode": "table", "M": None}] + [
        {"name": f"{mode}_M{m}", "mode": mode, "M": m} for mode in ("table", "virtual") for m in MS
    ]
    rows = {}
    for cfg in cfgs:
        t0 = time.time()
        g = ga4(cfg)
        r = {
            "n4_300mm": evaluate("n4", g, cfg, grid=(1.0, 1.0), n_psi=240, full=True),
            "n1_300mm": evaluate("n1_300", np.zeros(1), cfg, grid=(1.0, 1.0), n_psi=240, full=True),
        }
        psi2 = scan2(cfg)
        r["n2_150mm"] = evaluate("n2", [psi2], cfg, grid=(1.0, 1.0), n_psi=240, full=True) | {"psi2_best_deg": psi2}
        r["n1_150mm"] = evaluate("n1_150", np.zeros(1), cfg, grid=(1.0, 1.0), n_psi=240, full=True)
        rows[cfg["name"]] = r
        log(
            f"{cfg['name']}: n4 gain {r['n4_300mm']['mean_gain']:.1f} SLL {r['n4_300mm']['mean_sll']:.1f}; n2 gain {r['n2_150mm']['mean_gain']:.1f} SLL {r['n2_150mm']['mean_sll']:.1f} [{time.time() - t0:.0f}s]"
        )
        S["rows"] = rows
        (RES / "weakmod_summary.json").write_text(json.dumps(S, indent=1, default=float))
    p = ROOT / "docs/STATUS.md"
    s = p.read_text(encoding="latin-1")
    i = s.index("## Open next steps")
    row = "\n## Phase 29 weak-modulation synthesis test (2026-10-28) - MODEL_ONLY\n\nFixed modulation index M = 0.08 / 0.13 / 0.20 per sub-hologram, real slot table vs virtual continuous cell, 4 feeds 300 mm and the two-feed 150 mm case: `results/phase29/weakmod_summary.json`.\n\n"
    p.write_text(s[:i] + row + s[i:], encoding="latin-1", newline="")
    for a in (
        ["git", "add", "-A"],
        [
            "git",
            "commit",
            "-q",
            "-m",
            f"Phase 29 weak-modulation synthesis test: fixed modulation index, table vs virtual cell, 4 feeds 300 mm + two-feed 150 mm (MODEL_ONLY)\n\n{ATTR}",
        ],
        ["git", "push", "hub", "phase29-seven-hex"],
    ):
        r = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
        log(f"{a[1]} rc={r.returncode}")
    (RES / "weakmod_DONE").write_text("done\n")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "all":
        main()
