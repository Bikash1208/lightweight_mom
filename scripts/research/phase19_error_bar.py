"""Phase 19 D-bis: Monte Carlo error bar on the lens prediction from the measured neighbour-coupling error.

    .venv/bin/python scripts/research/phase19_error_bar.py

Needs results/transmitarray_dual30/{layout,summary,coupling_samples,coupling_superposition}.json.  Every scenario perturbs the
cell's co-pol exit field, E_c' = E_c (1 + eps_c), with d_n = T_n/T_c - 1 the isolated-response contrast to neighbour n:
  S   PRIMARY.  eps_c = K * mean(d_n) over the neighbours with |d_n| >= D_MIN (N of them), K bootstrapped as a COMPLEX number
      (measured phase kept) from exact-supercell samples with the same number N of differing neighbours:
      N=1 strips [A A B B]; N=2 alternating strips [A B]; N=4 2x2 checkerboards; N=3 pools N=2 and N=4.  x-neighbours use the
      E-normal samples, y-neighbours the E-parallel ones (the horn is x-polarised).  The supercell data constrain the phase
      (single-neighbour K is not a random phase), so it is used.
  R   SAME magnitudes as S, phase uniform random (the assumption the task asks to state when the data does not constrain it).
  R0  random phase with the earlier ORDER-0 relation (0.07@60, 0.2@110, 0.5@170 deg) used as a per-cell magnitude.
  Cn  NAIVE linear superposition of single-neighbour coupling over all neighbours -- refuted by the 2- and 4-neighbour supercells
      (it over-predicts), kept only to show how much that assumption matters.
Writes results/transmitarray_dual30/error_bar.json.
"""

import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.transmitarray import pattern as pt  # noqa: E402
from src.transmitarray.layout import LensLayout  # noqa: E402
from src.transmitarray.uncertainty import PatternOperator, neighbour_lists, percentiles  # noqa: E402

RES = ROOT / "results" / "transmitarray_dual30"
C0 = 299_792_458.0
N_DRAWS = 2000
SEED = 20261008
D_MIN = 0.6  # |T_n/T_c - 1| counted as a 'differing' neighbour (~35 deg at equal amplitude)
STEP_MIN = 45.0  # supercell samples used for K: only pairs with a real contrast, where d is not tiny


def load():
    layout = json.loads((RES / "layout.json").read_text())["cells"]
    summ = json.loads((RES / "summary.json").read_text())
    samples = json.loads((RES / "coupling_samples.json").read_text())["samples"]
    sup = json.loads((RES / "coupling_superposition.json").read_text())["pairs"]
    return layout, summ, samples, sup


def fields(layout, summ):
    inp = summ["inputs"]
    lay = LensLayout(inp["d_flats_mm"] * 1e-3, inp["pitch_mm"] * 1e-3, inp["focal_mm"] * 1e-3)
    x = np.array([c["x_mm"] for c in layout]) * 1e-3
    y = np.array([c["y_mm"] for c in layout]) * 1e-3
    th, _ph, r = lay.incidence(x, y)
    k0 = 2 * math.pi * inp["freq_hz"] / C0
    jl = [
        np.array([[complex(*c["J"][0]), complex(*c["J"][1])], [complex(*c["J"][2]), complex(*c["J"][3])]])
        for c in layout
    ]
    ex, ey, tau = pt.exit_field(lay, inp["q"], x, y, th, r, k0, jl)
    t_co = np.array([j[0, 0] for j in jl])
    return dict(
        lay=lay, x=x, y=y, th=th, r=r, k0=k0, ex=ex, ey=ey, tau=tau, t_co=t_co, q=inp["q"], pitch=inp["pitch_mm"] * 1e-3
    )


def median_eps_relation(samples, axis):
    """Measured per-cell |eps| vs |step|: bin medians -> interpolation table (anchored at 0)."""
    s = [(abs(q["step_deg"]), q["eps_abs"]) for q in samples if q["pol"] == axis]
    edges = [0, 30, 60, 90, 120, 150, 181]
    cx, cy = [0.0], [0.0]
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        v = [e for st, e in s if lo <= st < hi]
        if v:
            cx.append(float(np.mean([st for st, _ in s if lo <= st < hi])))
            cy.append(float(np.median(v)))
    return np.array(cx), np.array(cy)


def build_edges(layout, pitch_mm, t_co):
    nb = neighbour_lists([c["x_mm"] for c in layout], [c["y_mm"] for c in layout], pitch_mm)
    ci, cj, ax, ratio, step = [], [], [], [], []
    for i, lst in enumerate(nb):
        for j, axis in lst:
            ci.append(i)
            cj.append(j)
            ax.append(0 if axis == "x" else 1)
            ratio.append(t_co[j] / t_co[i] - 1.0)
            step.append(abs(float(np.degrees(np.angle(t_co[j] / t_co[i])))))
    return np.array(ci), np.array(cj), np.array(ax), np.array(ratio), np.array(step)


def k_pools(samples, sup):
    """Complex K = eps/d pools keyed by (N differing neighbours, axis)."""
    pools = {
        (1, "x"): [complex(*q["kappa"]) for q in samples if q["pol"] == "x_normal" and abs(q["step_deg"]) >= STEP_MIN],
        (1, "y"): [
            complex(*q["kappa"]) for q in samples if q["pol"] == "y_parallel" and abs(q["step_deg"]) >= STEP_MIN
        ],
        (2, "x"): [],
        (2, "y"): [],
        (4, "xy"): [],
    }
    for pair in sup:
        if abs(pair["step_deg"]) < STEP_MIN:
            continue
        for key, dest in (("alt:x_normal", (2, "x")), ("alt:y_parallel", (2, "y")), ("chk:x_normal", (4, "xy"))):
            for c in pair["cases"].get(key, []):
                pools[dest].append(complex(*c["kappa_eff_per_neighbour"]) * c["n_diff"])
    pools = {k: np.array(v) for k, v in pools.items()}
    pools[(3, "xy")] = np.concatenate([pools[(2, "x")], pools[(2, "y")], pools[(4, "xy")]])
    pools[(2, "xy")] = np.concatenate([pools[(2, "x")], pools[(2, "y")]])
    return pools


def cell_structure(layout, pitch_mm, t_co):
    """Per cell: list of (d_n, axis) over its neighbours."""
    nb = neighbour_lists([c["x_mm"] for c in layout], [c["y_mm"] for c in layout], pitch_mm)
    return [[(t_co[j] / t_co[i] - 1.0, ax) for j, ax in lst] for i, lst in enumerate(nb)]


def eps_model_S(rng, struct, pools, n, random_phase=False):
    eps = np.zeros((len(struct), n), complex)
    for i, nbs in enumerate(struct):
        diff = [(d, ax) for d, ax in nbs if abs(d) >= D_MIN]
        N = len(diff)
        if N == 0:
            continue
        axes = {ax for _d, ax in diff}
        key = (N, next(iter(axes))) if N == 1 else (N, next(iter(axes))) if (N == 2 and len(axes) == 1) else (N, "xy")
        pool = pools[key]
        kk = pool[rng.integers(0, len(pool), n)]
        if random_phase:
            kk = np.abs(kk) * np.exp(1j * rng.uniform(0, 2 * np.pi, n))
        eps[i] = kk * np.mean([d for d, _ax in diff])
    return eps


def eps_model_naive(rng, struct, pools, n):
    eps = np.zeros((len(struct), n), complex)
    for i, nbs in enumerate(struct):
        for d, ax in nbs:
            pool = pools[(1, ax)]
            eps[i] += pool[rng.integers(0, len(pool), n)] * d
    return eps


def eps_model_R0(rng, struct, n):
    var = np.array(
        [
            sum(
                np.interp(
                    abs(np.degrees(np.angle(1 + d))), [0, 21, 60, 110, 170, 180], [0, 0.012, 0.07, 0.21, 0.5, 0.5]
                )
                ** 2
                for d, _ax in nbs
            )
            for nbs in struct
        ]
    )
    return np.sqrt(var)[:, None] * np.exp(1j * rng.uniform(0, 2 * np.pi, (len(struct), n)))


def summarize(res):
    out = {
        "mean_gain_dbi": percentiles(res["mean_gain_dbi"]),
        "gain_beam_phi0_dbi": percentiles(res["gain_dbi"][:, 0]),
        "gain_beam_phi180_dbi": percentiles(res["gain_dbi"][:, 1]),
        "beam_theta_phi0_deg": percentiles(res["theta_deg"][:, 0]),
        "beam_theta_phi180_deg": percentiles(res["theta_deg"][:, 1]),
        "pointing_error_deg_worst_beam": percentiles(res["pointing_error_deg"].max(axis=1)),
        "sll_db": percentiles(res["sll_db"]),
        "level_difference_db_signed": percentiles(res["level_diff_db"]),
        "level_difference_db_abs": percentiles(np.abs(res["level_diff_db"])),
    }
    return out


def stability(res, key_fn, sizes):
    v = key_fn(res)
    return {str(k): percentiles(v[:k]) for k in sizes if k <= len(v)}


def main():
    layout, summ, samples, sup = load()
    f = fields(layout, summ)
    uinc = (f["x"] / f["r"], f["y"] / f["r"])
    op = PatternOperator(f["x"], f["y"], f["k0"], f["pitch"], uinc, summ["inputs"]["beam_deg"])
    nominal = op.run(f["ex"][:, None])
    committed = summ[summ["chosen"]]["pattern"]
    pools = k_pools(samples, sup)
    struct = cell_structure(layout, summ["inputs"]["pitch_mm"], f["t_co"])
    n_diff = np.array([sum(abs(d) >= D_MIN for d, _ax in nbs) for nbs in struct])
    edge_steps = np.array([abs(float(np.degrees(np.angle(1 + d)))) for nbs in struct for d, _ax in nbs])
    out = dict(
        n_draws=N_DRAWS,
        seed=SEED,
        d_min=D_MIN,
        step_min=STEP_MIN,
        labels={
            "nominal": "ISOLATED-CELL deterministic prediction (the committed numbers), recomputed on the Monte Carlo operator",
            "level_difference_identity": "0.00 dB nominal inter-beam level difference is an IDENTITY of the mirror-symmetric layout, not a performance result",
        },
        nominal_isolated_cell={
            "mean_gain_dbi": float(nominal["mean_gain_dbi"][0]),
            "gain_dbi": nominal["gain_dbi"][0].tolist(),
            "beam_theta_deg": nominal["theta_deg"][0].tolist(),
            "sll_db": float(nominal["sll_db"][0]),
            "level_difference_db": float(nominal["level_diff_db"][0]),
            "committed_fine_grid": dict(
                peak_gain_dbi=committed["beam_phi0"]["peak_gain_dbi"], sll_db=committed["sll_db"]
            ),
        },
        edge_phase_step_deg=dict(
            median=float(np.median(edge_steps)), p90=float(np.percentile(edge_steps, 90)), n_edges=int(len(edge_steps))
        ),
        differing_neighbour_count_histogram={str(k): int(np.sum(n_diff == k)) for k in range(5)},
        k_pool_sizes={f"{k[0]}:{k[1]}": int(len(v)) for k, v in pools.items()},
        k_pool_median={
            f"{k[0]}:{k[1]}": [
                float(np.abs(np.median(v.real) + 1j * np.median(v.imag))),
                float(np.degrees(np.angle(np.median(v.real) + 1j * np.median(v.imag)))),
            ]
            for k, v in pools.items()
        },
        scenarios={},
    )
    rng = np.random.default_rng(SEED)
    for sc in ("S", "R", "R0", "Cn"):
        if sc == "S":
            eps = eps_model_S(rng, struct, pools, N_DRAWS)
        elif sc == "R":
            eps = eps_model_S(rng, struct, pools, N_DRAWS, random_phase=True)
        elif sc == "R0":
            eps = eps_model_R0(rng, struct, N_DRAWS)
        else:
            eps = eps_model_naive(rng, struct, pools, N_DRAWS)
        res = op.run(f["ex"][:, None] * (1.0 + eps))
        out["scenarios"][sc] = dict(
            summary=summarize(res),
            eps_rms=float(np.sqrt(np.mean(np.abs(eps) ** 2))),
            eps_abs_median=float(np.median(np.abs(eps))),
            stability={
                "mean_gain_dbi": stability(res, lambda r: r["mean_gain_dbi"], (250, 500, 1000, 2000)),
                "sll_db": stability(res, lambda r: r["sll_db"], (250, 500, 1000, 2000)),
                "level_difference_db_abs": stability(res, lambda r: np.abs(r["level_diff_db"]), (250, 500, 1000, 2000)),
                "pointing_error_deg_worst_beam": stability(
                    res, lambda r: r["pointing_error_deg"].max(axis=1), (250, 500, 1000, 2000)
                ),
            },
        )
        print(
            sc,
            json.dumps(out["scenarios"][sc]["summary"]["mean_gain_dbi"]),
            json.dumps(out["scenarios"][sc]["summary"]["sll_db"]),
            flush=True,
        )
    (RES / "error_bar.json").write_text(json.dumps(out, indent=1))
    print("nominal", out["nominal_isolated_cell"])


if __name__ == "__main__":
    main()
