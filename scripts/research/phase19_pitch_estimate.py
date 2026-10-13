"""Phase 19 D-bis: estimate (not a re-design) what a smaller cell pitch buys against the neighbour-coupling error.

    .venv/bin/python scripts/research/phase19_pitch_estimate.py

Runs ONLY the dual-beam aperture-phase synthesis on the finer grid (no cell assignment, no layout is committed) and compares
the neighbour phase steps of the target transmission phase with the current 5 mm layout, then applies the measured
coupling magnitudes (|K| of the N-differing-neighbour pools) to both.  Writes results/transmitarray_dual30/pitch_estimate.json.
"""

import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from phase19_error_bar import D_MIN, k_pools  # noqa: E402

from src.transmitarray.design import radius_for_fill  # noqa: E402
from src.transmitarray.layout import LensLayout  # noqa: E402
from src.transmitarray.pipeline import DesignConfig, synthesize  # noqa: E402
from src.transmitarray.uncertainty import neighbour_lists  # noqa: E402

RES = ROOT / "results" / "transmitarray_dual30"
C0, F = 299_792_458.0, 28e9
D_M, FOCAL, MIN_FEAT_MM = 0.150, 0.180, 0.6
COS = math.cos(math.radians(22.5))


def target_steps(pitch_m, n_cells_across):
    lay = LensLayout(D_M, D_M / n_cells_across, FOCAL)
    cfg = DesignConfig(F, lay, 30.0, 10.0, None, 0.6e-3, 230.0, 0.8044, "")
    syn = synthesize(cfg)
    ph = np.exp(1j * np.radians(syn["base_phase_deg"]))  # unit-amplitude target transmission
    x_mm, y_mm = syn["x"] * 1e3, syn["y"] * 1e3
    nb = neighbour_lists(x_mm, y_mm, pitch_m * 1e3)
    d_all = [[(ph[j] / ph[i] - 1.0, ax) for j, ax in lst] for i, lst in enumerate(nb)]
    steps = [abs(np.degrees(np.angle(1 + d))) for lst in d_all for d, _ in lst]
    return lay, d_all, np.array(steps)


def eps_rms(d_all, pools):
    med = {k: abs(complex(np.median(v.real), np.median(v.imag))) for k, v in pools.items()}
    eps = []
    for lst in d_all:
        diff = [(d, ax) for d, ax in lst if abs(d) >= D_MIN]
        n = len(diff)
        if n == 0:
            eps.append(0.0)
            continue
        axes = {ax for _d, ax in diff}
        key = (n, next(iter(axes))) if (n == 1 or (n == 2 and len(axes) == 1)) else (n, "xy")
        eps.append(med[key] * abs(np.mean([d for d, _ in diff])))
    return np.array(eps)


def main():
    samples = json.loads((RES / "coupling_samples.json").read_text())["samples"]
    sup = json.loads((RES / "coupling_superposition.json").read_text())["pairs"]
    pools = k_pools(samples, sup)
    out = {}
    for name, n_across in (("current_p5mm", 30), ("lambda_over_3", 42)):
        p = D_M / n_across
        lay, d_all, steps = target_steps(p, n_across)
        eps = eps_rms(d_all, pools)
        fill_cap = 0.804
        r_min = MIN_FEAT_MM / (2 * COS)
        r_max = (p * 1e3 - MIN_FEAT_MM) / (2 * COS)
        r_cap_mm = radius_for_fill(1 - fill_cap, p * 1e3)
        out[name] = dict(
            pitch_mm=p * 1e3,
            pitch_over_lambda=p / (C0 / F),
            cells_across=n_across,
            n_cells=int(len(d_all)),
            neighbour_target_phase_step_deg=dict(median=float(np.median(steps)), p90=float(np.percentile(steps, 90))),
            fraction_of_edges_over_90deg=float(np.mean(steps > 90)),
            octagon_R_range_mm=[r_min, r_max],
            smallest_pocket_flat_to_flat_mm_at_fill_cap=2 * r_cap_mm * COS,
            pocket_aspect_ratio_at_fill_cap=25.0 / (2 * r_cap_mm * COS),
            largest_pillar_fill=2 * math.sqrt(2) * (r_max / (p * 1e3)) ** 2,
            eps_estimate=dict(median=float(np.median(eps)), rms=float(np.sqrt(np.mean(eps**2)))),
            pitch_is_exact_decimal=bool(abs(p * 1e3 * 1e6 - round(p * 1e3 * 1e6)) < 1e-9),
        )
        print(name, json.dumps(out[name]), flush=True)
    a, b = out["current_p5mm"], out["lambda_over_3"]
    out["ratios"] = dict(
        cell_count=b["n_cells"] / a["n_cells"],
        median_step=b["neighbour_target_phase_step_deg"]["median"] / a["neighbour_target_phase_step_deg"]["median"],
        eps_rms=b["eps_estimate"]["rms"] / a["eps_estimate"]["rms"],
        eps_median=b["eps_estimate"]["median"] / max(a["eps_estimate"]["median"], 1e-9),
    )
    out["note"] = (
        "ESTIMATE only: phase synthesis on the finer grid, unit-amplitude target phases, measured |K| medians; no cell assignment or re-design."
    )
    (RES / "pitch_estimate.json").write_text(json.dumps(out, indent=1))
    print(out["ratios"])


if __name__ == "__main__":
    main()
