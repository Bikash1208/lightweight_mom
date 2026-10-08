"""Phase 19 D-bis: which perturbation model reproduces the EXACT order spectra of the real 9-cell windows?

    .venv/bin/python scripts/research/phase19_row_models.py        (after phase19_row_validation.py)

For each window of results/transmitarray_dual30/row_validation.json compare the exact orders with the isolated-cell
prediction and with draws of model R (measured magnitudes, random phase) and model S (measured magnitudes AND phase):
spectrum error ||E - E_exact|| / ||E_iso|| and total propagating power relative to the isolated prediction.
Writes results/transmitarray_dual30/row_models.json.
"""

import json
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from phase19_error_bar import D_MIN, k_pools  # noqa: E402
from phase19_local_periodicity import isolated_txx  # noqa: E402
from phase19_row_validation import ORDERS, patch_orders  # noqa: E402

OUT = ROOT / "results" / "transmitarray_dual30"
N = 600


def main():
    rv = json.loads((OUT / "row_validation.json").read_text())
    samples = json.loads((OUT / "coupling_samples.json").read_text())["samples"]
    sup = json.loads((OUT / "coupling_superposition.json").read_text())["pairs"]
    pools = k_pools(samples, sup)
    cells = sorted({tuple(c) for w in rv["windows"] for c in w["cells"]})
    with ProcessPoolExecutor(max_workers=4) as ex:
        t_iso = dict(zip(cells, ex.map(isolated_txx, [k for k, _ in cells], [r for _, r in cells]), strict=True))
    rng = np.random.default_rng(7)
    out = []
    for w in rv["windows"]:
        t = np.array([t_iso[tuple(c)] for c in w["cells"]])
        e_ex = np.array([complex(*v) for v in w["orders_exact"]])
        e_iso = patch_orders(list(t))
        norm = np.linalg.norm(e_iso)
        pw0 = np.sum(np.abs(e_iso) ** 2)
        res = {}
        for name, rand in (("R_random_phase", True), ("S_measured_phase", False)):
            errs, pws = [], []
            for _ in range(N):
                eps = np.zeros(len(t), complex)
                for i in range(len(t)):
                    d = [t[j] / t[i] - 1.0 for j in (i - 1, (i + 1) % len(t)) if abs(t[j] / t[i] - 1.0) >= D_MIN]
                    if not d:
                        continue
                    pool = pools[(len(d), "x")]
                    k = pool[rng.integers(0, len(pool))]
                    if rand:
                        k = abs(k) * np.exp(1j * rng.uniform(0, 2 * np.pi))
                    eps[i] = k * np.mean(d)
                e_m = patch_orders(list(t * (1 + eps)))
                errs.append(np.linalg.norm(e_m - e_ex) / norm)
                pws.append(np.sum(np.abs(e_m) ** 2) / pw0)
            res[name] = dict(
                spectrum_error_vs_exact=dict(
                    p10=float(np.percentile(errs, 10)),
                    p50=float(np.percentile(errs, 50)),
                    p90=float(np.percentile(errs, 90)),
                ),
                power_over_isolated=dict(
                    p10=float(np.percentile(pws, 10)),
                    p50=float(np.percentile(pws, 50)),
                    p90=float(np.percentile(pws, 90)),
                ),
            )
        res["isolated_spectrum_error_vs_exact"] = float(np.linalg.norm(e_iso - e_ex) / norm)
        res["exact_power_over_isolated"] = float(np.sum(np.abs(e_ex) ** 2) / pw0)
        out.append(res)
        print(json.dumps(res), flush=True)
    (OUT / "row_models.json").write_text(json.dumps(dict(orders=ORDERS, n_draws=N, windows=out), indent=1))


if __name__ == "__main__":
    main()
