"""Phase 26 B2.4d: B2.4 Floquet gate with HARMONIC-DEPENDENT cell tensors Y_i(kt = beta + (n+q)K/2).

    python scripts/research/phase26_ramp_floquet_dispersive.py <out.json> [--mesh 24] [--n-max 3] [--design D,L,W ...]

Same three ramps as `phase26_ramp_floquet.py`; the coupling of harmonics n <- q uses the Bloch tensors of the
cells at the mean wavenumber of the pair. Harmonics beyond |kt| = pi/P are aliased in the single-cell Bloch
extraction, hence the modest n_max. Labels: SOLVER_INTERNAL, MODEL_ONLY.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase26_ramp_floquet as rf  # noqa: E402
import phase26_ramp_locality as rl  # noqa: E402
from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.slotted_cells import SlottedCircularPatch  # noqa: E402
from src.dispersion.tensor_floquet import (  # noqa: E402
    floquet_complex_beta_dispersive,
    floquet_complex_beta_tensor,
    staircase_fourier_tensor,
)

F, P, SUB, K0, GATE = rl.F, rl.P, rl.SUB, rl.K0, 0.01


def one(name, params, mesh_n, n_max, b_super=None):
    t0 = time.time()
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    cells = [rl.single_cell(*p) for p in params]
    bases = [RooftopBasis(c, mesh.cells_x, mesh.cells_y) for c in cells]
    n = len(params)
    geos = [SlottedCircularPatch(d * 1e-3, l * 1e-3, w * 1e-3, np.radians(a)) for d, l, w, a in params]
    betas = np.array([rl.cell_beta(*p, mesh_n) for p in params])
    if b_super is None:
        b_super, _ = rl.supercell_beta(geos, mesh_n, float(betas.mean()))
    K = 2.0 * np.pi / (n * P)
    y_by_m = {}
    for s in range(-2 * n_max, 2 * n_max + 1):
        kx = b_super + 0.5 * s * K
        tens = np.array([ts.bloch_sheet_tensor(c, F, mesh, kx, 0.0, basis=b) for c, b in zip(cells, bases)])
        y_by_m[s] = staircase_fourier_tensor(tens, n_max)
    ym0 = y_by_m[0][2 * n_max]
    b_eff, _ = ts.find_beta_tensor(ym0, 0.0, SUB, F)
    tens0 = np.array([ts.bloch_sheet_tensor(c, F, mesh, b_super, 0.0, basis=b) for c, b in zip(cells, bases)])
    b_const, _ = floquet_complex_beta_tensor(tens0, n * P, SUB, F, b_eff, n_max=n_max)
    b_disp, its = floquet_complex_beta_dispersive(y_by_m, n * P, SUB, F, b_eff, n_max)
    rec = {
        "params_D_L_W_angle": params,
        "beta_super_over_k0": b_super / K0,
        "n_max": n_max,
        "effective_medium_over_k0": b_eff / K0,
        "floquet_constant_tensor": {"beta_over_k0": b_const.real / K0, "alpha_np_per_m": -b_const.imag,
                                    "error_k0": (b_const.real - b_super) / K0},
        "floquet_harmonic_dependent": {"beta_over_k0": b_disp.real / K0, "alpha_np_per_m": -b_disp.imag,
                                       "error_k0": (b_disp.real - b_super) / K0, "iterations": its},
        "gate_k0": GATE,
        "seconds": time.time() - t0,
    }
    rec["gate"] = "PASS" if abs(rec["floquet_harmonic_dependent"]["error_k0"]) <= GATE else "FAIL"
    print(f"[{name}] super {b_super / K0:.4f}  const-tensor {b_const.real / K0:.4f}  dispersive {b_disp.real / K0:.4f} "
          f"(err {rec['floquet_harmonic_dependent']['error_k0']:+.4f}, alpha {-b_disp.imag:.2f}) -> {rec['gate']}  ({rec['seconds']:.0f}s)", flush=True)
    return rec


def main(out_path, mesh_n, designs, n_max, n_angle=8):
    res = {"labels": ["SOLVER_INTERNAL", "MODEL_ONLY", "NOT_HFSS_VALIDATED"], "mesh_per_cell": mesh_n, "ramps": {}}
    for d, l, w in designs:
        params = [(d, l, w, 180.0 * i / n_angle) for i in range(n_angle)]
        res["ramps"][f"D{d:g}_{l:g}x{w:g}"] = one(f"D{d:g} {l:g}x{w:g}", params, mesh_n, n_max)
        Path(out_path).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 24
    n_max = int(a[a.index("--n-max") + 1]) if "--n-max" in a else 3
    designs = [tuple(float(v) for v in a[i + 1].split(",")) for i, s in enumerate(a) if s == "--design"] or [(6.0, 4.8, 2.0)]
    main(a[0], mesh, designs, n_max)
