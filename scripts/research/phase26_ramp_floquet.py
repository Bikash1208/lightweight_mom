"""Phase 26 B2.4 (corrected gate): supercell ramp vs the FLOQUET solution of the modulated tensor sheet.

    python scripts/research/phase26_ramp_floquet.py <out.json> [--mesh 24] [--design D,L,W ...] [--n-angle 8]

The B2.4 locality gate compared the supercell `beta` with the mean of the local `beta_i`; a periodically
modulated surface has no reason to obey that (second-order shift in the modulation depth). The right
reference is the harmonic balance of the staircase sheet built from the LOCAL tensors (`tensor_floquet`).
Contrast numbers kept: local mean, effective medium (mean tensor).

Cell tensors are taken at one transverse wavenumber for every cell: `kt = beta_super` along x ("at the
operating wavenumber", the case a database tabulated vs kt is interpolated at) and, for contrast, each cell
at its own Bloch root `beta_i` ("own"). Gate: |beta_floquet - beta_super| <= 0.01 k0; alpha agreement is
reported (supercell alpha: one Newton step of the complex MoM tensor determinant, radiation into n = -1 and
the dielectric loss included -- an estimate, not a converged complex root).
Labels: SOLVER_INTERNAL, MODEL_ONLY, NOT_HFSS_VALIDATED.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase26_ramp_locality as rl  # noqa: E402
from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.sheet_model import _y_air, _y_slab  # noqa: E402
from src.dispersion.slotted_cells import CellRow, SlottedCircularPatch  # noqa: E402
from src.dispersion.tensor_floquet import floquet_complex_beta_tensor  # noqa: E402
from src.frequency import angular_frequency  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

F, P, SUB, K0 = rl.F, rl.P, rl.SUB, rl.K0
GATE_K0 = 0.01


def supercell_alpha(cells, mesh_n, beta_r):
    """One complex Newton step of det(Y_s + diag(A, B)) at the real constant `beta_r`; returns
    `(beta_complex, y_supercell)`."""
    n = len(cells)
    uc = UnitCell(RectangularLattice(n * P, P), CellRow(tuple(cells), P), SUB, Metal(), backing="pec")
    mesh = MeshSettings(n * mesh_n, mesh_n, max(n * mesh_n, mesh_n))
    basis = RooftopBasis(uc, mesh.cells_x, mesh.cells_y)
    omega = angular_frequency(F)

    def f(beta):
        y = ts.bloch_sheet_tensor(uc, F, mesh, beta, 0.0, basis=basis, orders=(n * mesh_n, mesh_n))
        a = _y_air(beta, omega, K0) + _y_slab(beta, omega, SUB.eps_complex, K0, SUB.thickness)
        b = ts._b_te(beta, omega, SUB.eps_complex, K0, SUB.thickness)
        return complex(np.linalg.det(y + np.array([[a, 0.0], [0.0, b]])))

    h = 0.003 * K0
    f0 = f(beta_r)
    fp = (f(beta_r + h) - f(beta_r - h)) / (2 * h)
    return complex(beta_r - f0 / fp)


def one(name, params, mesh_n):
    t0 = time.time()
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    geos = [SlottedCircularPatch(d * 1e-3, l * 1e-3, w * 1e-3, np.radians(a)) for d, l, w, a in params]
    cells = [rl.single_cell(*p) for p in params]
    n = len(params)
    betas = np.array([rl.cell_beta(*p, mesh_n) for p in params])
    b_super, _ = rl.supercell_beta(geos, mesh_n, float(betas.mean()))
    beta_c = supercell_alpha(geos, mesh_n, b_super)
    print(f"[{name}] local mean {betas.mean() / K0:.4f}  super {b_super / K0:.4f}  alpha_super~{-beta_c.imag:.2f} Np/m ({time.time() - t0:.0f}s)", flush=True)

    def tensors_at(kts):
        return np.array([ts.bloch_sheet_tensor(c, F, mesh, kt, 0.0) for c, kt in zip(cells, kts)])

    out = {
        "n_cells": n,
        "params_D_L_W_angle": params,
        "beta_local_over_k0": (betas / K0).tolist(),
        "beta_local_mean_over_k0": float(betas.mean() / K0),
        "beta_super_over_k0": b_super / K0,
        "alpha_super_estimate_np_per_m": float(-beta_c.imag),
        "beta_super_complex_estimate_over_k0_real": float(beta_c.real / K0),
        "gate_k0": GATE_K0,
        "models": {},
    }
    for tag, kts in (("at_operating_beta", [b_super] * n), ("at_own_beta", betas.tolist())):
        tens = tensors_at(kts)
        ym = tens.mean(axis=0)
        b_eff, _ = ts.find_beta_tensor(ym, 0.0, SUB, F)
        rec = {"effective_medium_beta_over_k0": b_eff / K0}
        for nm in (6, 10):
            bf, its = floquet_complex_beta_tensor(tens, n * P, SUB, F, b_eff, n_max=nm)
            rec[f"floquet_nmax{nm}"] = {
                "beta_over_k0": float(bf.real / K0),
                "alpha_np_per_m": float(-bf.imag),
                "beta_error_vs_super_k0": float((bf.real - b_super) / K0),
                "iterations": its,
            }
        e = abs(rec["floquet_nmax10"]["beta_error_vs_super_k0"])
        rec["gate"] = "PASS" if e <= GATE_K0 else "FAIL"
        out["models"][tag] = rec
        print(f"[{name}] {tag}: eff-med {b_eff / K0:.4f}  floquet {rec['floquet_nmax10']['beta_over_k0']:.4f} "
              f"(err {rec['floquet_nmax10']['beta_error_vs_super_k0']:+.4f} k0, alpha {rec['floquet_nmax10']['alpha_np_per_m']:.2f}) -> {rec['gate']}", flush=True)
    out["mean_of_local_error_k0"] = float((betas.mean() - b_super) / K0)
    out["seconds"] = time.time() - t0
    return out


def main(out_path, mesh_n, designs, n_angle):
    res = {"labels": ["SOLVER_INTERNAL", "MODEL_ONLY", "NOT_HFSS_VALIDATED"], "mesh_per_cell": mesh_n, "ramps": {}}
    for d, l, w in designs:
        params = [(d, l, w, 180.0 * i / n_angle) for i in range(n_angle)]
        res["ramps"][f"D{d:g}_{l:g}x{w:g}_slot_angle"] = one(f"D{d:g} {l:g}x{w:g}", params, mesh_n)
        Path(out_path).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 24
    designs = [tuple(float(v) for v in a[i + 1].split(",")) for i, s in enumerate(a) if s == "--design"] or [(6.0, 4.8, 2.0)]
    n_angle = int(a[a.index("--n-angle") + 1]) if "--n-angle" in a else 8
    main(a[0], mesh, designs, n_angle)
