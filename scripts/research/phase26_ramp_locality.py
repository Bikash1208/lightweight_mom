"""Phase 26 B2.4: locality check for CONTINUOUS modulation -- supercell ramps vs the adiabatic (local) prediction.

    python scripts/research/phase26_ramp_locality.py <out.json> [--mesh 24] [--design D,L,W] [--sizes Dmin,Dmax] [--quick]

Two ramps of the slotted-circle family, a surface wave travelling along x (the ramp direction):
  * `slot_angle`: 8 cells, slot angle 0 -> 157.5 deg in 22.5 deg steps at fixed size (0 -> 180 deg is exactly
    periodic: the cell at 180 deg is the cell at 0 deg);
  * `size`: 10 cells, circle diameter D_min -> D_max -> D_min (triangular, periodic) at fixed slot angle.
The supercell is ONE periodic MoM (period N P x P, `tensor_sheet.bloch_sheet_tensor` with the supercell's own
Floquet orders); its TM0-like constant `beta_super` is the root of the supercell's mean-sheet dispersion at
the Bloch wavenumber `kx0 = beta`. The local prediction is the adiabatic one: every cell's own `beta_i`
(single-cell Bloch solve, same mesh, propagation along x) and the accumulated phase `sum beta_i P`. Gate:
|beta_super N P - sum beta_i P| <= 15 deg.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.slotted_cells import CellRow, SlottedCircularPatch  # noqa: E402
from src.frequency import angular_frequency, free_space_wavenumber  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

F, P = 10.54e9, 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
K0 = free_space_wavenumber(F)
GATE_DEG = 15.0


def single_cell(d, length, width, angle_deg):
    geo = SlottedCircularPatch(d * 1e-3, length * 1e-3, width * 1e-3, np.radians(angle_deg))
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")


def cell_beta(d, length, width, angle_deg, mesh_n):
    b, _ = ts.bloch_beta(single_cell(d, length, width, angle_deg), F, MeshSettings(mesh_n, mesh_n, mesh_n), 0.0, n_grid=5, refine=1)
    return b


def supercell_beta(cells, mesh_n, beta_guess):
    """Root of the supercell's mean-sheet dispersion along x by secant on 3 Bloch wavenumbers."""
    n = len(cells)
    lat = RectangularLattice(n * P, P)
    uc = UnitCell(lat, CellRow(tuple(cells), P), SUB, Metal(), backing="pec")
    mesh = MeshSettings(n * mesh_n, mesh_n, max(n * mesh_n, mesh_n))
    basis = RooftopBasis(uc, mesh.cells_x, mesh.cells_y)
    omega = angular_frequency(F)
    orders = (n * mesh_n, mesh_n)

    def residual(beta):
        from src.dispersion.sheet_model import _y_air, _y_slab

        y = ts.bloch_sheet_tensor(uc, F, mesh, beta, 0.0, basis=basis, orders=orders)
        a = (_y_air(beta, omega, K0) + _y_slab(beta, omega, SUB.eps_complex, K0, SUB.thickness)).imag
        b = ts._b_te(beta, omega, SUB.eps_complex, K0, SUB.thickness).imag
        m = y.imag + np.array([[a, 0.0], [0.0, b]])
        return float(np.linalg.det(m)), y

    b0 = beta_guess
    b1 = beta_guess * 1.004
    f0, _ = residual(b0)
    f1, y1 = residual(b1)
    root = None
    for _ in range(4):
        if f1 == f0:
            break
        b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
        b0, f0 = b1, f1
        b1 = b2
        f1, y1 = residual(b1)
        if abs(b1 - b0) < 2e-4 * K0:
            root = b1
            break
    return float(b1 if root is None else root), y1


def run_ramp(name, params, mesh_n):
    t0 = time.time()
    betas = [cell_beta(*p, mesh_n) for p in params]  # a missing root raises: reported, never filled in
    print(f"[{name}] local betas/k0: {np.round(np.array(betas) / K0, 4).tolist()}  ({time.time() - t0:.0f}s)", flush=True)
    cells = [SlottedCircularPatch(d * 1e-3, l * 1e-3, w * 1e-3, np.radians(a)) for d, l, w, a in params]
    guess = float(np.mean(betas))
    b_super, _ = supercell_beta(cells, mesh_n, guess)
    b_check, _ = supercell_beta(cells, mesh_n, guess * 0.985)  # a second start: the root must not depend on it
    n = len(params)
    phase_local = float(np.degrees(np.sum(betas) * P))
    phase_super = float(np.degrees(b_super * n * P))
    err = phase_super - phase_local
    out = {
        "n_cells": n,
        "params_D_L_W_angle": params,
        "beta_local_over_k0": [b / K0 for b in betas],
        "beta_super_over_k0": b_super / K0,
        "beta_super_second_start_over_k0": b_check / K0,
        "mean_local_beta_over_k0": guess / K0,
        "phase_local_deg": phase_local,
        "phase_super_deg": phase_super,
        "accumulated_phase_error_deg": err,
        "beta_spread_over_k0": float((max(betas) - min(betas)) / K0),
        "gate_deg": GATE_DEG,
        "gate": "PASS" if abs(err) <= GATE_DEG else "FAIL",
        "seconds": time.time() - t0,
    }
    print(f"[{name}] super beta/k0 {b_super / K0:.4f}  mean local {guess / K0:.4f}  phase error {err:+.1f} deg  -> {out['gate']}", flush=True)
    return out


def main(out_path, mesh_n=24, quick=False, design=(5.6, 4.4, 1.2), size_range=(5.2, 5.8), n_angle=None, only_angle=False):
    d0, length, width = design
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "mesh_per_cell": mesh_n, "frequency_ghz": F / 1e9, "ramps": {}}
    n_a = n_angle or (4 if quick else 8)
    angle = [(d0, length, width, 180.0 * i / n_a) for i in range(n_a)]
    res["ramps"]["slot_angle_0_to_180"] = run_ramp("slot_angle", angle, mesh_n)
    Path(out_path).write_text(json.dumps(res, indent=1))
    if only_angle:
        res["design_D_L_W_mm"] = [d0, length, width]
        Path(out_path).write_text(json.dumps(res, indent=1))
        return
    n_s = 6 if quick else 10
    sizes = np.concatenate([np.linspace(size_range[0], size_range[1], n_s // 2 + 1)[:-1], np.linspace(size_range[1], size_range[0], n_s // 2 + 1)[:-1]])
    size = [(float(d), length, width, 0.0) for d in sizes]
    res["ramps"][f"size_triangle_{size_range[0]:g}_{size_range[1]:g}"] = run_ramp("size", size, mesh_n)
    res["design_D_L_W_mm"] = [d0, length, width]
    Path(out_path).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 24
    design = tuple(float(v) for v in a[a.index("--design") + 1].split(",")) if "--design" in a else (5.6, 4.4, 1.2)
    sizes = tuple(float(v) for v in a[a.index("--sizes") + 1].split(",")) if "--sizes" in a else (5.2, 5.8)
    n_angle = int(a[a.index("--n-angle") + 1]) if "--n-angle" in a else None
    main(a[0], mesh, quick="--quick" in a, design=design, size_range=sizes, n_angle=n_angle, only_angle="--only-angle" in a)
