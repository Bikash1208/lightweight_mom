"""Phase 26 B2.4b: internal cross-check -- database Bloch beta vs the UNIFORM supercell solve.

    python scripts/research/phase26_database_crosscheck.py <out.json> [--mesh 24] [--n-cells 4] [--design D,L,W]

For the item-C cells, the single-cell Bloch beta (`tensor_sheet.bloch_beta`: the database route) against the
same cell repeated `N` times in one periodic MoM along x (supercell orders `N m x m`), propagation along x
(phi = 0) and the diagonal (phi = 45 deg; Bloch wavevector beta (cos, sin) on the same supercell). Gate
0.005 k0. A uniform supercell must reproduce the cell exactly, so this tests the Floquet-order truncation,
the staircase masks and the root selection -- not physics. Labels: SOLVER_INTERNAL, NOT_HFSS_VALIDATED.
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
from src.frequency import angular_frequency  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

F, P, SUB, K0 = rl.F, rl.P, rl.SUB, rl.K0
GATE_K0 = 0.005


def uniform_supercell_beta(geo, n, mesh_n, phi, guess):
    uc = UnitCell(RectangularLattice(n * P, P), CellRow(tuple([geo] * n), P), SUB, Metal(), backing="pec")
    mesh = MeshSettings(n * mesh_n, mesh_n, max(n * mesh_n, mesh_n))
    basis = RooftopBasis(uc, mesh.cells_x, mesh.cells_y)
    omega = angular_frequency(F)
    p = np.array([np.cos(phi), np.sin(phi)])
    q = np.array([-np.sin(phi), np.cos(phi)])

    def res(beta):
        y = ts.bloch_sheet_tensor(uc, F, mesh, beta * p[0], beta * p[1], basis=basis, orders=(n * mesh_n, mesh_n))
        a = (_y_air(beta, omega, K0) + _y_slab(beta, omega, SUB.eps_complex, K0, SUB.thickness)).imag
        b = ts._b_te(beta, omega, SUB.eps_complex, K0, SUB.thickness).imag
        return float(np.linalg.det(y.imag + a * np.outer(p, p) + b * np.outer(q, q)))

    b0, b1 = guess, guess * 1.004
    f0, f1 = res(b0), res(b1)
    for _ in range(6):
        if f1 == f0:
            break
        b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
        b0, f0, b1 = b1, f1, b2
        f1 = res(b1)
        if abs(b1 - b0) < 1e-4 * K0:
            break
    return float(b1)


def main(out_path, mesh_n, n, design):
    d, l, w = design
    res = {"labels": ["SOLVER_INTERNAL", "NOT_HFSS_VALIDATED"], "mesh_per_cell": mesh_n, "n_cells": n, "gate_k0": GATE_K0, "cases": []}
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    for ang in (0.0, 45.0, 90.0):
        geo = SlottedCircularPatch(d * 1e-3, l * 1e-3, w * 1e-3, np.radians(ang))
        cell = rl.single_cell(d, l, w, ang)
        for phi_deg in (0.0, 45.0):
            t0 = time.time()
            phi = np.radians(phi_deg)
            b_db, _ = ts.bloch_beta(cell, F, mesh, phi, n_grid=5, refine=1)
            b_sc = uniform_supercell_beta(geo, n, mesh_n, phi, b_db)
            err = (b_sc - b_db) / K0
            rec = {"slot_angle_deg": ang, "propagation_deg": phi_deg, "database_beta_over_k0": b_db / K0,
                   "supercell_beta_over_k0": b_sc / K0, "error_k0": err, "gate": "PASS" if abs(err) <= GATE_K0 else "FAIL",
                   "seconds": time.time() - t0}
            res["cases"].append(rec)
            print(rec, flush=True)
            Path(out_path).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 24
    n = int(a[a.index("--n-cells") + 1]) if "--n-cells" in a else 4
    design = tuple(float(v) for v in a[a.index("--design") + 1].split(",")) if "--design" in a else (6.0, 4.8, 2.0)
    main(a[0], mesh, n, design)
