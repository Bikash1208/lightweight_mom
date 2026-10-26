"""Phase 26 B2.4e: SOLVER_SUPERCELL table -- beta and alpha of a slot-angle-rotating row vs modulation period
and depth (replaces the Floquet-of-local-tensors model, which fails the 0.01 k0 gate for the primary cell).

    python scripts/research/phase26_supercell_table.py <out.json> --slot L,W [--diam 6.0] [--cells 3,4,5,6,8] [--mesh 24]

One periodic MoM per (period, depth): N cells of the slotted circle, slot angle 180 i / N deg (i = 0..N-1)
along the wave direction x; the modulation vector is parallel to the propagation (gamma = 0), period N P,
K / k0 = lambda0 / (N P). Depth = slot size (anisotropy delta) at fixed diameter.

The complex root beta - j alpha is found by ANALYTIC CONTINUATION: `f(beta) = det(Y_s(beta) + diag(A, B))` is
analytic in beta (no branch point within 0.4 k0 of the root), is evaluated on 7 real points, interpolated
by a degree-6 polynomial in beta and its complex root nearest the real root taken; the same with a
different step and degree 5 give the stability. Residual check: |f| at the real part of the root against the
polynomial's own prediction. alpha includes dielectric loss: the uniform supercell of the slot-0 cell is
solved the same way as baseline (`alpha_uniform`) and `alpha_leak = alpha - alpha_uniform`.
Labels: SOLVER_SUPERCELL, MODEL_ONLY, NOT_HFSS_VALIDATED, gamma = 0 only.
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
LAMBDA0 = 299792458.0 / F


class Row:
    def __init__(self, geos, mesh_n):
        n = len(geos)
        self.n = n
        self.uc = UnitCell(RectangularLattice(n * P, P), CellRow(tuple(geos), P), SUB, Metal(), backing="pec")
        self.mesh = MeshSettings(n * mesh_n, mesh_n, max(n * mesh_n, mesh_n))
        self.basis = RooftopBasis(self.uc, self.mesh.cells_x, self.mesh.cells_y)
        self.orders = (n * mesh_n, mesh_n)
        self.omega = angular_frequency(F)

    def y(self, beta):
        return ts.bloch_sheet_tensor(self.uc, F, self.mesh, beta, 0.0, basis=self.basis, orders=self.orders)

    def f(self, beta):
        y = self.y(beta)
        a = _y_air(beta, self.omega, K0) + _y_slab(beta, self.omega, SUB.eps_complex, K0, SUB.thickness)
        b = ts._b_te(beta, self.omega, SUB.eps_complex, K0, SUB.thickness)
        return complex(np.linalg.det(y + np.array([[a, 0.0], [0.0, b]])))

    def f_imag(self, beta):
        y = self.y(beta)
        a = (_y_air(beta, self.omega, K0) + _y_slab(beta, self.omega, SUB.eps_complex, K0, SUB.thickness)).imag
        b = ts._b_te(beta, self.omega, SUB.eps_complex, K0, SUB.thickness).imag
        return float(np.linalg.det(y.imag + np.array([[a, 0.0], [0.0, b]])))

    def real_root(self, guess):
        b0, b1 = guess, guess * 1.004
        f0, f1 = self.f_imag(b0), self.f_imag(b1)
        for _ in range(8):
            if f1 == f0:
                break
            b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
            b0, f0, b1 = b1, f1, b2
            f1 = self.f_imag(b1)
            if abs(b1 - b0) < 1e-4 * K0:
                break
        return float(b1)

    def complex_root(self, beta_r, step_k0=0.004, m=3, deg=6):
        h = step_k0 * K0
        xs = beta_r + h * np.arange(-m, m + 1)
        vals = np.array([self.f(x) for x in xs])
        out = []
        for d, s in ((deg, 1.0), (deg - 1, 1.0)):
            poly = np.polynomial.Polynomial.fit((xs - beta_r) / h, vals, d)
            r = poly.roots()
            r = r[np.argmin(np.abs(r))]
            out.append(beta_r + h * complex(r))
        # residual: the true f at Re(root) against the polynomial's own value there
        poly = np.polynomial.Polynomial.fit((xs - beta_r) / h, vals, deg)
        br = out[0].real
        f_true = self.f(br)
        f_poly = poly((br - beta_r) / h)
        return out[0], abs(out[0] - out[1]), abs(f_true - f_poly) / max(abs(f_true), 1e-300)


def run_point(n, slot, diam, mesh_n):
    t0 = time.time()
    L, W = slot
    geos = [SlottedCircularPatch(diam * 1e-3, L * 1e-3, W * 1e-3, np.radians(180.0 * i / n)) for i in range(n)]
    row = Row(geos, mesh_n)
    cell0 = rl.single_cell(diam, L, W, 0.0)
    guess = rl.cell_beta(diam, L, W, 45.0, mesh_n)  # between the 1.04 and 1.17 extremes
    br = row.real_root(guess)
    bc, spread, resid = row.complex_root(br)
    # uniform baseline (loss only)
    uni = Row([SlottedCircularPatch(diam * 1e-3, L * 1e-3, W * 1e-3, 0.0)] * n, mesh_n)
    bu = uni.real_root(br * 0.97)
    bcu, _, _ = uni.complex_root(bu)
    k_mod = 2 * np.pi / (n * P)
    rec = {
        "n_cells": n, "period_mm": n * P * 1e3, "K_over_k0": k_mod / K0, "slot_L_W_mm": [L, W], "diameter_mm": diam,
        "beta_real_root_over_k0": br / K0, "beta_over_k0": bc.real / K0, "alpha_np_per_m": float(-bc.imag),
        "beta_uniform_over_k0": bcu.real / K0, "alpha_uniform_np_per_m": float(-bcu.imag),
        "alpha_leak_np_per_m": float(-bc.imag + bcu.imag),
        "root_spread_deg5_vs_deg6_over_k0": spread / K0, "poly_residual_rel": resid,
        "seconds": time.time() - t0,
    }
    print(rec, flush=True)
    return rec


def main(out, slot, diam, cells, mesh_n):
    res = {"labels": ["SOLVER_SUPERCELL", "MODEL_ONLY", "NOT_HFSS_VALIDATED", "GAMMA_0_COLLINEAR_ONLY"],
           "mesh_per_cell": mesh_n, "frequency_ghz": F / 1e9, "points": []}
    for n in cells:
        res["points"].append(run_point(n, slot, diam, mesh_n))
        Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    slot = tuple(float(v) for v in a[a.index("--slot") + 1].split(","))
    diam = float(a[a.index("--diam") + 1]) if "--diam" in a else 6.0
    cells = [int(v) for v in a[a.index("--cells") + 1].split(",")] if "--cells" in a else [3, 4, 5, 6, 8]
    mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 24
    main(a[0], slot, diam, cells, mesh)
