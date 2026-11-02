"""Phase 30 item 2: coarse cell screening of candidate geometries x substrates at 10.54 GHz (MODEL_ONLY).

    python scripts/research/phase30_cell_screen.py <out_dir> [--mesh 48] [--workers 44] [--quick]

Per cell (axes of the cell along x/y): normal-incidence tensor -> mean-sheet beta0; Bloch sheet tensor at kt = beta0
along x (Y_par, E along the cell's x axis) and along y (Y_perp) -> mean loading Y_bar, anisotropy delta; the
surface-wave reactances X = Im 1/(Y + Y_slab(beta_bar)), tensor modulation M_t = |X_par - X_perp|/|X_par + X_perp|
and anisotropy ratio X_hi/X_lo; first-order leakage alpha of a rotating-axis tensor of THIS cell's anisotropy
(`leakage_alpha_first_order_tensor`, theta = 17 deg); dielectric loss from the complex transverse-resonance root.
Edge power on a hexagon of effective radius R = 70 mm: uniform modulation exp(-2 alpha R); ideal uniform-illumination
taper bound 1/(1 + alpha R).  Conductor loss: NOT estimated (PEC model).
"""

import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.elliptical_cells import EllipticalPatch  # noqa: E402
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor  # noqa: E402
from src.dispersion.sheet_model import _y_air, _y_slab, find_beta_from_sheet_admittance  # noqa: E402
from src.frequency import angular_frequency, free_space_wavenumber  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

F = 10.54e9
K0 = free_space_wavenumber(F)
R_EDGE = 0.070
THETA = np.radians(17.0)
GAPS = (1.5e-3, 0.9e-3, 0.5e-3, 0.3e-3, 0.2e-3, 0.15e-3, 0.1e-3)  # edge-to-edge gap to the neighbour along the major axis [m]
# name: (eps_r, tan_delta (datasheet nominal), thickness m, pitch m)
SUBSTRATES = {
    "RT5880": (2.2, 0.0009, 1.575e-3, 6.4e-3),
    "RO3006": (6.15, 0.0020, 0.635e-3, 4.0e-3),
    "RO3010": (10.2, 0.0022, 0.635e-3, 3.2e-3),
}
RATIOS = (1.0, 0.92, 0.84, 0.75, 0.62, 0.48)  # minor / major


PITCHES = {"RT5880": (6.4e-3, 4.8e-3), "RO3006": (4.0e-3,), "RO3010": (3.2e-3,)}


def cell_specs(quick=False):
    """(substrate, family, major, minor/major, slot_L, slot_W, cross, exponent, pitch)."""
    specs = []
    for sub, pitches in PITCHES.items():
        for pitch in pitches:
            for gap in GAPS:
                major = pitch - gap
                for n, tag in ((2.0, "ellipse"), (4.0, "squircle4"), (8.0, "rect8")):
                    for r in RATIOS:
                        specs.append((sub, tag, major, r, 0.0, 0.0, False, n, pitch))
                for r in (1.0, 0.75):
                    for wf in (0.10, 0.18):  # slot width / pitch; slot length 0.62 major along the major axis
                        specs.append((sub, "squircle4+slot", major, r, 0.62 * major, wf * pitch, False, 4.0, pitch))
                        specs.append((sub, "ellipse+slot", major, r, 0.62 * major, wf * pitch, False, 2.0, pitch))
                specs.append((sub, "squircle4+cross", major, 1.0, 0.70 * major, 0.12 * pitch, True, 4.0, pitch))
                specs.append((sub, "circle+cross", major, 1.0, 0.70 * major, 0.12 * pitch, True, 2.0, pitch))
                specs.append((sub, "circle+slot", major, 1.0, 0.70 * major, 0.14 * pitch, False, 2.0, pitch))
    return specs[::7] if quick else specs


def _dielectric_alpha(y_sheet, substrate, beta):
    """|Im beta| [Np/m] of the complex root of the scalar transverse resonance with the lossy slab (secant)."""
    omega = angular_frequency(F)
    eps_c = substrate.eps_complex

    def g(b):
        return y_sheet + _y_air(b, omega, K0) + _y_slab(b, omega, eps_c, K0, substrate.thickness)

    b0, b1 = complex(beta), complex(beta) * (1 + 1e-6)
    g0, g1 = g(b0), g(b1)
    for _ in range(30):
        if g1 == g0:
            break
        b0, b1, g0 = b1, b1 - g1 * (b1 - b0) / (g1 - g0), g1
        g1 = g(b1)
        if abs(g1) < 1e-14:
            break
    return abs(b1.imag), abs(b1.real - beta) / beta


def evaluate(args):
    spec, mesh_n = args
    sub_name, kind, major, ratio, slot_l, slot_w, cross, expo, pitch = spec
    eps, tand, h, _ = SUBSTRATES[sub_name]
    sub = Substrate(sub_name.lower(), eps, tand, h)
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    out = {"substrate": sub_name, "family": kind, "pitch_mm": pitch * 1e3, "exponent": expo, "major_mm": major * 1e3, "ratio": ratio,
           "slot_L_mm": slot_l * 1e3, "slot_W_mm": slot_w * 1e3, "gap_mm": (pitch - major) * 1e3,
           "gap_mesh_cells": (pitch - major) / (pitch / mesh_n)}  # gap < 1.5 mesh cells is NOT resolved: rows invalid
    try:
        geo = EllipticalPatch(major, major * ratio, 0.0, slot_l, slot_w, 0.0, cross, expo)
        cell = UnitCell(RectangularLattice(pitch, pitch), geo, sub, Metal(), backing="pec")
        basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
        yn = ts.sheet_admittance_tensor_for_cell(cell, F, mesh)
        y0 = 1j * 0.5 * (yn[0, 0].imag + yn[1, 1].imag)
        beta0 = find_beta_from_sheet_admittance(y0, sub, F).beta_rad_per_m
        yp = ts.bloch_sheet_tensor(cell, F, mesh, beta0, 0.0, basis=basis)
        yq = ts.bloch_sheet_tensor(cell, F, mesh, 0.0, beta0, basis=basis)
        y_par, y_perp = yp[0, 0].imag, yq[1, 1].imag
        y_bar, delta = 0.5 * (y_par + y_perp), 0.5 * abs(y_par - y_perp)
        mode = find_beta_from_sheet_admittance(1j * y_bar, sub, F)
        beta = mode.beta_rad_per_m
        omega = angular_frequency(F)
        y_slab = _y_slab(beta, omega, sub.eps_complex, K0, sub.thickness)
        x_par = float(np.imag(1.0 / (1j * y_par + y_slab)))
        x_perp = float(np.imag(1.0 / (1j * y_perp + y_slab)))
        alpha = leakage_alpha_first_order_tensor(beta, 1j * y_bar, delta, THETA, sub, F)
        a_d, _ = _dielectric_alpha(1j * y_bar, sub, beta)
        out.update({
            "ok": True, "y_par_mS": y_par * 1e3, "y_perp_mS": y_perp * 1e3, "y_bar_mS": y_bar * 1e3, "delta_mS": delta * 1e3,
            "beta_over_k0": beta / K0, "x_bar_ohm": mode.x_s_ohm, "x_par_ohm": x_par, "x_perp_ohm": x_perp,
            "M_t": abs(x_par - x_perp) / abs(x_par + x_perp),
            "aniso_ratio": max(abs(x_par), abs(x_perp)) / max(min(abs(x_par), abs(x_perp)), 1e-9),
            "alpha_np_m": alpha, "alpha_over_k0": alpha / K0,
            "alpha_diel_np_m": a_d, "rad_over_diel": alpha / a_d if a_d > 0 else None,
            "edge_frac_uniform": float(np.exp(-2.0 * alpha * R_EDGE)),
            "edge_frac_ideal_taper_bound": float(1.0 / (1.0 + alpha * R_EDGE)),
            "offdiag_par_mS": float(abs(yp[0, 1]) * 1e3),
        })
    except Exception as e:  # noqa: BLE001 -- record, never fabricate
        out.update({"ok": False, "error": f"{type(e).__name__}: {e}"})
    return out


def main(out_dir, mesh_n=48, workers=44, quick=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    specs = [sp for sp in cell_specs(quick) if (sp[8] - sp[2]) >= 1.5 * sp[8] / mesh_n - 1e-12]  # resolved gaps only
    t0 = time.time()
    with Pool(workers) as pool:
        rows = []
        for i, r in enumerate(pool.imap(evaluate, [(s, mesh_n) for s in specs], chunksize=1)):
            rows.append(r)
            with open(out / 'cell_screen.jsonl', 'a') as fh:
                fh.write(json.dumps(r) + "\n")
            if i % 20 == 0:
                print(f"{i + 1}/{len(specs)} [{time.time() - t0:.0f}s]", flush=True)
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "BLOCH_SHEET_TENSOR_AT_BETA0"], "f_hz": F, "mesh": mesh_n,
           "R_edge_m": R_EDGE, "theta_deg": 17.0, "substrates": SUBSTRATES, "n": len(rows), "seconds": time.time() - t0, "rows": rows}
    (out / "cell_screen.json").write_text(json.dumps(res, indent=1))
    print("done", len(rows), "ok", sum(r.get("ok", False) for r in rows), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--mesh", type=int, default=48)
    ap.add_argument("--workers", type=int, default=44)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    main(a.out_dir, a.mesh, a.workers, a.quick)
