"""Phase 25 follow-up 5: the continuous, modest-modulation family -- square patches restricted to a
non-resonant size range, smooth sinusoidal phase modulation of moderate depth, same plane-wave locality check.

    python scripts/research/phase25_continuous_family.py <out.json>

6.4 mm cell, RT/duroid 5880 1.575 mm, PEC ground, 10.54 GHz, mesh 24 (0.267 mm; squares are even numbers of
mesh cells), 5x5 supercell (D = 32 mm = 1.12 lambda). Library: isolated patch-plane Gamma of every even-cell
square (same mesh and harmonic cut-off as the supercell solve). A block is built by mapping a target phase
phi_c + (depth/2) m(x, y), m in [-1, 1], to the nearest library square. Blocks:
  x_p1     m = cos(2 pi x / D)                  smooth, one period across the supercell, rows identical
  diag_p11 m = cos(2 pi (x + y) / D)            smooth, diagonal harmonic (period D/sqrt 2 = 22.6 mm)
  hex_crop m = cos(k_sw rho - k0 sin(th)(x cos ph + y sin ph)) over a 5x5 crop of hex150 port 1 (beam 30 deg,
            45 deg): the continuous (Fong-type) pattern of the real design's actual neighbourhood; NOT periodic,
            so the periodic supercell has one seam -- reported as such
The forward model is the isolated-cell Fourier sum, as in Part II of `phase25_forward_validation.py`.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase25_forward_validation as fv  # noqa: E402

MESH, N, PERIOD_MM = 24, 5, 6.4
F = fv.F0
PHI_C = 112.0  # centre of the modulation (deg), inside the non-resonant part of the curve
DEPTHS = (15.0, 30.0, 50.0)  # peak-to-peak phase swing, deg
MAX_NON_RESONANT_CELLS = 18  # 4.80 mm: slope <= ~30 deg/mm; 20 cells (5.33 mm) is the resonance onset (~70 deg/mm)


def library(sub, metal, p_m, f_matched):
    from src.phase12d.aperture_test import aperture_floquet_order

    del aperture_floquet_order
    h = PERIOD_MM / MESH
    lib = {0: (fv._empty_gamma(sub, F), 0.0)}
    for c in range(2, 21, 2):
        lib[c] = (fv._isolated_gamma(c * h, MESH, sub, metal, F, p_m, f_matched), c * h)
    return lib


def nearest_cells(phase_deg, lib, max_cells):
    keys = [c for c in lib if c <= max_cells]
    ph = np.array([np.degrees(np.angle(lib[c][0])) for c in keys])
    d = np.abs((ph[None, :] - np.asarray(phase_deg)[..., None] + 180.0) % 360.0 - 180.0)
    return np.array(keys)[np.argmin(d, axis=-1)]


def run_block(cells_map, lib, sub, metal, p_m):
    h = PERIOD_MM / MESH
    lw = [[(max(c * h, 1e-6), max(c * h, 1e-6)) for c in row] for row in cells_map]
    pw, n_unk = fv._mom_supercell(lw, MESH, sub, metal, F, p_m)
    gamma = np.vectorize(lambda c: lib[int(c)][0])(cells_map)
    xs = -N * p_m / 2 + (np.arange(N) + 0.5) * p_m
    xx, yy = np.meshgrid(xs, xs, indexing="xy")
    k0 = 2 * np.pi * F / 299792458.0
    rows = {}
    for m in range(-2, 3):
        for n in range(-2, 3):
            kx, ky = 2 * np.pi * m / (N * p_m), 2 * np.pi * n / (N * p_m)
            fwd = complex(np.mean(gamma * np.exp(1j * (kx * xx + ky * yy))))
            mom = complex(pw.e_reflected_x[2 + m, 2 + n])
            prop = kx * kx + ky * ky < k0 * k0
            w = (1.0 - (ky / k0) ** 2) / (np.sqrt(1.0 - (kx * kx + ky * ky) / k0**2)) if prop else 0.0
            rows[f"{m},{n}"] = {
                "propagating": bool(prop),
                "mom_abs": abs(mom),
                "forward_abs": abs(fwd),
                "mom_phase_deg": float(np.degrees(np.angle(mom))),
                "forward_phase_deg": float(np.degrees(np.angle(fwd))),
                "complex_error": abs(fwd - mom),
                "mom_power_fraction": float(pw.fraction[2 + m, 2 + n]),
                "forward_power_fraction": float(w * abs(fwd) ** 2),
            }
    prop_rows = [r for r in rows.values() if r["propagating"]]
    tot_m = sum(r["mom_power_fraction"] for r in prop_rows)
    tot_f = sum(r["forward_power_fraction"] for r in prop_rows)
    for r in prop_rows:
        r["power_error_pp"] = 100.0 * abs(r["forward_power_fraction"] / tot_f * tot_m - r["mom_power_fraction"])
    first = [rows[o] for o in ("1,0", "-1,0", "0,1", "0,-1", "1,1", "-1,-1") if o in rows]
    mod_orders = [
        k for k, r in rows.items() if k != "0,0" and r["mom_abs"] > 0.02 or r["forward_abs"] > 0.02 and k != "0,0"
    ]
    return {
        "n_unknowns": n_unk,
        "cells_map": np.asarray(cells_map).tolist(),
        "orders": rows,
        "worst_propagating_order_power_error_pp": max(r["power_error_pp"] for r in prop_rows),
        "forward_radiated_fraction": tot_f,
        "mom_radiated_fraction": tot_m,
        "modulated_orders_over_0p02": mod_orders,
        "max_modulated_order_amplitude_mom": max((rows[k]["mom_abs"] for k in rows if k != "0,0"), default=0.0),
        "max_modulated_order_amplitude_forward": max((rows[k]["forward_abs"] for k in rows if k != "0,0"), default=0.0),
        "specular_complex_error": rows["0,0"]["complex_error"],
        "first_orders": first,
    }


def main(out):
    from scipy.special import jv

    from src.dispersion.analytic_tm0 import tm0_beta
    from src.materials import Metal, Substrate
    from src.phase12d.aperture_test import aperture_floquet_order

    sub = Substrate("user", fv.EPS_R, fv.TAND, fv.H_MM * 1e-3)
    metal = Metal(model="pec")
    p_m = PERIOD_MM * 1e-3
    f_matched = aperture_floquet_order(N * MESH) // N
    t0 = time.time()
    lib = library(sub, metal, p_m, f_matched)
    curve = []
    for c, (g, side) in sorted(lib.items()):
        z = (1 + g) / (1 - g)
        curve.append(
            {
                "cells": c,
                "side_mm": side,
                "magnitude": abs(g),
                "phase_deg": float(np.degrees(np.angle(g))),
                "X_normalised": float(z.imag),
            }
        )
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"],
        "mesh_cells_per_6p4mm": MESH,
        "library_seconds": round(time.time() - t0, 1),
        "library": curve,
        "non_resonant_max_cells": MAX_NON_RESONANT_CELLS,
        "phi_c_deg": PHI_C,
        "cases": {},
    }
    Path(out).write_text(json.dumps(res, indent=1))
    xs = -N * p_m / 2 + (np.arange(N) + 0.5) * p_m
    xx, yy = np.meshgrid(xs, xs, indexing="xy")
    D = N * p_m
    # hex150 port-1 crop geometry
    hexres = fv._synthesize()
    port = hexres["ports"][0]
    gx, gy = np.asarray(port["x_m"]), np.asarray(port["y_m"])
    i0, j0 = 20, 20
    cx, cy = gx[i0 : i0 + N, j0 : j0 + N], gy[i0 : i0 + N, j0 : j0 + N]
    fx, fy = port["feed_x_mm"] * 1e-3, port["feed_y_mm"] * 1e-3
    k_sw = tm0_beta(F, fv.EPS_R, fv.H_MM * 1e-3)
    k0 = 2 * np.pi * F / 299792458.0
    th, ph = np.radians(30.0), np.radians(45.0)
    m_hex = np.cos(k_sw * np.hypot(cx - fx, cy - fy) - k0 * np.sin(th) * (cx * np.cos(ph) + cy * np.sin(ph)))
    patterns = {
        "x_p1": np.cos(2 * np.pi * xx / D),
        "diag_p11": np.cos(2 * np.pi * (xx + yy) / D),
        "hex_crop": m_hex,
    }
    # known answer: uniform block of the largest non-resonant square must be exact
    res["cases"]["uniform_4p8mm"] = run_block(np.full((N, N), MAX_NON_RESONANT_CELLS), lib, sub, metal, p_m)
    Path(out).write_text(json.dumps(res, indent=1))
    for depth in DEPTHS:
        max_cells = 20 if depth > 30.0 else MAX_NON_RESONANT_CELLS
        for name, m in patterns.items():
            cells = nearest_cells(PHI_C + 0.5 * depth * m, lib, max_cells)
            t1 = time.time()
            r = run_block(cells, lib, sub, metal, p_m)
            r["depth_deg_pp"] = depth
            r["max_cells_allowed"] = max_cells
            r["phase_realised_deg"] = np.vectorize(lambda c: float(np.degrees(np.angle(lib[int(c)][0]))))(
                cells
            ).tolist()
            # expected first-harmonic amplitude of an ideal sinusoidal phase modulation vs an ideal 180 deg 1-bit pair
            r["ideal_first_harmonic_J1"] = float(jv(1, np.radians(depth) / 2))
            r["ideal_1bit_180_first_harmonic"] = 2.0 / np.pi
            r["cost_vs_ideal_1bit_dB"] = float(20 * np.log10(r["ideal_first_harmonic_J1"] / (2.0 / np.pi)))
            res["cases"][f"depth{depth:g}_{name}"] = r
            Path(out).write_text(json.dumps(res, indent=1))
            print(
                f"[depth {depth:g} {name}] {time.time() - t1:.0f}s unk {r['n_unknowns']}: worst order {r['worst_propagating_order_power_error_pp']:.1f} pp, "
                f"max modulated |A| MoM {r['max_modulated_order_amplitude_mom']:.3f} fwd {r['max_modulated_order_amplitude_forward']:.3f}, "
                f"specular err {r['specular_complex_error']:.3f}",
                flush=True,
            )


if __name__ == "__main__":
    main(sys.argv[1])
