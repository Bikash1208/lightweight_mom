"""Phase 25 follow-up 1: is the large local-periodicity error of the 0.55 mm-gap 1-bit pair a mesh artefact?

    python scripts/research/phase25_block_locality.py <out.json> [mesh ...]

Small square|empty supercells (1x2 stripes, 2x2 checkerboard, 2x2 with one empty cell) of the design pair
(5.867 mm square = 22 cells at mesh 24 = 44 cells at mesh 48, gap 0.533 mm, grounded RT/duroid 5880 1.575 mm,
10.54 GHz), plane wave at normal incidence, solved as full periodic supercells. A supercell of period
12.8 mm < lambda radiates everything into the specular order, so the comparison is the complex E_x amplitude
of every Floquet order |m|,|n| <= 2 at the patch plane (the evanescent ones are the near-field the cells
couple through): forward model (Fourier sum of the ISOLATED cell phasors solved at the same mesh and
harmonic cut-off) against the MoM. The 1x2 stripe is a 2x2 supercell with identical rows (the solver's
supercell is square; the field is periodic with 6.4 mm in y either way).
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

F0, P_MM = fv.F0, fv.PERIOD
CASES = {
    "1x2_stripes": np.array([[1, 0], [1, 0]]),
    "2x2_checkerboard": np.array([[1, 0], [0, 1]]),
    "2x2_three_square_one_empty": np.array([[1, 1], [1, 0]]),
}
SIDE_CELLS = {24: 22, 48: 44}


def run(mesh, side_cells=None, cases=None, n=2, frequency_hz=F0, tag="", period_mm=P_MM):
    from src.materials import Metal, Substrate
    from src.phase12d.aperture_test import aperture_floquet_order

    substrate = Substrate("user", fv.EPS_R, fv.TAND, fv.H_MM * 1e-3)
    metal = Metal(model="pec")
    p_m = period_mm * 1e-3
    h = period_mm / mesh
    cells = side_cells or SIDE_CELLS[mesh]
    side = cells * h
    f_matched = aperture_floquet_order(n * mesh) // n
    g_empty = fv._empty_gamma(substrate, frequency_hz)
    g_sq = fv._isolated_gamma(side, mesh, substrate, metal, frequency_hz, p_m, f_matched)
    out = {
        "mesh_cells_per_cell": mesh,
        "square_side_mm": side,
        "gap_mm": period_mm - side,
        "period_mm": period_mm,
        "frequency_GHz": frequency_hz / 1e9,
        "isolated_empty": [abs(g_empty), float(np.degrees(np.angle(g_empty)))],
        "isolated_square": [abs(g_sq), float(np.degrees(np.angle(g_sq)))],
        "cases": {},
    }
    xs = -n * p_m / 2 + (np.arange(n) + 0.5) * p_m
    xx, yy = np.meshgrid(xs, xs, indexing="xy")
    for name, lev in (cases or CASES).items():
        t0 = time.time()
        lw = [[(side, side) if lev[i, j] else (1e-6, 1e-6) for j in range(n)] for i in range(n)]
        pw, n_unk = fv._mom_supercell(lw, mesh, substrate, metal, frequency_hz, p_m)
        gamma = np.where(lev == 1, g_sq, g_empty)
        rows = {}
        for m in range(-2, 3):
            for k in range(-2, 3):
                kx, ky = 2 * np.pi * m / (n * p_m), 2 * np.pi * k / (n * p_m)
                fwd = complex(np.mean(gamma * np.exp(1j * (kx * xx + ky * yy))))
                mom = complex(pw.e_reflected_x[2 + m, 2 + k])
                rows[f"{m},{k}"] = {
                    "mom": [abs(mom), float(np.degrees(np.angle(mom)))],
                    "forward": [abs(fwd), float(np.degrees(np.angle(fwd)))],
                    "complex_error": abs(fwd - mom),
                }
        out["cases"][name] = {
            "n_unknowns": n_unk,
            "seconds": round(time.time() - t0, 1),
            "orders": rows,
            "specular_complex_error": rows["0,0"]["complex_error"],
            "specular_amplitude_mom_forward": [rows["0,0"]["mom"][0], rows["0,0"]["forward"][0]],
            "worst_order_complex_error": max(r["complex_error"] for r in rows.values()),
            "mom_radiated_fraction": float(pw.fraction[pw.propagating].sum()),
        }
        c = out["cases"][name]
        print(
            f"[{tag}mesh{mesh} {name}] {c['seconds']}s unk {n_unk}: specular |A| MoM {c['specular_amplitude_mom_forward'][0]:.3f}"
            f" fwd {c['specular_amplitude_mom_forward'][1]:.3f}, specular err {c['specular_complex_error']:.3f},"
            f" worst order err {c['worst_order_complex_error']:.3f}",
            flush=True,
        )
    return out


if __name__ == "__main__":
    result = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "meshes": {}}
    for mesh in [int(a) for a in sys.argv[2:]] or [24, 48]:
        result["meshes"][f"mesh{mesh}"] = run(mesh)
        Path(sys.argv[1]).write_text(json.dumps(result, indent=1))
