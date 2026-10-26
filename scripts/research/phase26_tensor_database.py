"""Phase 26 B: tensor-cell database (L, W, rotation) -> 2x2 sheet admittance, and its coverage map.

    python scripts/research/phase26_tensor_database.py <out_dir> [--mesh 48] [--quick]

Rotated rectangular PEC patches on the 6.4 mm cell, grounded RT/duroid 5880 1.575 mm, 10.54 GHz, normal-
incidence MoM -> `tensor_sheet.sheet_admittance_tensor_from_gamma`. Only the NON-RESONANT range: the
rotated rectangle's bounding box stays <= 4.8 mm in x and y (gap to the neighbouring cell's box >= 1.6 mm,
the Phase 25 locality criterion). Rotations 0..90 deg are solved; 90..180 deg follow from the lattice mirror
symmetry (`Y(-a) = M Y(a) M`, `M = diag(1, -1)`, exact for a centred pattern on this mesh).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402
from src.geometry import Geometry, UnitCell  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

F, P = 10.54e9, 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
BBOX_MAX_MM = 4.8
LENGTHS_MM = [1.6, 2.4, 3.2, 4.0, 4.8]
WIDTHS_MM = [0.8, 1.6, 2.4, 3.2, 4.0, 4.8]
ANGLES = [0, 15, 30, 45, 60, 75, 90]


def bbox_mm(length, width, angle_deg):
    c, s = abs(np.cos(np.radians(angle_deg))), abs(np.sin(np.radians(angle_deg)))
    return length * c + width * s, length * s + width * c


def solve_tensor(length_mm, width_mm, angle_deg, mesh_n):
    geo = Geometry("rectangle", length_mm * 1e-3, width_mm * 1e-3, rotation=np.radians(angle_deg))
    cell = UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")
    return ts.sheet_admittance_tensor_for_cell(cell, F, MeshSettings(mesh_n, mesh_n, mesh_n))


def mirror(y):
    m = np.diag([1.0, -1.0])
    return m @ y @ m


def entry(length, width, angle, y, with_beta=True):
    yi, delta, axis, asym = ts.decompose(y)
    out = {
        "L_mm": length,
        "W_mm": width,
        "angle_deg": angle,
        "y_iso_im_S": float(yi.imag),
        "y_iso_re_S": float(yi.real),
        "delta_S": delta,
        "axis_deg": float(np.degrees(axis)),
        "asym_S": float(asym),
        "y_xx": [float(y[0, 0].real), float(y[0, 0].imag)],
        "y_yy": [float(y[1, 1].real), float(y[1, 1].imag)],
        "y_xy": [float(y[0, 1].real), float(y[0, 1].imag)],
    }
    if with_beta:
        k0 = free_space_wavenumber(F)
        try:
            b_par, _ = ts.find_beta_tensor(y, axis, SUB, F)
            b_perp, _ = ts.find_beta_tensor(y, axis + np.pi / 2, SUB, F)
            b_45, c45 = ts.find_beta_tensor(y, axis + np.pi / 4, SUB, F)
            out.update(beta_par_over_k0=b_par / k0, beta_perp_over_k0=b_perp / k0, beta_45_over_k0=b_45 / k0, hybrid_45=c45)
        except Exception as exc:  # no bound mode: recorded, never invented
            out["beta_error"] = type(exc).__name__
    return out


def main(out_dir, mesh_n=48, quick=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows = []
    lengths = LENGTHS_MM[1::2] if quick else LENGTHS_MM
    widths = WIDTHS_MM[1::2] if quick else WIDTHS_MM
    for length in lengths:
        for width in widths:
            if width > length:
                continue
            solved = {}
            for a in ANGLES:
                bx, by = bbox_mm(length, width, a)
                if max(bx, by) > BBOX_MAX_MM + 1e-9:
                    continue
                if width == length and a not in (0, 45):
                    continue  # a square is isotropic at every angle but 45 deg changes the gap
                solved[a] = solve_tensor(length, width, a, mesh_n)
            for a, y in solved.items():
                rows.append(entry(length, width, a, y))
                if 0 < a < 90 and length != width:
                    rows.append(entry(length, width, 180 - a, mirror(y)))
            print(f"L={length} W={width}: angles {sorted(solved)}  [{time.time() - t0:.0f}s]", flush=True)
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "TENSOR_SHEET_FROM_NORMAL_INCIDENCE"],
        "frequency_ghz": F / 1e9,
        "period_mm": P * 1e3,
        "mesh": mesh_n,
        "bbox_max_mm": BBOX_MAX_MM,
        "rows": rows,
        "seconds": time.time() - t0,
    }
    (out_dir / "tensor_database.json").write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    a = sys.argv[1:]
    mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 48
    main(a[0], mesh, quick="--quick" in a)
