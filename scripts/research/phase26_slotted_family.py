"""Phase 26 B2.2/B2.3: the slotted-patch cell family at the operating point (Bloch kt > k0), its tensor
(mean loading, anisotropy), and the error of the normal-incidence sheet model for these cells.

    python scripts/research/phase26_slotted_family.py <out_dir> [--mesh 48] [--quick]

For every cell (circular patch D with a centred rectangular slot L x W, slot axis along x; bare circles; slotted
squares as the comparison): `tensor_sheet.bloch_beta` along phi = 0 (along the slot), 45, 90 deg; the Bloch
sheet tensor at the fixed operating wavenumber kt = 1.05 k0 along and across the slot (-> mean loading
Y_bar = (Y_par + Y_perp)/2 and anisotropy delta = |Y_par - Y_perp|/2); the normal-incidence tensor model's
beta for the same directions (its error). A staircase-rotation check on a few cells.
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
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor  # noqa: E402
from src.dispersion.sheet_model import find_beta_from_sheet_admittance  # noqa: E402
from src.dispersion.slotted_cells import SlottedCircularPatch, SlottedSquarePatch  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

F, P = 10.54e9, 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
K0 = free_space_wavenumber(F)
KT_OP = 1.05 * K0
R_MEAN = 0.14


def make_cell(kind, size, length, width, angle_deg=0.0):
    ang = np.radians(angle_deg)
    if kind == "circle":
        geo = SlottedCircularPatch(size * 1e-3, length * 1e-3, width * 1e-3, ang)
    else:
        geo = SlottedSquarePatch(size * 1e-3, length * 1e-3, width * 1e-3, ang)
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")


def cell_list(quick=False):
    cells = []
    for d in (5.4, 5.8, 6.0):
        cells.append(("circle", d, 0.0, 0.0))
    slots = [(4.6, 1.2), (5.0, 1.6), (5.0, 2.0), (5.4, 1.2), (5.4, 1.6)]
    for d in (5.4, 5.6, 5.8, 6.0):
        for length, width in slots:
            chord = 2.0 * np.sqrt((d / 2) ** 2 - (width / 2) ** 2)
            if length <= chord - 0.25:
                cells.append(("circle", d, length, width))
    for s in (5.4, 5.8):
        for length, width in ((4.6, 1.2), (5.0, 1.6)):
            cells.append(("square", s, length, width))
    return cells[::3] if quick else cells


def evaluate(kind, size, length, width, mesh_n):
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    cell = make_cell(kind, size, length, width)
    basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    out = {"kind": kind, "size_mm": size, "slot_L_mm": length, "slot_W_mm": width, "mesh": mesh_n}
    y_norm = ts.sheet_admittance_tensor_for_cell(cell, F, mesh)
    out["normal_incidence_y_xx_yy_mS"] = [y_norm[0, 0].imag * 1e3, y_norm[1, 1].imag * 1e3]
    bloch, normal = {}, {}
    for phi in (0.0, 45.0, 90.0):
        try:
            b, _ = ts.bloch_beta(cell, F, mesh, np.radians(phi), n_grid=5, refine=1)
            bloch[phi] = b / K0
        except Exception:
            bloch[phi] = None
        try:
            normal[phi] = ts.find_beta_tensor(y_norm, np.radians(phi), SUB, F)[0] / K0
        except Exception:
            normal[phi] = None
    out["beta_bloch_over_k0"] = {str(int(k)): v for k, v in bloch.items()}
    out["beta_normal_incidence_over_k0"] = {str(int(k)): v for k, v in normal.items()}
    out["beta_error_normal_incidence_over_k0"] = {
        str(int(k)): (None if bloch[k] is None or normal[k] is None else normal[k] - bloch[k]) for k in bloch
    }
    y_par = ts.bloch_sheet_tensor(cell, F, mesh, KT_OP, 0.0, basis=basis)
    y_perp = ts.bloch_sheet_tensor(cell, F, mesh, 0.0, KT_OP, basis=basis)
    yp, yq = y_par[0, 0].imag, y_perp[1, 1].imag
    out["y_par_mS"], out["y_perp_mS"] = yp * 1e3, yq * 1e3
    out["y_bar_mS"], out["delta_mS"] = 0.5 * (yp + yq) * 1e3, 0.5 * abs(yp - yq) * 1e3
    out["y_offdiag_par_mS"] = float(abs(y_par[0, 1]) * 1e3)
    # conversion estimate at this operating point (first order; the depth caveat of the required map applies)
    beta = find_beta_from_sheet_admittance(1j * 0.5 * (yp + yq), SUB, F).beta_rad_per_m
    alpha = leakage_alpha_first_order_tensor(beta, 1j * 0.5 * (yp + yq), 0.5 * abs(yp - yq), np.radians(30.0), SUB, F)
    out["beta_mean_sheet_over_k0"] = beta / K0
    out["alpha_np_per_m"] = alpha
    out["eta_conv_hexagon_est"] = float(1.0 - np.exp(-2.0 * alpha * R_MEAN))
    return out


def rotation_check(kind, size, length, width, mesh_n, angles=(30.0, 45.0, 60.0, 90.0)):
    """Joint rotation of cell and wave: Y(theta-slot, k along x) = R(theta) Y(0-slot, k at -theta) R(theta)^T."""
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    c0 = make_cell(kind, size, length, width, 0.0)
    b0 = RooftopBasis(c0, mesh.cells_x, mesh.cells_y)
    rows = []
    for a in angles:
        th = np.radians(a)
        ca = make_cell(kind, size, length, width, a)
        ya = ts.bloch_sheet_tensor(ca, F, mesh, KT_OP, 0.0)
        y0 = ts.bloch_sheet_tensor(c0, F, mesh, KT_OP * np.cos(-th), KT_OP * np.sin(-th), basis=b0)
        pred = ts.rotate_tensor(y0, th)
        rows.append({"angle_deg": a, "relative_error": float(np.linalg.norm(ya - pred) / np.linalg.norm(ya))})
    return rows


def main(out_dir, mesh_n=48, quick=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows = []
    for kind, size, length, width in cell_list(quick):
        r = evaluate(kind, size, length, width, mesh_n)
        rows.append(r)
        print(
            f"{kind} {size} slot {length}x{width}: ybar {r['y_bar_mS']:.2f} delta {r['delta_mS']:.2f} "
            f"beta_bloch {r['beta_bloch_over_k0']} err_normal {r['beta_error_normal_incidence_over_k0']} "
            f"[{time.time() - t0:.0f}s]",
            flush=True,
        )
        (out / "slotted_family.json").write_text(json.dumps({"mesh": mesh_n, "rows": rows}, indent=1))
    rot = {}
    for spec in (("circle", 5.8, 5.0, 2.0), ("circle", 5.6, 4.6, 1.2)):
        rot[f"{spec[0]}_{spec[1]}_{spec[2]}x{spec[3]}"] = rotation_check(*spec, mesh_n)
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "BLOCH_SHEET_TENSOR_AT_OPERATING_KT"], "mesh": mesh_n, "kt_operating_over_k0": 1.05, "rows": rows, "rotation_check": rot, "seconds": time.time() - t0}
    (out / "slotted_family.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv[1:]
    mesh = int(a[a.index("--mesh") + 1]) if "--mesh" in a else 48
    main(a[0], mesh, quick="--quick" in a)
