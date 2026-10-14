"""Phase 26 B: extraction validation and coverage map of the tensor-cell database.

    python scripts/research/phase26_tensor_coverage.py results/phase26

Reads `tensor_database.json`; writes `tensor_coverage.json` and `tensor_coverage.png`. Also re-solves three cells at
meshes 32/48/64/96 (staircase scatter of the rotated mask, the dominant numerical error of this cell family).
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase26_tensor_database as db  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor  # noqa: E402
from src.dispersion.sheet_model import find_beta_from_sheet_admittance  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402

R_MEAN_M = 0.14  # mean distance from the hexagon's centre to its edge (apothem 130, circumradius 150 mm)


def main(out_dir):
    out = Path(out_dir)
    data = json.loads((out / "tensor_database.json").read_text())
    rows = data["rows"]
    k0 = free_space_wavenumber(db.F)
    fam = {}
    for r in rows:
        fam.setdefault((r["L_mm"], r["W_mm"]), []).append(r)
    res = {"labels": data["labels"], "mesh": data["mesh"], "families": [], "validation": {}}

    # ---- extraction validation ------------------------------------------------------------
    recip = max(abs(complex(*r["y_xy"]) - complex(*r["y_xy"])) for r in rows)  # placeholder 0 (xy == yx stored once)
    asym = max(r["asym_S"] / max(abs(r["y_iso_im_S"]), 1e-12) for r in rows)
    axis_err, equiv = [], []
    for (length, width), rs in fam.items():
        rs = sorted(rs, key=lambda r: r["angle_deg"])
        base = next((r for r in rs if r["angle_deg"] == 0), None)
        if base is None or length == width:
            continue
        y0 = np.array([[complex(*base["y_xx"]), complex(*base["y_xy"])], [complex(*base["y_xy"]), complex(*base["y_yy"])]])
        for r in rs:
            a = r["angle_deg"]
            y = np.array([[complex(*r["y_xx"]), complex(*r["y_xy"])], [complex(*r["y_xy"]), complex(*r["y_yy"])]])
            d = (r["axis_deg"] - a + 90.0) % 180.0 - 90.0
            axis_err.append(abs(d))
            equiv.append(ts.rotation_equivariance_error(y, y0, np.radians(a)))
    res["validation"] = {
        "reciprocity_max_asym_over_yiso": asym,
        "rotation_axis_error_deg_max": float(max(axis_err)),
        "rotation_axis_error_deg_rms": float(np.sqrt(np.mean(np.square(axis_err)))),
        "rotation_equivariance_error_rel_max": float(max(equiv)),
        "rotation_equivariance_error_rel_median": float(np.median(equiv)),
        "note": "equivariance error includes the real square-lattice dependence of the neighbour gap on the rotation angle "
        "(the isotropic part grows towards 45 deg), not only staircase noise",
    }

    # ---- mesh scatter ---------------------------------------------------------------------
    scatter = {}
    for length, width, a in ((4.0, 2.4, 30.0), (4.0, 2.4, 45.0), (3.2, 1.6, 60.0)):
        row = []
        for n in (32, 48, 64, 96):
            y = db.solve_tensor(length, width, a, n)
            yi, delta, axis, _ = ts.decompose(y)
            row.append({"mesh": n, "y_iso_mS": float(yi.imag * 1e3), "delta_mS": delta * 1e3, "axis_deg": float(np.degrees(axis))})
        scatter[f"L{length}_W{width}_a{a:g}"] = row
    res["validation"]["mesh_scatter"] = scatter

    # ---- coverage -------------------------------------------------------------------------
    sub = db.SUB
    for (length, width), rs in sorted(fam.items()):
        rs = sorted(rs, key=lambda r: r["angle_deg"])
        yi = np.array([r["y_iso_im_S"] for r in rs])
        de = np.array([r["delta_S"] for r in rs])
        ybar, dbar = float(np.mean(yi)), float(np.mean(de))
        beta = find_beta_from_sheet_admittance(1j * ybar, sub, db.F).beta_rad_per_m
        theta = np.radians(30.0)
        alpha = leakage_alpha_first_order_tensor(beta, 1j * ybar, dbar, theta, sub, db.F)
        eta = 1.0 - np.exp(-2.0 * alpha * R_MEAN_M)
        full_rotation = bool(np.hypot(length, width) <= db.BBOX_MAX_MM + 1e-9)  # every axis angle fits the window
        need = {}
        for target in (0.5, 0.9):
            a_req = -np.log(1.0 - target) / (2.0 * R_MEAN_M)
            need[f"delta_for_{int(target * 100)}pct_mS"] = float(dbar * np.sqrt(a_req / alpha) * 1e3) if alpha > 0 else None
        res["families"].append(
            {
                "L_mm": length,
                "W_mm": width,
                "n_angles": len(rs),
                "y_iso_mean_mS": ybar * 1e3,
                "y_iso_spread_pct": float(100 * (yi.max() - yi.min()) / max(abs(ybar), 1e-12)),
                "delta_mean_mS": dbar * 1e3,
                "delta_spread_pct": float(100 * (de.max() - de.min()) / max(dbar, 1e-12)) if dbar > 1e-9 else None,
                "beta_over_k0_mean_sheet": beta / k0,
                "alpha_cp_leakage_np_per_m": alpha,
                "eta_conv_over_0p14m": float(eta),
                "full_rotation_in_window": full_rotation,
                **need,
            }
        )
    best = max(res["families"], key=lambda f: f["eta_conv_over_0p14m"])
    res["best_family_conversion"] = best
    full = [f for f in res["families"] if f["full_rotation_in_window"]]
    res["best_full_rotation_family"] = max(full, key=lambda f: f["delta_mean_mS"])
    # ---- beyond the non-resonant window: axis-aligned elongated patches (no rotation freedom claimed) ----
    beyond = []
    for length in (4.8, 5.2, 5.6, 6.0):
        for width in (0.8, 1.6, 2.4, 3.2):
            y = db.solve_tensor(length, width, 0.0, data["mesh"])
            yi, delta, _ax, _ = ts.decompose(y)
            try:
                beta = find_beta_from_sheet_admittance(1j * yi.imag, sub, db.F).beta_rad_per_m
                alpha = leakage_alpha_first_order_tensor(beta, 1j * yi.imag, delta, np.radians(30.0), sub, db.F)
            except Exception:
                beta, alpha = float("nan"), float("nan")
            beyond.append(
                {
                    "L_mm": length,
                    "W_mm": width,
                    "gap_to_neighbour_mm": 6.4 - length,
                    "y_xx_mS": float(y[0, 0].imag * 1e3),
                    "y_yy_mS": float(y[1, 1].imag * 1e3),
                    "delta_mS": delta * 1e3,
                    "alpha_np_per_m": alpha,
                    "eta_conv_over_0p14m": float(1.0 - np.exp(-2.0 * alpha * R_MEAN_M)),
                }
            )
    res["beyond_nonresonant_axis_aligned"] = beyond
    (out / "tensor_coverage.json").write_text(json.dumps(res, indent=1))

    # ---- picture --------------------------------------------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for r in rows:
        ax[0].scatter(r["y_iso_im_S"] * 1e3, r["delta_S"] * 1e3, c=r["L_mm"], cmap="viridis", vmin=1.6, vmax=4.8, s=14)
    ax[0].scatter([0.0], [0.0], c="k", marker="s")
    ax[0].annotate("empty", (0.0, 0.0), textcoords="offset points", xytext=(5, 5))
    ax[0].set_xlabel("isotropic part  Im(Y_iso)  [mS]")
    ax[0].set_ylabel("anisotropy  delta  [mS]")
    ax[0].set_title("reachable (Y_iso, delta); every point rotatable to any axis")
    for f in res["families"]:
        ax[1].scatter(f["delta_mean_mS"], 100 * f["eta_conv_over_0p14m"], c="C1", s=18)
    ax[1].set_xlabel("anisotropy delta [mS]")
    ax[1].set_ylabel("CP-leakage conversion over 0.14 m [%]")
    ax[1].set_title("first-order conversion reachable (non-resonant range)")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "tensor_coverage.png", dpi=130)
    print(json.dumps(res["validation"], indent=1))
    print("best conversion family:", best)
    print("best full-rotation family:", res["best_full_rotation_family"])
    for b in beyond:
        print("beyond", b["L_mm"], b["W_mm"], "delta %.2f mS alpha %.3f eta %.3f" % (b["delta_mS"], b["alpha_np_per_m"], b["eta_conv_over_0p14m"]))
    y = sorted(res["families"], key=lambda f: -f["delta_mean_mS"])[:5]
    for f in y:
        print(f["L_mm"], f["W_mm"], "delta %.3f mS  yiso %.3f mS  alpha %.3f  eta %.4f  beta/k0 %.4f spread y %.1f%% d %s" % (
            f["delta_mean_mS"], f["y_iso_mean_mS"], f["alpha_cp_leakage_np_per_m"], f["eta_conv_over_0p14m"], f["beta_over_k0_mean_sheet"], f["y_iso_spread_pct"], f["delta_spread_pct"]))


if __name__ == "__main__":
    main(sys.argv[1])
