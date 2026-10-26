"""Phase 26 B2 step 1: the (mean loading, anisotropy) region that converts >= 85 % of the launched surface wave
inside the 150 mm hexagon -- the TARGET for the tensor cell family.

    python scripts/research/phase26_required_map.py results/phase26

Mean sheet `Y_bar = j y mS` (sets beta/k0 through `sheet_model`), rotating-axis tensor anisotropy `delta`
(`leakage_model.leakage_alpha_first_order_tensor`, beam at theta 30 deg), conversion averaged over the hexagon's
azimuths with `leaky_forward.polygon_ray_distance`. Launch fraction from `pin_launch` (loaded slab, the current
pin, both PinLength readings). Directivity of the illuminated hexagon aperture from the continuous-phase scalar
aperture sum with the `rho^-1/2 exp(-alpha rho)` illumination (an estimate; item C replaces it).
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from src.dispersion.leakage_model import leakage_alpha_first_order_tensor  # noqa: E402
from src.dispersion.pin_launch import Launcher, PinLaunch  # noqa: E402
from src.dispersion.sheet_model import _y_slab, find_beta_from_sheet_admittance  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402
from src.hologram.leaky_forward import polygon_ray_distance  # noqa: E402
from src.materials import Substrate  # noqa: E402

F = 10.54e9
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
H = 1.575e-3
THETA0 = np.radians(30.0)
TARGET = 0.85
HEX_R = 0.15  # circumradius [m]


def hexagon(r=HEX_R):
    a = np.radians(np.arange(0, 360, 60) + 30.0)
    return np.stack([r * np.cos(a), r * np.sin(a)], axis=1)


def eta_and_directivity(alpha, poly, n_psi=360):
    psi = np.linspace(-np.pi, np.pi, n_psi, endpoint=False)
    rho_edge = polygon_ray_distance(poly, (0.0, 0.0), psi)
    eta = float(np.mean(1.0 - np.exp(-2.0 * alpha * rho_edge)))
    # aperture directivity of the illuminated polygon (scalar, continuous phase, obliquity cos(theta0))
    lam = 299792458.0 / F
    n = 160
    x = np.linspace(-HEX_R, HEX_R, n)
    xg, yg = np.meshgrid(x, x)
    inside = np.ones_like(xg, dtype=bool)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        inside &= ((b[0] - a[0]) * (yg - a[1]) - (b[1] - a[1]) * (xg - a[0])) >= -1e-12
    rho = np.maximum(np.hypot(xg, yg), 0.0032)
    a_n = (rho**-0.5 * np.exp(-alpha * rho))[inside]
    cell = (x[1] - x[0]) ** 2
    d = 4 * np.pi / lam**2 * (a_n.sum() ** 2 / np.sum(a_n**2)) * cell * np.cos(THETA0)
    return eta, float(10 * np.log10(d))


def main(out_dir):
    out = Path(out_dir)
    k0 = free_space_wavenumber(F)
    poly = hexagon()
    ys = np.arange(0.0, 16.01, 1.0)  # mS
    ds = np.concatenate([np.arange(0.0, 4.0, 0.5), np.arange(4.0, 20.1, 1.0)])
    grid = []
    launch_cache = {}
    for y in ys:
        beta = find_beta_from_sheet_admittance(1j * y * 1e-3, SUB, F).beta_rad_per_m
        y_down = 1j * y * 1e-3  # the sheet; depth metric below uses the TM0 total at beta
        if y not in launch_cache:
            launch_cache[y] = {
                name: PinLaunch(F, 2.2, H, Launcher(top), n_seg=16, y_sheet_s=1j * y * 1e-3).split()["launch_fraction"]
                for name, top in (("from_ground", 110 * 25.4e-6 - H), ("above_surface", 110 * 25.4e-6))
            }
        for d in ds:
            alpha = leakage_alpha_first_order_tensor(beta, y_down, d * 1e-3, THETA0, SUB, F)
            eta, d_ap = eta_and_directivity(alpha, poly)
            lf = launch_cache[y]["from_ground"]
            grid.append(
                {
                    "y_mS": float(y),
                    "delta_mS": float(d),
                    "beta_over_k0": beta / k0,
                    "alpha_np_per_m": alpha,
                    "eta_conv": eta,
                    "launch_fraction": lf,
                    "d_ap_dbi": d_ap,
                    "gain_estimate_dbi": d_ap + 10 * np.log10(max(eta * lf, 1e-12)),
                }
            )
    req = []
    for y in ys:
        col = [g for g in grid if g["y_mS"] == y]
        ok = [g for g in col if g["eta_conv"] >= TARGET]
        dreq = min((g["delta_mS"] for g in ok), default=None)
        beta_y = col[0]["beta_over_k0"] * k0
        y_slab = _y_slab(beta_y, 2 * np.pi * F, SUB.eps_complex, k0, H)
        depth = None if dreq is None else float(dreq * 1e-3 / abs(1j * y * 1e-3 + y_slab))
        req.append(
            {
                "modulation_depth_delta_over_Ydown": depth,
                "y_mS": float(y),
                "beta_over_k0": col[0]["beta_over_k0"],
                "delta_required_mS": min(g["delta_mS"] for g in ok) if ok else None,
                "launch_fraction_from_ground": launch_cache[y]["from_ground"],
                "launch_fraction_above_surface": launch_cache[y]["above_surface"],
            }
        )
    best_gain = max(grid, key=lambda g: g["gain_estimate_dbi"] if g["eta_conv"] >= TARGET else -1e9)
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"],
        "target_conversion": TARGET,
        "required_delta_vs_mean_loading": req,
        "best_gain_in_region": best_gain,
        "rows": grid,
    }
    (out / "required_map.json").write_text(json.dumps(res, indent=1))
    for r in req:
        print("y %4.1f mS beta/k0 %.4f  delta_req %s depth %s launch %.3f / %.3f" % (
            r["y_mS"], r["beta_over_k0"], r["delta_required_mS"], r["modulation_depth_delta_over_Ydown"] and round(r["modulation_depth_delta_over_Ydown"], 2), r["launch_fraction_from_ground"], r["launch_fraction_above_surface"]))
    print("best", best_gain)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        yy = sorted({g["y_mS"] for g in grid})
        dd = sorted({g["delta_mS"] for g in grid})
        z = np.array([[next(g["eta_conv"] for g in grid if g["y_mS"] == y and g["delta_mS"] == d) for y in yy] for d in dd])
        gz = np.array([[next(g["gain_estimate_dbi"] for g in grid if g["y_mS"] == y and g["delta_mS"] == d) for y in yy] for d in dd])
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
        c = ax[0].contourf(yy, dd, z, levels=[0, 0.1, 0.25, 0.5, 0.75, 0.85, 0.95, 1.0], cmap="viridis")
        ax[0].contour(yy, dd, z, levels=[TARGET], colors="w")
        ax[0].set_xlabel("mean loading  Im(Y_bar) [mS]")
        ax[0].set_ylabel("anisotropy delta [mS]")
        ax[0].set_title("converted fraction over the hexagon")
        fig.colorbar(c, ax=ax[0])
        c2 = ax[1].contourf(yy, dd, gz, levels=np.arange(8, 26, 2), cmap="magma")
        ax[1].contour(yy, dd, z, levels=[TARGET], colors="w")
        ax[1].set_xlabel("mean loading [mS]")
        ax[1].set_title("gain estimate [dBi] (launch fraction from the pin model)")
        fig.colorbar(c2, ax=ax[1])
        fig.tight_layout()
        fig.savefig(out / "required_map.png", dpi=130)
    except ImportError:
        pass


if __name__ == "__main__":
    main(sys.argv[1])
