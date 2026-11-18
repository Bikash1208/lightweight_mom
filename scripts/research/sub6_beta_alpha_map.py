"""Map HFSS PEC-lid TM beta pairs to open-air first-order tensor leakage.

The finite PEC lid is used only for the beta-to-sheet inversion. The physical
open-air slab is then used for the mean TM wave and radiating harmonic. A closed
PEC lid has zero far-field leakage; this map estimates leakage after its removal.
"""

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator

from src.constants import EPS0
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor
from src.dispersion.sheet_model import _y_air, _y_slab, find_beta_from_sheet_admittance
from src.frequency import angular_frequency, free_space_wavenumber
from src.materials import Substrate


F = 5.8e9
SUB = Substrate("RT5880", 2.2, 0.0009, 3.175e-3)
LID_H = 40e-3  # actual air_top: 40 mm above top copper; sheet-plane offset <= 0.035 mm
THETA = math.radians(30)
K0 = free_space_wavenumber(F)
OMEGA = angular_frequency(F)


def lid_air_admittance(beta):
    kappa = np.sqrt(np.asarray(beta, dtype=float) ** 2 - K0**2)
    return 1j * OMEGA * EPS0 / (kappa * np.tanh(kappa * LID_H))


def sheet_b_lid(beta_over_k0):
    beta = np.asarray(beta_over_k0, dtype=float) * K0
    return -np.imag(lid_air_admittance(beta) + _y_slab(beta, OMEGA, SUB.eps_complex, K0, SUB.thickness))


def alpha_direct(bx, by):
    b_x, b_y = float(sheet_b_lid(bx)), float(sheet_b_lid(by))
    b_bar, delta = (b_x + b_y) / 2, abs(b_x - b_y) / 2
    beta_mean = find_beta_from_sheet_admittance(1j * b_bar, SUB, F).beta_rad_per_m
    alpha = leakage_alpha_first_order_tensor(beta_mean, 1j * b_bar, delta, THETA, SUB, F)
    return {"B_x_mS": b_x * 1e3, "B_y_mS": b_y * 1e3,
            "B_mean_mS": b_bar * 1e3, "delta_mS": delta * 1e3,
            "beta_mean_open_over_k0": beta_mean / K0, "alpha_np_per_m": alpha,
            "alpha_over_k0": alpha / K0}


def run(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    axis = np.linspace(1.005, 1.32, 316)
    bx, by = np.meshgrid(axis, axis, indexing="xy")
    bs_x, bs_y = sheet_b_lid(bx), sheet_b_lid(by)
    bbar = (bs_x + bs_y) / 2
    delta = abs(bs_x - bs_y) / 2
    # The open-air mean beta and alpha/delta^2 depend only on B_mean. Interpolate
    # those 1-D quantities, then check the resulting map against direct calls.
    bvals = np.linspace(float(bbar.min()), float(bbar.max()), 1200)
    beta_open = np.empty_like(bvals)
    coeff = np.empty_like(bvals)
    for i, b in enumerate(bvals):
        beta = find_beta_from_sheet_admittance(1j * b, SUB, F).beta_rad_per_m
        beta_open[i] = beta / K0
        coeff[i] = leakage_alpha_first_order_tensor(beta, 1j * b, 1e-3, THETA, SUB, F) / (1e-3**2 * K0)
    alpha_over_k0 = PchipInterpolator(bvals, coeff)(bbar) * delta**2
    beta_mean = PchipInterpolator(bvals, beta_open)(bbar)
    points = {
        "N2_E": (1.083, 1.0675),
        "ST": (1.099, 1.0776),
    }
    point_results = {}
    for name, pair in points.items():
        direct = alpha_direct(*pair)
        interp = float(PchipInterpolator(bvals, coeff)((direct["B_mean_mS"] * 1e-3)) * (direct["delta_mS"] * 1e-3)**2)
        if abs(interp - direct["alpha_over_k0"]) > 1e-5:
            raise AssertionError(f"map interpolation disagrees at {name}")
        point_results[name] = {"beta_x_over_k0": pair[0], "beta_y_over_k0": pair[1], **direct}
    # More checks near both requested contour levels, not only the two HFSS pairs.
    for target in (0.010, 0.019):
        idx = np.unravel_index(np.argmin(abs(alpha_over_k0 - target)), alpha_over_k0.shape)
        direct = alpha_direct(float(bx[idx]), float(by[idx]))["alpha_over_k0"]
        if abs(float(alpha_over_k0[idx]) - direct) > 1e-5:
            raise AssertionError(f"map interpolation disagrees near contour {target}")
    with (out / "beta_alpha_map.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["beta_x_over_k0", "beta_y_over_k0", "B_x_mS", "B_y_mS", "B_mean_mS", "delta_mS", "beta_mean_open_over_k0", "alpha_over_k0"])
        for values in zip(bx.flat, by.flat, bs_x.ravel() * 1e3, bs_y.ravel() * 1e3,
                          bbar.ravel() * 1e3, delta.ravel() * 1e3, beta_mean.flat, alpha_over_k0.flat):
            w.writerow([f"{v:.9g}" for v in values])
    fig, ax = plt.subplots(figsize=(8, 7), layout="constrained")
    mesh = ax.pcolormesh(bx, by, np.minimum(alpha_over_k0, 0.05), shading="auto", cmap="viridis", vmin=0, vmax=0.05)
    cs = ax.contour(bx, by, alpha_over_k0, levels=[0.010, 0.019], colors=["white", "#ff6b57"], linewidths=2)
    ax.clabel(cs, fmt={0.010: "0.010", 0.019: "0.019"}, fontsize=10)
    for name, row in point_results.items():
        ax.plot(row["beta_x_over_k0"], row["beta_y_over_k0"], "o", markersize=7, markeredgecolor="black", markerfacecolor="yellow")
        ax.annotate(name, (row["beta_x_over_k0"], row["beta_y_over_k0"]), xytext=(7, 6), textcoords="offset points", color="white", fontsize=10, weight="bold")
    ax.plot(axis, axis, "--", color="0.7", lw=1)
    ax.set(xlim=(axis.min(), axis.max()), ylim=(axis.min(), axis.max()), xlabel=r"HFSS $\beta_x/k_0$", ylabel=r"HFSS $\beta_y/k_0$",
           title=r"RT5880 3.175 mm, 5.8 GHz, 30°: inferred open-air $\alpha/k_0$")
    ax.set_aspect("equal")
    fig.colorbar(mesh, ax=ax, label=r"$\alpha/k_0$ (MODEL_ONLY, clipped at 0.05)")
    fig.text(0.5, 0.005, "HFSS beta inversion includes PEC lid at air_top = 40 mm; leakage evaluation removes lid.", ha="center", fontsize=9)
    fig.savefig(out / "beta_alpha_map.png", dpi=180)
    plt.close(fig)
    summary = {"label": "MODEL_ONLY / NOT_HFSS_VALIDATED_LEAKAGE", "frequency_hz": F,
               "substrate": {"eps_r": SUB.eps_r, "tan_delta": SUB.tan_delta, "thickness_mm": SUB.thickness * 1e3},
               "lid_height_mm_for_beta_inversion": LID_H * 1e3, "theta_deg": math.degrees(THETA),
               "map_axis_min_max_step": [float(axis.min()), float(axis.max()), float(axis[1] - axis[0])],
               "contours_alpha_over_k0": [0.010, 0.019], "points": point_results,
               "interpretation": "PEC lid used to infer sheet B from HFSS beta; open air used to evaluate radiating first-order tensor alpha. Actual PEC-lidded cell has zero far-field radiative leakage.",
               "limitations": "Rotating-axis first-order tensor potential; not uniform-cell leakage or HFSS measured attenuation."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    import sys
    run(sys.argv[1])
