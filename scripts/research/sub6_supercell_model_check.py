"""MODEL_ONLY frequency-specific beta->B->tensor-alpha check; no AEDT calls."""

import csv
import json
import math
from pathlib import Path

import numpy as np

from src.constants import EPS0
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor
from src.dispersion.sheet_model import _y_slab, find_beta_from_sheet_admittance
from src.frequency import angular_frequency, free_space_wavenumber
from src.materials import Substrate


ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"
SUB = Substrate("RT5880", 2.2, 0.0009, 3.175e-3)
P_M = 0.0116
C0 = 299792458.0
LID_M = 0.040
ANGLES = (25, 30, 35, 40)


def read_modes(cell, direction, phase):
    if direction == "across" and phase == 85:
        file = R / f"{cell}_PEC_modes.json"
    else:
        file = R / f"{cell}_PEC_{direction}_{phase}_air20.035_modes.json"
    return json.loads(file.read_text(encoding="utf-8"))["modes"]


def beta_pair(cell, freq_ghz, x_modes=(1, 1), y_modes=(1, 1)):
    def axis(direction, indices):
        a = read_modes(cell, direction, 85)[indices[0] - 1]
        b = read_modes(cell, direction, 95)[indices[1] - 1]
        fa, fb = a["f_ghz"], b["f_ghz"]
        ba = 85 * C0 / (360 * P_M * fa * 1e9)
        bb = 95 * C0 / (360 * P_M * fb * 1e9)
        value = ba + (freq_ghz - fa) * (bb - ba) / (fb - fa)
        return {"beta_over_k0": value, "f_endpoints_ghz": [fa, fb],
                "mode_indices": list(indices), "bracketed": fa <= freq_ghz <= fb}
    return {"x": axis("across", x_modes), "y": axis("along", y_modes)}


def estimate(freq_ghz, pair, theta_deg):
    f = freq_ghz * 1e9
    k0 = free_space_wavenumber(f)
    omega = angular_frequency(f)

    def sheet_b(beta_ratio):
        beta = beta_ratio * k0
        kappa = np.sqrt(beta * beta - k0 * k0)
        y_air_pec = 1j * omega * EPS0 / (kappa * np.tanh(kappa * LID_M))
        return float(-(y_air_pec + _y_slab(beta, omega, SUB.eps_complex,
                                          k0, SUB.thickness)).imag)

    bx = sheet_b(pair["x"]["beta_over_k0"])
    by = sheet_b(pair["y"]["beta_over_k0"])
    bmean = (bx + by) / 2
    delta = abs(bx - by) / 2
    beta_mean = find_beta_from_sheet_admittance(1j * bmean, SUB, f).beta_rad_per_m
    alpha = leakage_alpha_first_order_tensor(beta_mean, 1j * bmean,
                                             delta, math.radians(theta_deg), SUB, f)
    return {"B_x_mS": bx * 1e3, "B_y_mS": by * 1e3,
            "B_mean_mS": bmean * 1e3, "delta_mS": delta * 1e3,
            "beta_mean_open_over_k0": beta_mean / k0,
            "alpha_np_m": alpha}


def main():
    cases = [
        ("ANCHOR", "main", (1, 1), (1, 1), 5.8, None),
        ("ANCHOR", "alternate_y", (1, 1), (1, 2), 5.8, None),
        ("B2p2_D3p0", "main", (1, 1), (1, 1), 5.8, None),
        ("B2p2_D3p6", "x_low", (1, 1), (1, 1), 5.8, None),
        ("B2p2_D3p6", "x_high", (2, 2), (1, 1), 5.8, None),
    ]
    for f, measured in ((5.594682258439221, 1.1119),
                        (5.681668085433831, 1.4529971755240307),
                        (5.752876246511, 1.8311)):
        cases.append(("ANCHOR", "main_at_supercell_f", (1, 1), (1, 1), f, measured))
    rows = []
    for cell, branch, xm, ym, f, measured in cases:
        pair = beta_pair(cell, f, xm, ym)
        for theta in ANGLES:
            est = estimate(f, pair, theta)
            rows.append({"label": "MODEL_ONLY", "cell": cell, "branch": branch,
                         "f_ghz": f, "theta_deg": theta,
                         "beta_x_over_k0": pair["x"]["beta_over_k0"],
                         "beta_y_over_k0": pair["y"]["beta_over_k0"],
                         "x_bracketed": pair["x"]["bracketed"],
                         "y_bracketed": pair["y"]["bracketed"],
                         **est, "hfss_alpha_np_m": measured,
                         "model_over_hfss": est["alpha_np_m"] / measured if measured else None})
    out = {"method": "40mm PEC-lid TM beta inversion, open-air mean-sheet tensor leakage",
           "source": "own AEDT 2025.2 lidded 85/95 degree rows; linear beta/k0(f) interpolation/extrapolation",
           "rows": rows}
    (R / "model_check.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    with (R / "model_check.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        if row["theta_deg"] == 30:
            print(row["cell"], row["branch"], f"f={row['f_ghz']:.6f}",
                  f"beta={row['beta_x_over_k0']:.6f}/{row['beta_y_over_k0']:.6f}",
                  f"alpha={row['alpha_np_m']:.6f}",
                  f"ratio={row['model_over_hfss']}",
                  f"bracket={row['x_bracketed']}/{row['y_bracketed']}")


if __name__ == "__main__":
    main()
