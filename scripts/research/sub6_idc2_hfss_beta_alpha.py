"""Invert user-supplied HFSS TM beta pairs into a MODEL_ONLY 30-degree alpha.

No capacitor, MoM geometry solve, or HFSS session is used. The beta values are
treated as measured inputs; their precision limits the reported alpha.
"""

import json
import math
from pathlib import Path

import numpy as np

from src.dispersion.leakage_model import leakage_alpha_first_order_tensor
from src.dispersion.sheet_model import _y_air, _y_slab, find_beta_from_sheet_admittance
from src.frequency import angular_frequency, free_space_wavenumber
from src.materials import Substrate

F = 5.8e9
SUB = Substrate("RT5880", 2.2, 0.0009, 3.175e-3)
K0 = free_space_wavenumber(F)
OMEGA = angular_frequency(F)
R = .150
INPUTS = {
    "MN_R5p6_W0p4_N2_T1p2_D3p0_E": {"beta_x_over_k0": 1.083, "beta_y_over_k0": 1.0675,
                                      "rounding_x": .0005, "rounding_y": .00005},
    "ST_R5p6_L6p72_W3p6": {"beta_x_over_k0": 1.099, "beta_y_over_k0": 1.0776,
                             "rounding_x": .0005, "rounding_y": .00005},
}


def sheet_b_from_beta(beta_over_k0):
    beta = beta_over_k0 * K0
    return float(-(_y_air(beta, OMEGA, K0) + _y_slab(beta, OMEGA, SUB.eps_complex, K0, SUB.thickness)).imag)


def estimate(beta_x_over_k0, beta_y_over_k0):
    bx, by = sheet_b_from_beta(beta_x_over_k0), sheet_b_from_beta(beta_y_over_k0)
    ybar = .5 * (bx + by)
    delta = .5 * abs(bx - by)
    beta_mean = find_beta_from_sheet_admittance(1j * ybar, SUB, F).beta_rad_per_m
    alpha = leakage_alpha_first_order_tensor(beta_mean, 1j * ybar, delta, math.radians(30), SUB, F)
    slab = _y_slab(beta_mean, OMEGA, SUB.eps_complex, K0, SUB.thickness)
    xs = [float(np.imag(1 / (1j * b + slab))) for b in (bx, by)]
    forward = [find_beta_from_sheet_admittance(1j * b, SUB, F).beta_rad_per_m / K0 for b in (bx, by)]
    if max(abs(forward[0] - beta_x_over_k0), abs(forward[1] - beta_y_over_k0)) > 1e-8:
        raise AssertionError("forward beta check failed after grounded-slab inversion")
    return {"sheet_admittance_x_mS": bx * 1e3, "sheet_admittance_y_mS": by * 1e3,
            "mean_sheet_susceptance_mS": ybar * 1e3, "delta_sheet_susceptance_mS": delta * 1e3,
            "beta_mean_over_k0": beta_mean / K0, "beta_forward_over_k0": forward,
            "surface_reactance_x_ohm": xs[0], "surface_reactance_y_ohm": xs[1],
            "M_t": abs(xs[0] - xs[1]) / abs(xs[0] + xs[1]),
            "alpha_np_per_m": alpha, "alpha_over_k0": alpha / K0,
            "ideal_conversion_R150": 1 - math.exp(-2 * alpha * R),
            "target_alpha_np_per_m": math.log(2) / (2 * R),
            "alpha_fraction_of_50pct_target": alpha / (math.log(2) / (2 * R))}


def run(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = {}
    for name, inp in INPUTS.items():
        bx, by = inp["beta_x_over_k0"], inp["beta_y_over_k0"]
        row = estimate(bx, by)
        corner = [estimate(bx + sx * inp["rounding_x"], by + sy * inp["rounding_y"])["alpha_np_per_m"]
                  for sx in (-1, 1) for sy in (-1, 1)]
        rows[name] = {"input_beta_over_k0": [bx, by], "input_rounding_half_step": [inp["rounding_x"], inp["rounding_y"]],
                      **row, "alpha_rounding_range_np_per_m": [min(corner), max(corner)]}
    n2, st = rows.keys()
    ratio = rows[n2]["alpha_np_per_m"] / rows[st]["alpha_np_per_m"]
    output = {"label": "MODEL_ONLY / NOT_HFSS_VALIDATED_LEAKAGE", "input_source": "user-supplied HFSS beta/k0 pairs",
              "frequency_hz": F, "k0_rad_per_m": K0, "theta_deg": 30,
              "substrate": {"eps_r": SUB.eps_r, "tan_delta": SUB.tan_delta, "thickness_m": SUB.thickness, "backing": "PEC"},
              "method": "For each real HFSS beta, B_sheet=-Im[Y_air_TM(beta)+Y_grounded_slab_TM(beta)]; mean B and half-difference delta; solve TM beta of mean jB, then leakage_alpha_first_order_tensor(beta_mean,jB,delta,theta=30).",
              "rows": rows, "N2_over_ST_alpha_ratio": ratio,
              "limits": ["Only TM phase constants are supplied; sheet conductance and real cell loss cannot be inferred.",
                         "Alpha is first-order radiating potential for a spatially rotated tensor hologram, not leakage of a uniform IDC cell.",
                         "No beta(f,r) data, so this single-frequency inversion cannot predict a radius-dependent flat resonance."]}
    (outdir / "result.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    lines = ["# Edge-to-edge two-tooth IDC: HFSS-beta inversion (MODEL_ONLY)", "",
             "Input β/k₀ values were supplied by the user for the 5.8 GHz HFSS cells. No HFSS instance or geometry solve was used. The grounded RT5880 slab is εr=2.2, tanδ=0.0009, h=3.175 mm over PEC; the 0.254 mm layer below the continuous ground is omitted. The separate x/y real TM phase constants are inverted with the same grounded-slab transverse-resonance relation used in `sheet_model.py`. The mean sheet beta is solved again from the mean susceptance before calling the first-order tensor leakage model at θ=30°.", "",
             "| Cell | β/k₀ x/y input | B_sheet x/y (mS) | B_mean / delta (mS) | β_mean/k₀ | M_t | α (Np/m) | α/k₀ | Ideal converted fraction R150 |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for name, row in rows.items():
        x, y = row["input_beta_over_k0"]
        lines.append(f"| {name} | {x:.4f}/{y:.4f} | {row['sheet_admittance_x_mS']:.4f}/{row['sheet_admittance_y_mS']:.4f} | {row['mean_sheet_susceptance_mS']:.4f}/{row['delta_sheet_susceptance_mS']:.4f} | {row['beta_mean_over_k0']:.5f} | {row['M_t']:.4f} | {row['alpha_np_per_m']:.4f} | {row['alpha_over_k0']:.6f} | {100*row['ideal_conversion_R150']:.2f}% |")
    lines += ["", f"N2_E / ST alpha ratio = **{ratio:.3f}**. N2_E alpha is **{rows[n2]['alpha_fraction_of_50pct_target']:.3f}×** the 50% target (2.3105 Np/m, alpha/k0≈0.0190). Thus this model does not support 50% conversion at R150. The ideal converted fraction assumes uniform leakage over all 150 mm and no other loss.", "",
              "Input rounding alone gives N2_E alpha approximately " + f"{rows[n2]['alpha_rounding_range_np_per_m'][0]:.4f}–{rows[n2]['alpha_rounding_range_np_per_m'][1]:.4f} Np/m" + f" and ST {rows[st]['alpha_rounding_range_np_per_m'][0]:.4f}–{rows[st]['alpha_rounding_range_np_per_m'][1]:.4f} Np/m. These are sensitivity ranges, not HFSS uncertainty bounds. Forward inversion reproduces both input beta pairs to <1e-8.", "",
              "The first-order alpha assumes a rotating-axis tensor modulation. A uniform periodic IDC cell has no such designed radiating harmonic. TM beta alone cannot determine sheet conductance, conductor loss, or a TE eigenmode. The observed flat resonance at p=14.5 mm and r=7.1 mm cannot be predicted or scaled in radius from these two 5.8 GHz beta pairs; that requires frequency-dependent dispersion at multiple radii for this IDC geometry. The earlier seven-leg meander model's branch end near 6.05 GHz is a different cell and is not evidence for the IDC radius trend. No radius dependence was fitted or invented.", "",
              "## HANDOFF", "", "- STATUS: Priority-one alpha estimate complete from user-supplied HFSS beta pairs.",
              "- Cell JSON: `result.json`; IDC coordinates are user-supplied, no geometry solve.",
              "- β/k₀ along/across: input x/y in table; mean beta solved from mean sheet susceptance.",
              "- α, M_t, loss ratio: table; N2/ST alpha ratio 0.723. Sheet conductance/loss unavailable from beta alone.",
              "- HFSS or MODEL_ONLY: HFSS beta inputs; derived alpha MODEL_ONLY.",
              "- Next request: HFSS beta(f,r) or eigenmode data at two radii to test resonance motion; slot-width/pitch screen after this gate.", ""]
    (outdir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return output


if __name__ == "__main__":
    import sys
    result = run(sys.argv[1])
    for name, row in result["rows"].items():
        print(name, "alpha", row["alpha_np_per_m"], "alpha/k0", row["alpha_over_k0"])
