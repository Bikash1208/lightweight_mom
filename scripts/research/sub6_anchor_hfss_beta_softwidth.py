"""Invert supplied HFSS anchor mode tracks and compare R150 sector soft widths."""

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from scripts.research.sub6_beta_alpha_map import alpha_direct
from scripts.research.sub6_meander_aperture import estimate

F0_GHZ = 5.8
P_M = 11.6e-3
C0 = 299792458.0
ALPHA_REQUESTED = 1.33


def beta_at_5p8(points):
    (phase0, f0, q0), (phase1, f1, q1) = points
    b0 = phase0 * C0 / (360 * P_M * f0 * 1e9)
    b1 = phase1 * C0 / (360 * P_M * f1 * 1e9)
    beta = b0 + (F0_GHZ - f0) * (b1 - b0) / (f1 - f0)
    return {"input_phase_deg_frequency_ghz_Q": points, "beta_points_over_k0": [b0, b1],
            "beta_at_5p8_over_k0": beta, "bracketed": f0 <= F0_GHZ <= f1,
            "method": "linear interpolation/extrapolation in (frequency, beta/k0)"}


def run(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tracks = {
        "across_x": beta_at_5p8([(85, 5.553, 2038), (95, 5.925, 1590)]),
        "along_y_bracketed_95_low": beta_at_5p8([(85, 5.743, 3694), (95, 6.140, 1755)]),
        "along_y_bracketed_95_high": beta_at_5p8([(85, 5.743, 3694), (95, 6.287, 2586)]),
        "along_y_high_95_low_extrapolated": beta_at_5p8([(85, 5.942, 1889), (95, 6.140, 1755)]),
        "along_y_high_extrapolated": beta_at_5p8([(85, 5.942, 1889), (95, 6.287, 2586)]),
    }
    bx = tracks["across_x"]["beta_at_5p8_over_k0"]
    by = tracks["along_y_bracketed_95_low"]["beta_at_5p8_over_k0"]
    by_alt = tracks["along_y_bracketed_95_high"]["beta_at_5p8_over_k0"]
    inversions = {
        "bracketed_95_low_pair": alpha_direct(bx, by),
        "bracketed_95_high_pair": alpha_direct(bx, by_alt),
        "approx_y_1p065": alpha_direct(bx, 1.065),
        "approx_y_1p070": alpha_direct(bx, 1.070),
    }
    y_for_133 = brentq(lambda y: alpha_direct(bx, y)["alpha_np_per_m"] - ALPHA_REQUESTED, 1.065, 1.070)
    inversions["y_needed_for_1p33"] = {"beta_y_over_k0": y_for_133, **alpha_direct(bx, y_for_133)}
    output_path = outdir / "hfss_beta_softwidth.json"
    previous = json.loads(output_path.read_text(encoding="utf-8")) if output_path.exists() else {}
    widths = previous.get("soft_width_cases", {}) if previous.get("soft_width_alpha_input_np_per_m") == ALPHA_REQUESTED else {}
    for width in (10, 20, 30, 40):
        if str(width) not in widths:
            widths[str(width)] = estimate(.150, False, angle_step_deg=1, alpha_np_m=ALPHA_REQUESTED, soft_boundary_mm=width)
    alpha_bracketed = inversions["bracketed_95_low_pair"]["alpha_np_per_m"]
    exact_widths = previous.get("bracketed_track_width_cases", {}) if abs(previous.get("bracketed_track_alpha_np_per_m", -1) - alpha_bracketed) < 1e-10 else {}
    for width in (20, 30, 40):
        if str(width) not in exact_widths:
            exact_widths[str(width)] = estimate(.150, False, angle_step_deg=1, alpha_np_m=alpha_bracketed, soft_boundary_mm=width)
    out = {"label": "MODEL_ONLY / NOT_HFSS_VALIDATED", "frequency_ghz": F0_GHZ,
           "period_mm": P_M * 1e3, "raw_mode_source": "user-supplied HFSS frequencies and Q",
           "tracks": tracks, "lid_correction_mm": 40, "inversions": inversions,
           "soft_width_alpha_input_np_per_m": ALPHA_REQUESTED,
           "soft_width_alpha_status": "user-requested approximate HFSS-pair map value; bracketed raw track maps to a different alpha",
           "soft_width_cases": widths, "bracketed_track_alpha_np_per_m": alpha_bracketed,
           "bracketed_track_width_cases_source": "along_y_bracketed_95_low",
           "bracketed_track_width_cases": exact_widths}
    output_path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    xs = [10, 20, 30, 40]
    ys = [widths[str(w)]["beams"][0]["sll_db"] for w in xs]
    fig, ax = plt.subplots(figsize=(6.5, 4.5), layout="constrained")
    ax.plot(xs, ys, "o-", label="A RHCP; B LHCP is symmetric")
    ax.axhline(-15, color="black", linestyle="--", label="SLL target")
    ax.set(xlabel="Soft sector boundary width (mm)", ylabel="Worst SLL (dB)",
           title="R150, alpha = 1.33 Np/m, 1° upper-hemisphere scan")
    ax.grid(alpha=.25)
    ax.legend()
    fig.savefig(outdir / "softwidth_sll.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"beta_x": bx, "beta_y_bracketed_95_low": by, "beta_y_bracketed_95_high": by_alt,
                      "alpha_bracketed_95_low": inversions["bracketed_95_low_pair"]["alpha_np_per_m"],
                      "alpha_bracketed_95_high": inversions["bracketed_95_high_pair"]["alpha_np_per_m"],
                      "beta_y_for_alpha_1p33": y_for_133,
                      "sll_by_width_alpha_1p33": dict(zip(xs, ys)),
                      "sll_by_width_bracketed_alpha": {w: exact_widths[str(w)]["beams"][0]["sll_db"] for w in (20, 30, 40)}}, indent=2), flush=True)


if __name__ == "__main__":
    import sys
    run(sys.argv[1])
