"""MODEL_ONLY R150 two-sector aperture with radial leakage grading."""

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from scripts.research.sub6_meander_aperture import C0, F, K0, THETA_BEAM, estimate

R = 0.150
STEP = 0.005
SOFT = 0.010
FEED_FLOOR = 0.0058
ALPHA0 = 2.368319675555133  # corrected seven-leg baseline


def samples():
    g = np.arange(-R + STEP / 2, R, STEP)
    xx, yy = np.meshgrid(g, g)
    keep = xx**2 + yy**2 <= R**2
    x, y = xx[keep], yy[keep]
    rho = np.hypot(x, y)
    right = 0.5 * (1 + np.tanh(x / SOFT))
    return x, y, rho, np.sqrt(right), np.sqrt(1 - right)


def profile(rho, grade, power, family="outer_taper"):
    """Nonnegative alpha(r) capped by anchor alpha0, 0<=grade<=1."""
    if family == "outer_taper":
        alpha = ALPHA0 * (1 - grade * (rho / R)**power)
        integrated = ALPHA0 * (rho - grade * rho**(power + 1) / ((power + 1) * R**power))
        mean_factor = 1 - grade / (power + 1)
    elif family == "inner_taper":
        alpha = ALPHA0 * (1 - grade + grade * (rho / R)**power)
        integrated = ALPHA0 * ((1 - grade) * rho + grade * rho**(power + 1) / ((power + 1) * R**power))
        mean_factor = 1 - grade * power / (power + 1)
    else:
        raise ValueError(family)
    # sqrt(alpha) converts local remaining SW power to radiated aperture field.
    radial = np.sqrt(2 * alpha / np.maximum(rho, FEED_FLOOR)) * np.exp(-integrated)
    radial /= np.sqrt(np.mean(radial**2))
    total_integral = ALPHA0 * R * mean_factor
    return radial, 1 - math.exp(-2 * total_integral)


def field(x, y, radial, sector, theta_signed_deg):
    a = np.radians(theta_signed_deg)
    u = np.sin(a)
    phase = np.exp(1j * K0 * (u[:, None] - math.sin(THETA_BEAM)) * x[None, :])
    return (phase @ (radial * sector)) * np.cos(a)


def cut_metrics(x, y, rho, right, grade, power, family="outer_taper"):
    radial, conversion = profile(rho, grade, power, family)
    ang = np.arange(-60, 60.0001, .25)
    power_cut = abs(field(x, y, radial, right, ang))**2
    db = 10 * np.log10(power_cut / power_cut.max())
    main = int(np.argmax(power_cut))
    inner = find_peaks(power_cut[:main])[0]
    inner = inner[ang[inner] >= 0]
    first_inner = int(inner[-1]) if len(inner) else None
    minima = find_peaks(-power_cut)[0]
    minima_between = [float(ang[k]) for k in minima if first_inner is not None and first_inner < k < main]
    return {"family": family, "grade": grade, "power": power, "conversion": conversion,
            "main_deg": float(ang[main]), "inner_sidelobe_deg": float(ang[first_inner]) if first_inner is not None else None,
            "inner_sidelobe_db": float(db[first_inner]) if first_inner is not None else None,
            "separating_minima_deg": minima_between, "angles_deg": ang, "cut_db": db}


def main(out_dir="results/sub6_meander_model/corrected"):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sphere = json.loads((out_dir / "anchor_radial_grade_sphere.json").read_text(encoding="utf-8"))
    x, y, rho, right, left = samples()
    cuts = {key: cut_metrics(x, y, rho, right, grade, 1, "inner_taper")
            for key, grade in (("baseline", 0), ("grade_0p92", .92), ("grade_1p00", 1))}
    conversion_base = cuts["baseline"]["conversion"]
    rows = {}
    for key, cut in cuts.items():
        row = {k: v for k, v in cut.items() if k not in ("angles_deg", "cut_db")}
        row["conversion_cost_percentage_points"] = 100 * (conversion_base - cut["conversion"])
        row["alpha_center_np_per_m"] = ALPHA0 * (1 - cut["grade"])
        row["alpha_edge_np_per_m"] = ALPHA0
        row["upper_hemisphere_1deg"] = sphere[key]
        rows[key] = row
    result = {"label": "MODEL_ONLY / NOT_HFSS_VALIDATED", "anchor_alpha_np_per_m": ALPHA0,
              "anchor_source": "corrected seven-leg baseline first-order rotating-axis leakage potential",
              "aperture_radius_mm": R * 1e3, "sample_step_mm": STEP * 1e3,
              "soft_sector_boundary_mm": SOFT * 1e3, "feed_floor_mm": FEED_FLOOR * 1e3,
              "grade_family": "alpha(r)=alpha_anchor*(1-g+g*r/R), 0<=g<=1; radiated field amplitude proportional to sqrt(2*alpha(r)/rho)*exp(-integral_0^r alpha(s)ds)",
              "ideal_conversion": "1-exp(-2*integral_0^R alpha(r)dr); excludes feed, loss, absorber and finite board effects",
              "boresight_lobe": "Fine 0.25-degree x-plane cut confirms baseline first inner secondary local maximum at 1.25 degrees toward each beam, separated from main peak at 28.5 degrees by a shallow local minimum at 6.5 degrees, not a true null. The LHCP cut is the mirror image.",
              "rows": rows}
    (out_dir / "anchor_radial_grade_result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for key, label in (("baseline", "constant alpha"), ("grade_0p92", "linear grade g=0.92"), ("grade_1p00", "linear grade g=1.00")):
        c = cuts[key]
        ax.plot(c["angles_deg"], c["cut_db"], label=label)
    ax.axhline(-15, color="black", linestyle="--", linewidth=1, label="-15 dB target")
    ax.set(xlim=(-5, 45), ylim=(-35, 2), xlabel="Signed theta on A beam's x-plane (degrees)",
           ylabel="RHCP power relative to each case's main peak (dB)",
           title="R150 two-sector: boresight inner sidelobe and radial leakage grade")
    ax.grid(alpha=.25)
    ax.legend()
    fig.savefig(out_dir / "anchor_radial_grade_cut.png", dpi=180)
    plt.close(fig)
    print(json.dumps({k: {"conversion": v["conversion"], "cost_pp": v["conversion_cost_percentage_points"],
                          "cut_sll_db": v["inner_sidelobe_db"], "sphere_sll_db": v["upper_hemisphere_1deg"]["beams"][0]["sll_db"]}
                      for k, v in rows.items()}, indent=2), flush=True)


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "results/sub6_meander_model/corrected")
