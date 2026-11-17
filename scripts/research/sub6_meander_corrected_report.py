"""Analyze corrected seven-vertical-leg rows without running another solver."""

import json
import math
from pathlib import Path

import numpy as np
from src.dispersion.sheet_model import _y_slab
from src.frequency import angular_frequency, free_space_wavenumber
from src.materials import Substrate

SUB = Substrate("rt5880", 2.2, .0009, .003175)
F = 5.8e9
R = .150
TARGET_ALPHA = math.log(2) / (2 * R)


def _complex_matrix(pairs):
    return np.array([[complex(*z) for z in row] for row in pairs])


def conversion(row):
    alpha = row.get("alpha_radiating_potential_np_per_m")
    mt = row["anisotropy"].get("M_t_principal")
    if alpha is None or mt is None or alpha <= 0:
        return {"error": "no positive first-order alpha estimate"}
    y = _complex_matrix(row["sheet_y_siemens"])
    delta = row["anisotropy"]["delta_mS"] * 1e-3
    scale = math.sqrt(TARGET_ALPHA / alpha)
    beta_mean = .5 * (row["along_x"]["hybrid_TM_like"]["beta_over_k0"]
                     + row["across_y"]["hybrid_TM_like"]["beta_over_k0"]) * free_space_wavenumber(F)
    slab = _y_slab(beta_mean, angular_frequency(F), SUB.eps_complex, free_space_wavenumber(F), SUB.thickness)
    ybar = .5 * (y[0, 0].imag + y[1, 1].imag)
    eig = [ybar - scale * delta, ybar + scale * delta]
    xs = [float(np.imag(1 / (1j * v + slab))) for v in eig]
    mt_scaled = abs(xs[0] - xs[1]) / max(abs(xs[0] + xs[1]), 1e-30)
    return {"alpha_radiating_potential_np_per_m": alpha,
            "converted_fraction_R150_ideal": 1 - math.exp(-2 * alpha * R),
            "converted_fraction_active_R120_ideal": 1 - math.exp(-2 * alpha * .120),
            "alpha_for_50pct_R150_np_per_m": TARGET_ALPHA,
            "M_t_current": mt, "M_t_for_50pct_formal_square_law": mt * scale,
            "M_t_for_50pct_fixed_mean_tensor": mt_scaled,
            "delta_mS_for_50pct_fixed_mean_tensor": delta * scale * 1e3,
            "threshold_note": "Formal local model only: alpha scales as tensor-delta squared at fixed mean sheet and beta. M_t can be nonmonotonic or singular when principal reactances have opposite signs; no physical geometry or HFSS confirmation is implied."}


def analyze(folder):
    folder = Path(folder)
    names = ("baseline_7leg_w0p4", "7leg_w0p6", "9leg_w0p4", "5leg_w0p4")
    rows = []
    for name in names:
        path = folder / f"{name}_5.80GHz_n88.json"
        if path.exists():
            row = json.loads(path.read_text(encoding="utf-8"))
            rows.append({"geometry": name, "conversion": conversion(row), "tm_x": row["along_x"]["TM_projection"],
                         "tm_y": row["across_y"]["TM_projection"], "te_x": row["along_x"]["TE_projection"],
                         "te_y": row["across_y"]["TE_projection"], "M_t": row["anisotropy"].get("M_t_principal")})
    disp = []
    for path in folder.glob("baseline_7leg_w0p4_*GHz_n88.json"):
        row = json.loads(path.read_text(encoding="utf-8"))
        if row.get("status") != "COMPLETE":
            continue
        disp.append({"f_ghz": row["frequency_hz"] / 1e9,
                     "tm_x": row["along_x"]["TM_projection"]["roots_beta_rad_per_m"],
                     "te_x": row["along_x"]["TE_projection"]["roots_beta_rad_per_m"],
                     "tm_y": row["across_y"]["TM_projection"]["roots_beta_rad_per_m"],
                     "te_y": row["across_y"]["TE_projection"]["roots_beta_rad_per_m"]})
    disp.sort(key=lambda r: r["f_ghz"])
    slopes = {}
    for axis in ("x", "y"):
        track = [(d["f_ghz"], d[f"tm_{axis}"][0]) for d in disp if len(d[f"tm_{axis}"]) == 1]
        slopes[axis] = []
        for (f0, b0), (f1, b1) in zip(track, track[1:]):
            if f1 - f0 > .55:
                continue
            db_df = (b1 - b0) / ((f1 - f0) * 1e9)
            vg_c = (2 * math.pi / 299792458.0) / db_df if db_df != 0 else None
            slopes[axis].append({"f_mid_ghz": .5 * (f0 + f1), "f_interval_ghz": [f0, f1],
                                 "vg_over_c": vg_c, "flattening": vg_c is not None and 0 < vg_c < .1})
    result = {"label": "MODEL_ONLY / NOT_HFSS_VALIDATED", "R150_m": R,
              "active_R120_assumption": "30 mm absorber ring within R150; only for sensitivity",
              "rows": rows, "dispersion": disp, "tm_group_velocity_intervals": slopes}
    (folder / "analysis.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    finalist_names = ("baseline_7leg_w0p4", "7leg_w0p6", "9leg_w0p4")
    finalists = []
    for row in rows:
        if row["geometry"] not in finalist_names:
            continue
        src = json.loads((folder / f"{row['geometry']}_5.80GHz_n88.json").read_text(encoding="utf-8"))
        p = src["parameters_si"]
        finalists.append({"name": row["geometry"], "label": "MODEL_ONLY / NOT_HFSS_VALIDATED",
                          "period_mm": 11.6, "disc_radius_mm": p["diameter"] * 500,
                          "leg_count": p["leg_count"], "leg_pitch_mm": p["leg_pitch"] * 1000,
                          "leg_depth_mm": p["leg_depth"] * 1000, "slot_width_mm": p["slot_width"] * 1000,
                          "orientation_deg": 0,
                          "beta_over_k0_TM_x": row["tm_x"]["roots_beta_over_k0"],
                          "beta_over_k0_TM_y": row["tm_y"]["roots_beta_over_k0"],
                          "beta_over_k0_TE_x": row["te_x"]["roots_beta_over_k0"],
                          "beta_over_k0_TE_y": row["te_y"]["roots_beta_over_k0"],
                          "M_t": row["M_t"], **row["conversion"]})
    (folder / "finalists.json").write_text(json.dumps({"label": "MODEL_ONLY / NOT_HFSS_VALIDATED",
        "selection": "corrected seven-leg proposal plus slot-width and leg-count neighbors; no HFSS ranking",
        "cell_parameters": finalists}, indent=2) + "\n", encoding="utf-8")
    lines = ["# Corrected seven-vertical-leg meander study (MODEL_ONLY)", "",
             "The user's baseline is seven slot legs at x = -3.6, -2.4, -1.2, 0, 1.2, 2.4, 3.6 mm, each spanning y = -1.5 to +1.5 mm, 0.4 mm wide, with alternating end bridges. The disc radius is 5.6 mm; period 11.6 mm. The continuous cut clears the rim and its interior contains no copper. The substrate is grounded RT5880, 3.175 mm. `along_x` means propagation along global +x (across the leg pitch); `across_y` is global +y (parallel to the legs). The old horizontal-run screen is superseded.", "",
             "The earlier β used a TM-like hybrid tensor root, not a pure TM mode. Below are separate projected TM and TE roots of the bound-mode transverse-resonance equations. Off-diagonal tensor terms couple them in the physical mode; see each row's coupling ratio and hybrid character. A blank TE entry means no bound TE root in the searched interval.", "",
             "| Geometry | TM β/k₀ x | TM β/k₀ y | TE β/k₀ x | TE β/k₀ y | M_t principal | α rad potential (Np/m) | Ideal converted fraction R150 | Formal M_t for 50% |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        c = row["conversion"]
        def roots(key):
            return ", ".join(f"{x:.4f}" for x in row[key]["roots_beta_over_k0"]) or "—"
        lines.append(f"| {row['geometry']} | {roots('tm_x')} | {roots('tm_y')} | {roots('te_x')} | {roots('te_y')} | {row['M_t']:.3f} | {c.get('alpha_radiating_potential_np_per_m', float('nan')):.3f} | {100*c.get('converted_fraction_R150_ideal', float('nan')):.1f}% | {c.get('M_t_for_50pct_fixed_mean_tensor', float('nan')):.3f} |")
    lines += ["", "Converted fraction uses 1-exp(-2 α R) with R=0.150 m and the first-order *rotating-axis modulation potential* α. The required α for 50% is ln(2)/(2R)=2.3105 Np/m. The M_t numbers scale tensor anisotropy at fixed mean sheet and β; they are local mathematical thresholds, not design guarantees. M_t can exceed 1 when principal reactances have opposite signs. No finite-conductor loss, feed loss, spatially varying coupling, sector boundary, or HFSS validation is included.", "",
              "The R150 two-sector ideal array-factor worst local maximum is RHCP at θ≈2°, φ=180° and its symmetric LHCP partner at φ=0°, -13.62 dB. With ρ^-1/2 exp(-αρ) illumination, α from the corrected baseline, a 5.8 mm feed-floor radius and a 10 mm tanh sector transition, the worst local maximum is -13.32 dB near boresight; directivity is 18.84 dBi per beam, AR≈0.90 dB and XPD≈25.7 dB. The soft boundary is one explicit modeling choice, not a fabricated feed design. Full array-factor data are in `aperture_r150_rerun.json`.", ""]
    if len(disp) > 1:
        lines += ["## TM/TE dispersion", "", "| Frequency (GHz) | TM β x (rad/m) | TM β y (rad/m) | TE β x (rad/m) | TE β y (rad/m) |", "|---:|---:|---:|---:|---:|"]
        for d in disp:
            def val(key):
                return ", ".join(f"{v:.2f}" for v in d[key]) or "—"
            lines.append(f"| {d['f_ghz']:.2f} | {val('tm_x')} | {val('tm_y')} | {val('te_x')} | {val('te_y')} |")
        lines += ["", "Finite-difference v_g/c by frequency interval is recorded in `analysis.json`. A flattening claim requires a continuous tracked TM branch and small positive v_g/c; missing or multiple roots are flagged, not interpolated through.", ""]
    lines += ["## HANDOFF", "", "- STATUS: Corrected geometry and 5.8 GHz screen complete; dispersion/mesh status in result files.",
              "- Cell JSON: `finalists.json`.", "- β/k₀ along/across: projected TM and TE table above; hybrid values in row JSON.",
              "- α, M_t, loss ratio: α and M_t above; loss comparison remains first-order MODEL_ONLY.",
              "- HFSS or MODEL_ONLY: MODEL_ONLY; no switching claim.",
              "- Next request: HFSS cell dispersion and board-level CP check after mesh and TM-branch review.", ""]
    (folder / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return result


if __name__ == "__main__":
    import sys
    r = analyze(sys.argv[1])
    for row in r["rows"]:
        print(row["geometry"], row["conversion"])
    print("dispersion frequencies", [x["f_ghz"] for x in r["dispersion"]])
