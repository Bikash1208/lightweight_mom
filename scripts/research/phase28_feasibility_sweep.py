"""Phase 28 item 2: the feasibility sweep of a coded monopole-fed hologram on the Phase 26 forward model (`src/hologram/coded_feasibility.py`).

    python scripts/research/phase28_feasibility_sweep.py <out.json> [--workers 24] [--quick]

Run at 1 BLAS thread. Writes ONE json (the report `phase28_feasibility_report.py` is generated from it alone): anchor check against the committed
Phase 26 single-hexagon prediction, the pitch x aperture x steering x polarisation x bit-depth cases, the switching-contrast derivation, and the
control-hardware counts. Everything is MODEL_ONLY (idealised cells, the 6.4 mm / 10.54 GHz supercell table assumed at every pitch).
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from src.dispersion.leakage_model import floquet_complex_beta, two_level_fourier  # noqa: E402
from src.hologram import coded_feasibility as cf  # noqa: E402
from src.hologram import tensor_cp_design as t26  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable  # noqa: E402

PITCHES_MM = [4.0, 5.0, 6.4, 8.0]
VARIANTS = [(k, b) for k in ("cp", "lp", "lpx") for b in (None, 1, 2)]
CONSTANTS = {
    "frequency_ghz": 10.54,
    "frequency_ghz_not_available": 10.0,
    "gain_loss_tolerance_db": 1.0,
    "peak_refine_step_deg": 0.25,
    "coarse_grid_theta_phi_deg": [1.0, 2.0],
    "sidelobe_exclusion_cone_deg": 15.0,
    "ar_xpd_bound_db": cf.BOUND_DB,
    "launch_fraction_bare": cf.LAUNCH_BARE,
    "eta_design": cf.ETA_DESIGN,
    "taper_c": cf.TAPER_C,
    "diode_forward_current_ma": 10.0,
    "diode_forward_current_low_ma": 1.0,
    "diode_forward_voltage_v": 0.89,
    "supply_rail_v": 3.3,
    "shift_register_clock_mhz": 10.0,
    "table_extrapolated_flag_fraction": 0.05,
}


def _name(kind, bits):
    return f"{kind}_{'cont' if bits is None else f'{bits}bit'}"


def _specs(quick):
    out = []
    steer150 = [(0.0, 0.0)] + [(t, p) for t in (10.0, 20.0, 30.0, 40.0, 50.0) for p in (0.0, 45.0, 225.0)]
    steer225 = [(0.0, 0.0), (30.0, 225.0), (50.0, 225.0)]
    if quick:
        steer150, steer225 = [(0.0, 0.0), (30.0, 225.0)], [(30.0, 225.0)]
    for pitch in PITCHES_MM:
        for t, p in steer150:
            out.append({"pitch_mm": pitch, "aperture_r_mm": 150.0, "theta0_deg": t, "phi0_deg": p})
        for t, p in steer225:
            out.append({"pitch_mm": pitch, "aperture_r_mm": 225.0, "theta0_deg": t, "phi0_deg": p})
    return out


def run_case(spec):
    table = SupercellTable.from_files(t26.TABLE_GLOB)
    pitch = spec["pitch_mm"] * 1e-3
    poly = cf.hexagon_m(spec["aperture_r_mm"] * 1e-3)
    x, y, d, march = cf.design_cells(table, poly, pitch, spec["theta0_deg"], spec["phi0_deg"])
    green = cf.green_function()
    beam = (spec["theta0_deg"], spec["phi0_deg"])
    eta = march["eta_conv"]
    row = dict(spec)
    row.update(
        {
            "n_cells": int(x.size),
            "eta_conv": float(eta),
            "table_extrapolated_fraction": cf.table_extrapolation_fraction(march, table),
            "variants": {},
            "bias": {},
        }
    )
    m_lp = cf.lp_modulation(d, x, y, *beam)
    for kind, bits in VARIANTS:
        if kind == "cp":
            idx, axis = cf.cp_states(d, bits)
            jx, jy = cf.cp_currents(d, axis)
        else:
            idx, mq = cf.lp_states(m_lp, bits)
            jx, jy = (cf.lp_currents if kind == "lp" else cf.lpx_currents)(d, mq)
        r = cf.analyse(jx, jy, x, y, green, pitch, "cp" if kind == "cp" else "lp", beam)
        r["gain_dbi"] = cf.launch_gain_dbi(r["d_ap_dbi"], eta)
        row["variants"][_name(kind, bits)] = r
        if bits is not None and kind != "lpx":
            frac = cf.state_bit_stats(idx, bits)
            b = cf.bias_budget(x.size, bits, frac)
            b["fraction_bits_forward"] = frac
            row["bias"][_name(kind, bits)] = b
    return row


def run_anchor():
    table = SupercellTable.from_files(t26.TABLE_GLOB)
    x, y, d, march = cf.design_cells(table, cf.hexagon_m(0.150), 6.4e-3, 30.0, 225.0)
    got = cf.analyse(d["jx"], d["jy"], x, y, cf.green_function(), 6.4e-3, "cp", (30.0, 225.0))
    ref = json.loads((ROOT / "results/phase26/tensor_cp_hexagon.json").read_text(encoding="utf-8"))["selected"]
    m = ref["metrics_ap"]
    rows = [
        ("d_ap_rhcp_dbi", got["d_ap_coarse_grid_dbi"], m["directivity_ap_rhcp_dbi"], 1e-6),
        (
            "gain_launch_bare_dbi",
            cf.launch_gain_dbi(got["d_ap_coarse_grid_dbi"], march["eta_conv"]),
            ref["gain_rhcp_dbi_launch_bare"],
            1e-6,
        ),
        ("ar_at_peak_db", got["ar_at_coarse_peak_db"], m["ar_at_peak_db"], 1e-6),
        ("xpd_at_peak_db", got["xpd_at_coarse_peak_db"], m["xpd_at_peak_db"], 1e-6),
        ("rhcp_sll_db", got["sll_coarse_grid_db"], m["rhcp_sll_db"], 1e-6),
        ("eta_conv", march["eta_conv"], ref["eta_conv"], 1e-12),
    ]
    out = {
        "source": "results/phase26/tensor_cp_hexagon.json (selected design: hand +1, eta_design 0.9, taper 0.5)",
        "rows": [],
    }
    for name, mine, theirs, tol in rows:
        out["rows"].append(
            {
                "quantity": name,
                "this_pipeline": mine,
                "phase26": theirs,
                "abs_diff": abs(mine - theirs),
                "tolerance": tol,
                "within_tolerance": abs(mine - theirs) <= tol,
            }
        )
    out["refined_d_ap_dbi"] = got["d_ap_dbi"]
    out["refinement_gain_over_coarse_db"] = got["d_ap_dbi"] - got["d_ap_coarse_grid_dbi"]
    out["refined_pointing_error_deg"] = got["pointing_error_deg"]
    out["all_within_tolerance"] = all(r["within_tolerance"] for r in out["rows"])
    return out


def _gain_with_cap(args):
    cap, theta0, phi0 = args
    table = (
        cf.capped_table(SupercellTable.from_files(t26.TABLE_GLOB), cap)
        if cap is not None
        else SupercellTable.from_files(t26.TABLE_GLOB)
    )
    x, y, d, march = cf.design_cells(table, cf.hexagon_m(0.150), 6.4e-3, theta0, phi0)
    r = cf.analyse(d["jx"], d["jy"], x, y, cf.green_function(), 6.4e-3, "cp", (theta0, phi0))
    return {
        "alpha_cap_np_per_m": cap,
        "gain_dbi": cf.launch_gain_dbi(r["d_ap_dbi"], march["eta_conv"]),
        "d_ap_dbi": r["d_ap_dbi"],
        "eta_conv": march["eta_conv"],
        "beta_over_k0_mean": float(np.mean(d["beta_over_k0"][d["alpha"] > 0])),
        "alpha_max_used_np_per_m": float(d["alpha"].max()),
    }


def run_contrast(pool, steerings=((0.0, 0.0), (30.0, 225.0), (50.0, 225.0))):
    table = SupercellTable.from_files(t26.TABLE_GLOB)
    amax = float(np.nanmax(table.alpha))
    caps = [float(c) for c in np.round(np.geomspace(0.25, amax, 20), 4)]
    out = {"table_alpha_max_np_per_m": amax, "table_alpha_min_np_per_m": float(np.nanmin(table.alpha)), "steerings": []}
    for theta0, phi0 in steerings:
        jobs = [(None, theta0, phi0)] + [(c, theta0, phi0) for c in caps]
        res = list(pool.map(_gain_with_cap, jobs))
        base = res[0]
        curve = res[1:]
        for r in curve:
            r["loss_vs_uncapped_db"] = base["gain_dbi"] - r["gain_dbi"]
        # smallest cap whose loss is <= 1 dB: interpolate in log(cap), then re-evaluate at the interpolated cap
        tol = CONSTANTS["gain_loss_tolerance_db"]
        cap_star = None
        for lo, hi in zip(curve[:-1], curve[1:], strict=True):
            if lo["loss_vs_uncapped_db"] > tol >= hi["loss_vs_uncapped_db"]:
                f = (lo["loss_vs_uncapped_db"] - tol) / (lo["loss_vs_uncapped_db"] - hi["loss_vs_uncapped_db"])
                cap_star = float(
                    math.exp(
                        math.log(lo["alpha_cap_np_per_m"])
                        + f * (math.log(hi["alpha_cap_np_per_m"]) - math.log(lo["alpha_cap_np_per_m"]))
                    )
                )
                break
        if cap_star is None and curve[0]["loss_vs_uncapped_db"] <= tol:
            cap_star = curve[0]["alpha_cap_np_per_m"]
        entry = {
            "theta0_deg": theta0,
            "phi0_deg": phi0,
            "uncapped": base,
            "curve": curve,
            "alpha_cap_for_1db_np_per_m": cap_star,
        }
        if cap_star is not None:
            chk = _gain_with_cap((cap_star, theta0, phi0))
            chk["loss_vs_uncapped_db"] = base["gain_dbi"] - chk["gain_dbi"]
            entry["check_at_cap_star"] = chk
            beta_over_k0 = base["beta_over_k0_mean"]
            entry["contrast"] = {f"bits_{b}": cf.contrast_from_alpha(cap_star, beta_over_k0, theta0, b) for b in (1, 2)}
            entry["floquet_check"] = floquet_check(cap_star, beta_over_k0, theta0, entry["contrast"]["bits_1"])
        out["steerings"].append(entry)
    return out


def floquet_check(alpha_target, beta_over_k0, theta0_deg, contrast_1bit):
    """Exact truncated-Floquet complex beta of a binary (duty 0.5) sheet with the first-order full-swing `dB` for the TM (zeta = 0) case, collinear."""
    k0 = cf.free_space_wavenumber(cf.F)
    ybar = 1j * contrast_1bit["y_mean_ms"] * 1e-3
    beta = contrast_1bit["beta_over_k0_mean"] * k0
    from src.dispersion.leakage_model import leakage_alpha_first_order

    th = np.radians(theta0_deg)
    k_s = leakage_alpha_first_order(beta, ybar, 1.0, th, 0.0, cf.SUB, cf.F)
    g = contrast_1bit["scalar"]["first_harmonic_gain_g"]
    db = 4.0 * math.sqrt(alpha_target / k_s) / g
    y0, y1 = ybar - 0.5j * db, ybar + 0.5j * db
    period = 2 * np.pi / (beta - k0 * math.sin(th))
    with np.errstate(all="ignore"):
        b_exact, _ = floquet_complex_beta(two_level_fourier(y0, y1, 6), period, 0.0, cf.SUB, cf.F, beta)
    if not np.isfinite(
        b_exact
    ):  # broadside: the n = -1 harmonic sits at the Bragg condition, the Floquet balance is singular
        return {
            "zeta_deg": 0.0,
            "d_b_full_swing_ms": db * 1e3,
            "alpha_first_order_target_np_per_m": alpha_target,
            "alpha_floquet_exact_np_per_m": None,
            "ratio_exact_over_first_order": None,
            "beta_over_k0_floquet": None,
            "note": "NOT_AVAILABLE: the Floquet balance has no finite root at broadside (open stop band)",
        }
    return {
        "zeta_deg": 0.0,
        "d_b_full_swing_ms": db * 1e3,
        "alpha_first_order_target_np_per_m": alpha_target,
        "alpha_floquet_exact_np_per_m": float(-b_exact.imag),
        "ratio_exact_over_first_order": float(-b_exact.imag / alpha_target),
        "beta_over_k0_floquet": float(b_exact.real / k0),
    }


def diode_terminal():
    """Terminal impedance swing of the repo's diode models at the design frequency (datasheet values, PROVISIONAL): NOT the cell-level contrast."""
    from src.ris import diodes

    out = {}
    for name, on, off in (
        ("SMP1345-079LF", diodes.smp1345_on("10mA"), diodes.smp1345_off()),
        ("SMP1345-079LF_1mA", diodes.smp1345_on("1mA"), diodes.smp1345_off()),
        ("MADP-000907-14020", diodes.madp_000907_on(), diodes.madp_000907_off()),
    ):
        z_on, z_off = on.impedance(cf.F), off.impedance(cf.F)
        out[name] = {
            "z_on_ohm": [z_on.real, z_on.imag],
            "z_off_ohm": [z_off.real, z_off.imag],
            "delta_x_off_minus_on_ohm": z_off.imag - z_on.imag,
            "r_on_ohm": z_on.real,
        }
    out["note"] = (
        "MADP-000907-14020 package inductance is NOT_AVAILABLE from the datasheet (L_s = 0 here); datasheet values are PROVISIONAL"
    )
    return out


def derive(rows):
    """Quantisation loss vs the same-family continuous design, scan loss vs broadside, per case; and the design-point summary."""
    key = lambda r: (r["pitch_mm"], r["aperture_r_mm"], r["theta0_deg"], r["phi0_deg"])  # noqa: E731
    by = {key(r): r for r in rows}
    for r in rows:
        b = by.get((r["pitch_mm"], r["aperture_r_mm"], 0.0, 0.0))
        for name, v in r["variants"].items():
            cont = r["variants"][name.split("_")[0] + "_cont"]
            v["quantisation_loss_db"] = cont["gain_dbi"] - v["gain_dbi"]
            v["scan_loss_db"] = None if b is None else b["variants"][name]["gain_dbi"] - v["gain_dbi"]
    summary = []
    for r in rows:
        if r["aperture_r_mm"] == 150.0 and r["theta0_deg"] == 30.0 and r["phi0_deg"] == 225.0:
            for name, v in r["variants"].items():
                summary.append(
                    {
                        "pitch_mm": r["pitch_mm"],
                        "variant": name,
                        "gain_dbi": v["gain_dbi"],
                        "quantisation_loss_db": v["quantisation_loss_db"],
                        "sll_db": v["sll_db"],
                        "pointing_error_deg": v["pointing_error_deg"],
                    }
                )
    return summary


def main(out_path, workers=24, quick=False):
    import provenance as prov

    t0 = time.time()
    specs = _specs(quick)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        print(f"{len(specs)} cases on {workers} workers", flush=True)
        rows = []
        for i, r in enumerate(pool.map(run_case, specs), 1):
            rows.append(r)
            print(
                f"  {i}/{len(specs)} p={r['pitch_mm']} R={r['aperture_r_mm']} th={r['theta0_deg']} ph={r['phi0_deg']} n={r['n_cells']}  ({time.time() - t0:.0f} s)",
                flush=True,
            )
        anchor = run_anchor()
        contrast = run_contrast(pool, ((30.0, 225.0),) if quick else ((0.0, 0.0), (30.0, 225.0), (50.0, 225.0)))
    summary = derive(rows)
    doc = {
        "label": "PHASE28_ITEM2_FEASIBILITY",
        "labels": [
            "MODEL_ONLY",
            "IDEALISED_CELLS",
            "NOT_HFSS_VALIDATED",
            "SUPERCELL_TABLE_6P4MM_10P54GHZ_ASSUMED_AT_ALL_PITCHES",
        ],
        "constants": CONSTANTS,
        "frequency_10p0_ghz": "NOT_AVAILABLE_FROM_CURRENT_DATA (the supercell table exists at 10.54 GHz only)",
        "anchor": anchor,
        "cases": rows,
        "design_point_summary": summary,
        "contrast": contrast,
        "diode_terminal": diode_terminal(),
        "seconds": time.time() - t0,
    }
    Path(out_path).write_text(json.dumps(doc, indent=1, sort_keys=True), encoding="utf-8", newline="\n")
    rec = prov.finish(
        prov.provenance({"tool": "phase28_feasibility_sweep.py", "n_cases": len(rows), "quick": quick}, workers=workers)
    )
    Path(str(out_path).replace(".json", ".provenance.json")).write_text(
        json.dumps(rec, indent=1, sort_keys=True), encoding="utf-8", newline="\n"
    )
    print("wrote", out_path, f"{time.time() - t0:.0f} s; anchor within tolerance: {anchor['all_within_tolerance']}")


if __name__ == "__main__":
    args = sys.argv[1:]
    w = int(args[args.index("--workers") + 1]) if "--workers" in args else 24
    main(args[0], workers=w, quick="--quick" in args)
