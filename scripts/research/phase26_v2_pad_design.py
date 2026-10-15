"""Phase 26 V2: graded launch pad around the pin of the single-hex tensor-CP board + predictions (written BEFORE the HFSS run).

    python scripts/research/phase26_v2_pad_design.py <out.json>

Change w.r.t. the current run (the ONLY one): cells with centre radius r < 30 mm are replaced by isotropic
cross-slotted circles (D 6.0 mm, both slots L 4.8 mm, width W(r) graded) -- 4-fold symmetric, no leakage, no phase
design -- and the feed-cell keep-out becomes the one feed cell (r < 3.2 mm; the lattice cells at 6.4 mm are
now filled). Every cell with r >= 30 mm is the current design's, unchanged (slot sizes and angles). The reference
phase gains a constant offset from the pad's own beta only (circular pad): the pattern is unaffected.
Loading grade: Ybar(r) = Ybar_high for r <= 15 mm, linear in r to the hologram's local mean at r = 30 mm; W(r)
from the Ybar(W) curve of the cross cell at kt = 1.06 k0 (`v2_pad_loading.json`/this script).

Isotropic loading available in the stable family is limited: an unslotted D 6.0 circle has 4.4 mS (LESS than the
hologram's 5.4-7.8 mS slotted cells); the cross-slotted circle 4.8 x 2.0 (bridge 0.43 mm) gives ~7.5 mS.
Predictions: gain, S11, grazing+back, LHCP share; launch fraction is bracketed (the model's sheet is infinite).
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase26_launcher_redesign as lr  # noqa: E402
import phase26_tensor_cp_hexagon as cp  # noqa: E402
import phase26_v2_pad_loading as pl  # noqa: E402
from src.dispersion.slotted_cells import CrossSlottedCircularPatch  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays  # noqa: E402

R_PAD, R_HIGH = 30.0, 15.0
L_CROSS = 4.8
W_GRID = (0.8, 1.0, 1.2, 1.6, 2.0)
YBAR_SLOT = {"4.4x1.2": 5.43, "4.6x1.6": 6.52, "4.8x1.6": 6.71, "4.8x2": 7.84}  # v2_pad_loading.json, kt 1.06 k0, 8 angles


def main(out):
    cells = json.loads((ROOT / "results/phase26/tensor_cp_hexagon_cells.json").read_text())
    x, y = np.array(cells["x_mm"]), np.array(cells["y_mm"])
    L, W, A = np.array(cells["slot_L_mm"]), np.array(cells["slot_W_mm"]), np.array(cells["slot_angle_deg"])
    r = np.hypot(x, y)
    # hologram local mean loading at the pad edge: cells in 30..50 mm
    ann = (r >= 30.0) & (r < 50.0)
    ylocal = float(np.mean([YBAR_SLOT[f"{l:g}x{w:g}"] for l, w in zip(L[ann], W[ann], strict=True)]))
    ycurve = {w: pl.ybar(CrossSlottedCircularPatch(6.0e-3, L_CROSS * 1e-3, w * 1e-3), 1.06) for w in W_GRID}
    print("local mean (30-50 mm) %.2f mS; cross Ybar(W):" % ylocal, {k: round(v, 2) for k, v in ycurve.items()}, flush=True)
    ws = np.array(W_GRID)
    ys = np.array([ycurve[w] for w in W_GRID])
    order = np.argsort(ys)
    y_high = float(ys.max())
    y_end = float(np.clip(ylocal, ys.min(), ys.max()))

    def w_of_y(v):
        return float(np.interp(v, ys[order], ws[order]))

    # pad lattice: the same 6.4 mm lattice, all cells with r < 30 mm except the feed cell
    n = 12
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(cp.PITCH * 1e3 * k, cp.PITCH * 1e3 * k)
    xp, yp = xg.ravel(), yg.ravel()
    rp = np.hypot(xp, yp)
    keep = (rp >= 3.2) & (rp < R_PAD)
    xp, yp, rp = xp[keep], yp[keep], rp[keep]
    yb = np.where(rp <= R_HIGH, y_high, y_high + (y_end - y_high) * (rp - R_HIGH) / (R_PAD - R_HIGH))
    wp = np.array([round(w_of_y(v) / 0.05) * 0.05 for v in yb])
    pad = [{"x_mm": float(a), "y_mm": float(b), "slot_L_mm": L_CROSS, "slot_W_mm": float(w), "slot_angle_deg": 0.0, "cross": 1}
           for a, b, w in zip(xp, yp, wp, strict=True)]
    keep_h = r >= R_PAD
    holo = [{"x_mm": float(a), "y_mm": float(b), "slot_L_mm": float(l), "slot_W_mm": float(w), "slot_angle_deg": float(an), "cross": 0}
            for a, b, l, w, an in zip(x[keep_h], y[keep_h], L[keep_h], W[keep_h], A[keep_h], strict=True)]
    allc = pad + holo
    Path(out).with_name(Path(out).stem + "_cells.json").write_text(json.dumps(allc))
    print("pad cells", len(pad), "hologram cells (unchanged)", len(holo), "removed hologram cells", int((~keep_h).sum()), flush=True)

    # ---- predictions -------------------------------------------------------------------------------
    table = SupercellTable.from_files(cp.TABLE_GLOB)
    poly = cp.hexagon()
    m = march_rays(table, poly, (0.0, 0.0), cp.K0, cp.THETA0, cp.PHI0, 0.9, taper_c=0.5, n_psi=720)
    idx = int(round(R_PAD * 1e-3 / m["drho"]))
    w30 = m["w"][:, idx]
    eta_old = float(np.mean(1.0 - m["w_edge"]))
    eta_new = float(np.mean(1.0 - m["w_edge"] / np.maximum(w30, 1e-12)))
    xs, ys_ = np.array([c["x_mm"] for c in holo]) * 1e-3, np.array([c["y_mm"] for c in holo]) * 1e-3
    d = cell_design(m, table, xs, ys_, (0.0, 0.0), cp.K0, cp.THETA0, cp.PHI0, hand=+1)
    # unchanged cells: same slots/angles as the stored design (consistency check), currents of r >= 30 mm only
    same = bool(np.allclose(np.degrees(d["slot_angle_rad"]), [c["slot_angle_deg"] for c in holo], atol=1e-6))
    green = cp.LayeredGreenFunction(cp.SUB, cp.F, backing="pec")
    theta, phi = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 2.0)
    et, ep = cp.pattern(d, xs, ys_, green, theta, phi)
    mt = cp.metrics(et, ep, theta, phi)
    # original (full aperture) for the ratio
    ca = SupercellTable.from_files(cp.TABLE_GLOB)
    xo, yo = cp.lattice_cells(poly)
    do = cell_design(m, ca, xo, yo, (0.0, 0.0), cp.K0, cp.THETA0, cp.PHI0, hand=+1)
    eo, po = cp.pattern(do, xo, yo, green, theta, phi)
    mo = cp.metrics(eo, po, theta, phi)
    # launch bracket from the calibrated feed model (top 6.0 mm above the surface, effective radius 1.5 mm)
    cal = lr.calibrate()
    z0, lf0 = lr.zin(lr.F0, lr.TOP_BASE, 0.0, 0.0, r_pin=lr.R_PIN)
    z1, lf1 = lr.zin(lr.F0, lr.TOP_BASE, 0.0, y_high, r_pin=lr.R_PIN)
    zl, lfl = lr.zin(lr.F0, lr.TOP_BASE, 0.0, 0.5 * y_high, r_pin=lr.R_PIN)
    c = cal["r_scale"]
    s11_now, s11_v2 = lr.s11_db(z0, c), lr.s11_db(z1, c)
    hfss = {"gain_dbi": 23.03, "main_beam_8deg": 0.271, "grazing_plus_back": 0.264 + 0.269, "lhcp_total": 0.310,
            "launch_current_model": lf0, "eta_old": eta_old}
    rows = {}
    for tag, lf in (("lower_launch_unchanged", lf0), ("central_half_of_the_infinite_sheet_gain", 0.5 * (lf0 + lf1)), ("upper_infinite_sheet", lf1)):
        dg = 10 * np.log10(lf / lf0) + 10 * np.log10(eta_new / eta_old) + (mt["directivity_ap_rhcp_dbi"] - mo["directivity_ap_rhcp_dbi"])
        main = hfss["main_beam_8deg"] * (lf * eta_new) / (lf0 * eta_old) * (mt["within_8deg_fraction"] / mo["within_8deg_fraction"])
        non_main = 1.0 - main
        base_non_main = 1.0 - hfss["main_beam_8deg"]
        rows[tag] = {
            "launch_fraction": float(lf), "gain_dbi": float(hfss["gain_dbi"] + dg), "delta_gain_db": float(dg),
            "main_beam_8deg_share": float(main),
            "grazing_plus_back_share": float(non_main * hfss["grazing_plus_back"] / base_non_main),
            "lhcp_total_share": float(non_main * hfss["lhcp_total"] / base_non_main),
        }
    res = {
        "labels": ["PREDICTION_NOT_HFSS_VALIDATED", "MODEL_ONLY", "HFSS_ANCHORED_TO_CURRENT_RUN"],
        "pad": {"R_high_mm": R_HIGH, "R_pad_mm": R_PAD, "L_cross_mm": L_CROSS, "ybar_high_mS": y_high, "ybar_hologram_local_mean_mS": ylocal,
                "ybar_cross_vs_W_mS_kt1.06": {str(k_): v for k_, v in ycurve.items()}, "W_at_pin_mm": float(w_of_y(y_high)),
                "W_at_edge_mm": float(w_of_y(y_end)), "n_pad_cells": len(pad), "n_hologram_cells_unchanged": len(holo),
                "n_hologram_cells_removed": int((~keep_h).sum()), "unslotted_D6_ybar_mS": 4.36},
        "hologram_cells_identical_to_current_design": same,
        "conversion": {"current": eta_old, "with_pad_no_radiation_inside_30mm": eta_new},
        "aperture_pattern_with_hole": {k: mt[k] for k in ("directivity_ap_rhcp_dbi", "xpd_at_peak_db", "ar_at_peak_db", "ar_max_over_minus3db_beam_db", "rhcp_sll_db", "within_8deg_fraction", "lhcp_fraction")},
        "aperture_pattern_current": {k: mo[k] for k in ("directivity_ap_rhcp_dbi", "xpd_at_peak_db", "ar_at_peak_db", "rhcp_sll_db", "within_8deg_fraction")},
        "feed_model": {"top_above_surface_mm": lr.TOP_BASE * 1e3, "r_pin_effective_mm": lr.R_PIN * 1e3, "calibration": cal,
                       "s11_db_current_est": float(s11_now), "s11_db_v2_est": float(s11_v2), "s11_db_half_loaded_est": float(lr.s11_db(zl, c)),
                       "z_current_ohm": [z0.real, z0.imag], "z_v2_ohm": [z1.real, z1.imag]},
        "current_run_hfss": hfss,
        "predictions": rows,
    }
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("conversion", "predictions")}, indent=1, default=float))
    print("S11 est now %.1f (HFSS -15.1), V2 %.1f" % (s11_now, s11_v2))


if __name__ == "__main__":
    main(sys.argv[1])
