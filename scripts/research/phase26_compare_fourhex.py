"""Phase 26: compare the user's HFSS run of board A (`hfss_exports/hex4_tensor_cp_static`) with `predictions.json`.

    python scripts/research/phase26_compare_fourhex.py --s4p board.s4p
        --gain 3=gain_port3.csv --gain 2=gain_port2.csv [--params 3=params3.txt --params 2=params2.txt]
        [--pred predictions.json] [--out report.json]

No modelling: every number is read from the exports or `predictions.json`. Per driven port with a gain CSV (HFSS
realized-gain RHCP/LHCP, 1 deg): peak RHCP gain, pointing, XPD, AR at peak / worst over the -3 dB beam, RHCP SLL, the
strongest LHCP lobe outside 9 deg of the beam (level re the RHCP peak, direction) and the LHCP level read at the
PREDICTED neighbour-lobe direction (max within 6 deg), against the table of `predictions.json` (s = 1 and s = 1.7).
From the `.s4p` at 10.54 GHz: S_ii and the five adjacent-pair couplings against the raw/calibrated estimates.
Named states: GAIN_WITHIN_PREDICTED_RANGE / GAIN_OUTSIDE_PREDICTED_RANGE, LHCP_LOBE_AT_OR_BELOW_PREDICTION /
LHCP_LOBE_ABOVE_PREDICTION, COUPLING_BELOW_RAW_ESTIMATE / COUPLING_ABOVE_RAW_ESTIMATE, NOT_PROVIDED, NO_PREDICTION.
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
PRED = ROOT / "hfss_exports/hex4_tensor_cp_static/predictions.json"
PHI = {1: 45.0, 2: 135.0, 3: 225.0, 4: 315.0}
GAIN_TOL_DB = 1.0
F_DESIGN = 10.54
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]

_spec = importlib.util.spec_from_file_location("p26_single", Path(__file__).resolve().parent / "phase26_compare_single_hex.py")
single = importlib.util.module_from_spec(_spec)
sys.modules["p26_single"] = single
_spec.loader.exec_module(single)


def _kv(items):
    out = {}
    for it in items or []:
        k, v = it.split("=", 1)
        out[int(k)] = v
    return out


def port_report(port, gain_path, params_path, pred):
    single.BEAM = (30.0, PHI[port])
    g = single.read_gain_csv(gain_path, F_DESIGN)
    a, (tg, pg, _ar) = single.analyse(g)
    th, ph, pl, pr = g["theta"], g["phi"], g["pl"], g["pr"]
    tgrid, pgrid = np.meshgrid(th, ph)
    peak = float(np.nanmax(pr))
    d_beam = single._angle_to(tgrid, pgrid, *single.BEAM)
    ok = np.isfinite(pl) & (d_beam > 9.0) & (tgrid <= 90.0)
    k = np.unravel_index(np.argmax(np.where(ok, pl, -1.0)), pl.shape)
    lobe = {"level_rel_rhcp_peak_db": float(10 * np.log10(pl[k] / peak)), "theta_phi_deg": [float(tgrid[k]), float(pgrid[k])], "deg_from_beam": float(d_beam[k])}
    rep = {"measured": a, "strongest_lhcp_lobe_outside_9deg": lobe, "freq_ghz": g["freq_ghz"]}
    if params_path:
        rep["antenna_params"] = single.read_params(params_path)
    p = pred.get("per_port_alone_by_edge_power_scale")
    if not p:
        rep["checks"] = {"state": "NO_PREDICTION"}
        return rep
    p1, p17 = p["s=1.0"][str(port)], p["s=1.7"][str(port)]
    lo, hi = min(p1["realized_gain_rhcp_dbi_anchored"], p17["realized_gain_rhcp_dbi_anchored"]) - GAIN_TOL_DB, max(p1["realized_gain_rhcp_dbi_anchored"], p17["realized_gain_rhcp_dbi_anchored"]) + GAIN_TOL_DB
    gdb = a["peak_rhcp_gain_dbi"]
    pl_theta_phi = p1["check_item_neighbour_lhcp_lobe"]["theta_phi_deg"]
    near = single._angle_to(tgrid, pgrid, *pl_theta_phi) <= 6.0
    at_pred = float(10 * np.log10(np.nanmax(np.where(near & np.isfinite(pl), pl, 0.0)) / peak))
    pred_lvl = (p1["check_item_neighbour_lhcp_lobe"]["rel_to_main_rhcp_peak_db"], p17["check_item_neighbour_lhcp_lobe"]["rel_to_main_rhcp_peak_db"])
    rep["checks"] = {
        "peak_gain": {"hfss_dbi": gdb, "predicted_range_dbi": [lo, hi], "state": "GAIN_WITHIN_PREDICTED_RANGE" if lo <= gdb <= hi else "GAIN_OUTSIDE_PREDICTED_RANGE"},
        "pointing_error_deg": a["pointing_error_deg"],
        "xpd_db": {"hfss": a["xpd_at_peak_db"], "model_s1": p1["xpd_at_peak_db"]},
        "ar_peak_db": {"hfss": a["ar_at_peak_db"], "model_s1": p1["ar_at_peak_db"], "worst_in_beam_hfss": a["ar_max_over_minus3db_beam_db"],
                       "worst_in_beam_model_s1": p1["ar_worst_in_minus3db_beam_db"]},
        "neighbour_lhcp_lobe": {"hfss_at_predicted_direction_db": at_pred, "predicted_direction_theta_phi_deg": pl_theta_phi, "predicted_db_s1_s1p7": list(pred_lvl),
                                "hfss_strongest_lhcp_outside_9deg": lobe,
                                "state": "LHCP_LOBE_AT_OR_BELOW_PREDICTION" if at_pred <= max(pred_lvl) + 1.0 else "LHCP_LOBE_ABOVE_PREDICTION"},
    }
    return rep


def coupling_report(s4p, pred):
    from src.phase9.hfss_reference import read_touchstone

    net = read_touchstone(s4p)
    f = net.frequency_hz / 1e9
    i = int(np.argmin(np.abs(f - F_DESIGN)))
    db = lambda a, b: float(20 * np.log10(max(abs(net.component(a, b)[i]), 1e-12)))  # noqa: E731
    out = {"f_nearest_ghz": float(f[i]), "s_ii_db": {str(p): db(p, p) for p in range(1, 5)}, "pairs": {}}
    est = pred.get("coupling", {}).get("adjacent_pairs", {})
    for a, b in PAIRS:
        v = {"s_ij_db": db(a, b), "s_ji_db": db(b, a), "worst_over_sweep_db": float(np.max(20 * np.log10(np.maximum(np.abs(net.component(a, b)), 1e-12))))}
        e = est.get(f"{a}-{b}")
        if e:
            v["raw_estimate_db"], v["calibrated_estimate_db"] = e["raw_estimate_db"], e.get("calibrated_estimate_db")
            v["state"] = "COUPLING_BELOW_RAW_ESTIMATE" if v["s_ij_db"] <= e["raw_estimate_db"] else "COUPLING_ABOVE_RAW_ESTIMATE"
        else:
            v["state"] = "NO_PREDICTION"
        out["pairs"][f"{a}-{b}"] = v
    out["not_adjacent_2-4_db"] = db(2, 4)
    return out


def compare(args):
    pred = json.loads(Path(args.pred).read_text())
    rep = {"pred": str(args.pred), "ports": {}, "coupling": "NOT_PROVIDED"}
    gains, params = _kv(args.gain), _kv(args.params)
    for port, path in gains.items():
        rep["ports"][str(port)] = port_report(port, path, params.get(port), pred)
    if args.s4p:
        rep["coupling"] = coupling_report(args.s4p, pred)
    return rep


def markdown(rep):
    L = ["# hex4_tensor_cp_static: HFSS vs prediction", ""]
    for p, r in rep["ports"].items():
        c = r["checks"]
        if "peak_gain" not in c:
            L.append(f"* port {p}: {c['state']}")
            continue
        a = r["measured"]
        L += [f"## Port {p} alone (beam 30, {PHI[int(p)]:g})", "",
              f"* peak RHCP {c['peak_gain']['hfss_dbi']:.2f} dBi vs predicted {c['peak_gain']['predicted_range_dbi'][0]:.1f}..{c['peak_gain']['predicted_range_dbi'][1]:.1f}: **{c['peak_gain']['state']}**; "
              f"pointing error {a['pointing_error_deg']:.1f} deg",
              f"* XPD {c['xpd_db']['hfss']:.1f} dB (model {c['xpd_db']['model_s1']:.1f}); AR {c['ar_peak_db']['hfss']:.2f} dB at peak / {c['ar_peak_db']['worst_in_beam_hfss']:.2f} worst in beam "
              f"(model {c['ar_peak_db']['model_s1']:.2f} / {c['ar_peak_db']['worst_in_beam_model_s1']:.2f}); RHCP SLL {a['rhcp_sll_db']:.1f} dB",
              f"* neighbour LHCP lobe: HFSS {c['neighbour_lhcp_lobe']['hfss_at_predicted_direction_db']:.1f} dB at the predicted direction "
              f"{c['neighbour_lhcp_lobe']['predicted_direction_theta_phi_deg']} (predicted {c['neighbour_lhcp_lobe']['predicted_db_s1_s1p7'][0]:.1f} / {c['neighbour_lhcp_lobe']['predicted_db_s1_s1p7'][1]:.1f} dB): "
              f"**{c['neighbour_lhcp_lobe']['state']}**; strongest LHCP outside 9 deg: {r['strongest_lhcp_lobe_outside_9deg']['level_rel_rhcp_peak_db']:.1f} dB at "
              f"{r['strongest_lhcp_lobe_outside_9deg']['theta_phi_deg']}",
              f"* budget: main beam <= 8 deg {a['budget_of_radiated']['main_beam_le_8deg_total']:.2f}, grazing {a['budget_of_radiated']['grazing_70_to_90']:.2f}, "
              f"back {a['budget_of_radiated']['back_hemisphere_gt_90']:.2f}, LHCP {a['budget_of_radiated']['lhcp_total_all_angles']:.2f}", ""]
    cp_ = rep["coupling"]
    if cp_ == "NOT_PROVIDED":
        L.append("Coupling: NOT_PROVIDED")
    else:
        L += ["## Coupling at %.2f GHz" % cp_["f_nearest_ghz"], "", "S_ii: " + ", ".join(f"{k}: {v:.1f} dB" for k, v in cp_["s_ii_db"].items()), "",
              "| pair | S_ij dB | worst over sweep | raw estimate | calibrated | state |", "|---|---|---|---|---|---|"]
        for k, v in cp_["pairs"].items():
            L.append(f"| {k} | {v['s_ij_db']:.1f} | {v['worst_over_sweep_db']:.1f} | {v.get('raw_estimate_db', float('nan')):.1f} | {v.get('calibrated_estimate_db') or float('nan'):.1f} | {v['state']} |")
        L.append(f"\n2-4 (not adjacent): {cp_['not_adjacent_2-4_db']:.1f} dB")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s4p")
    ap.add_argument("--gain", action="append")
    ap.add_argument("--params", action="append")
    ap.add_argument("--pred", default=str(PRED))
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    rep = compare(args)
    txt = markdown(rep)
    print(txt)
    if args.out:
        Path(args.out).write_text(json.dumps(rep, indent=1, default=float))
        Path(args.out).with_suffix(".md").write_text(txt)
    return rep


if __name__ == "__main__":
    main()
