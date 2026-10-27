"""Phase 28 item 2 at two frequencies: what changes between the 10.54 GHz sweep (history) and the 10.0 GHz sweep (the PI's operating frequency).

    python scripts/research/phase28_compare_frequencies.py <json_10p54> <json_10p0> <out_dir>

Writes `<out_dir>/frequency_comparison.json` (every compared quantity at both frequencies and the difference) and `<out_dir>/FREQUENCY_COMPARISON.md` generated from that
JSON only (`tests/test_phase28_compare_frequencies.py` checks that every number in the markdown is a value of the comparison JSON).
"""

import json
import sys
from pathlib import Path


def _case(d, pitch, theta, phi, r=150.0):
    for c in d["cases"]:
        if (c["pitch_mm"], c["aperture_r_mm"], c["theta0_deg"], c["phi0_deg"]) == (pitch, r, theta, phi):
            return c
    raise KeyError((pitch, theta, phi))


def _steer(d, theta):
    return next(s for s in d["contrast"]["steerings"] if s["theta0_deg"] == theta)


def quantities(d):
    dp = d["settings"]["design_point"]
    cs = _case(d, dp["pitch_mm"], dp["theta0_deg"], dp["phi0_deg"])
    st = _steer(d, dp["theta0_deg"])
    q = {
        "frequency_ghz": d["constants"]["frequency_ghz"],
        "table_alpha_max_np_per_m": d["contrast"]["table_alpha_max_np_per_m"],
        "table_alpha_min_np_per_m": d["contrast"].get(
            "table_alpha_min_core_columns_np_per_m", d["contrast"]["table_alpha_min_np_per_m"]
        ),
        "supported_scan_theta_deg": d["settings"]["supported_scan_theta_deg"],
        "design_point_eta_conv": cs["eta_conv"],
        "design_point_outside_table_fraction": cs["table_extrapolated_fraction"],
        "design_point_cp_continuous_gain_dbi": cs["variants"]["cp_cont"]["gain_dbi"],
        "design_point_cp_1bit_loss_db": cs["variants"]["cp_1bit"]["quantisation_loss_db"],
        "design_point_cp_2bit_loss_db": cs["variants"]["cp_2bit"]["quantisation_loss_db"],
        "design_point_cp_3bit_loss_db": cs["variants"]["cp_3bit"]["quantisation_loss_db"],
        "alpha_cap_for_1db_np_per_m": st["alpha_cap_for_1db_np_per_m"],
        "tensor_delta_ms": st["contrast"]["bits_2"]["tensor"]["delta_ms"],
        "tensor_dx_aniso_ohm": st["contrast"]["bits_2"]["tensor"]["dx_aniso_ohm"],
        "mean_beta_over_k0": st["contrast"]["bits_2"]["beta_over_k0_mean"],
        "mean_loading_ms": st["contrast"]["bits_2"]["y_mean_ms"],
    }
    for ce in st["combined_2bit"]:
        tag = f"combined_{ce['tolerance_db']:g}db"
        q[tag + "_alpha_cap_np_per_m"] = ce["alpha_cap_np_per_m"]
        q[tag + "_dx_aniso_ohm"] = None if ce["contrast"] is None else ce["contrast"]["tensor"]["dx_aniso_ohm"]
        q[tag + "_delta_ms"] = None if ce["contrast"] is None else ce["contrast"]["tensor"]["delta_ms"]
    for p in sorted({c["pitch_mm"] for c in d["cases"]}):
        q[f"cp_2bit_loss_db_pitch_{p:g}"] = _case(d, p, dp["theta0_deg"], dp["phi0_deg"])["variants"]["cp_2bit"][
            "quantisation_loss_db"
        ]
    for t in sorted(
        {
            c["theta0_deg"]
            for c in d["cases"]
            if c["pitch_mm"] == dp["pitch_mm"] and c["aperture_r_mm"] == dp["aperture_r_mm"]
        }
    ):
        phi = 0.0 if t == 0.0 else d["settings"]["scan_azimuth_deg"]
        q[f"outside_table_fraction_theta_{t:g}"] = _case(d, dp["pitch_mm"], t, phi)["table_extrapolated_fraction"]
        q[f"cp_2bit_gain_dbi_theta_{t:g}"] = _case(d, dp["pitch_mm"], t, phi)["variants"]["cp_2bit"][
            "gain_headline_dbi"
        ]
    return q


def compare(a, b):
    qa, qb = quantities(a), quantities(b)
    rows = []
    for k in qa:
        va, vb = qa[k], qb.get(k)
        diff = None if (va is None or vb is None or not isinstance(va, (int, float))) else vb - va
        rows.append({"quantity": k, "at_10p54": va, "at_10p0": vb, "difference": diff})
    return {"label": "PHASE28_FREQUENCY_COMPARISON", "rows": rows}


def markdown(c):
    lines = [
        "# Phase 28 item 2: 10.0 GHz against 10.54 GHz",
        "",
        "Generated from `frequency_comparison.json` only. Both sweeps use the same code (`src/hologram/coded_feasibility.py`) and the same cases; they differ in the frequency, the supercell table (rebuilt at the new frequency with the same method)",
        "and, as an ASSUMPTION carried over unchanged, the launch fraction (the 10.54 GHz Phase 26 value; it shifts every absolute gain equally and no loss). `difference` = 10.0 GHz value minus 10.54 GHz value.",
        "",
        "| quantity | 10.54 GHz | 10.0 GHz | difference |",
        "|---|---|---|---|",
    ]

    def fmt(v):
        if v is None:
            return "n/a"
        return f"{v:.4g}" if isinstance(v, float) else str(v)

    for r in c["rows"]:
        lines.append(f"| {r['quantity']} | {fmt(r['at_10p54'])} | {fmt(r['at_10p0'])} | {fmt(r['difference'])} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    a, b, out = (Path(x) for x in sys.argv[1:4])
    comp = compare(json.loads(a.read_text(encoding="utf-8")), json.loads(b.read_text(encoding="utf-8")))
    out.mkdir(parents=True, exist_ok=True)
    (out / "frequency_comparison.json").write_text(json.dumps(comp, indent=1), encoding="utf-8", newline="\n")
    (out / "FREQUENCY_COMPARISON.md").write_text(markdown(comp), encoding="utf-8", newline="\n")
    print("wrote", out)
