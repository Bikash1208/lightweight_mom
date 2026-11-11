"""Analyse Option 1 low-Q eigenmode branches and line effects."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import cell_param_pyaedt as cp
from opt1_cell import configure
from opt1_geometry import all_geometries
from opt1_queue import cases

configure()
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "sub6_cell_param" / "lab_B_opt1"
TARGET_F = 5.8


def eligible(mode: dict) -> bool:
    return 60 < mode["q"] < 1000


def branches(rows: list[dict]) -> dict:
    by_phase = {r["phase_deg"]: r for r in rows}
    if set(by_phase) != {80, 90, 100}:
        raise RuntimeError("Missing phase rows")
    out = []
    for middle in by_phase[90]["modes"]:
        if not eligible(middle):
            continue
        low, amb_low = cp._track(middle, [m for m in by_phase[80]["modes"] if eligible(m)])
        high, amb_high = cp._track(middle, [m for m in by_phase[100]["modes"] if eligible(m)])
        item = {"mode90": middle["mode"], "f90_ghz": middle["f_ghz"], "q90": middle["q"],
                "distance_to_5p8_ghz": abs(middle["f_ghz"] - TARGET_F),
                "ambiguous_track": amb_low or amb_high}
        if low is None or high is None:
            item.update(classification="UNTRACKED", beta_over_k0=None, beta_status="UNTRACKED")
        else:
            slope = (high["f_ghz"] - low["f_ghz"]) / 20
            klass = "RES" if abs(slope) < 0.005 else "SW" if abs(slope) > 0.02 else "HYBRID"
            beta, beta_status = cp.beta_over_k0([80, 90, 100],
                                                 [low["f_ghz"], middle["f_ghz"], high["f_ghz"]], 11.6)
            item.update(f80_ghz=low["f_ghz"], q80=low["q"], f100_ghz=high["f_ghz"], q100=high["q"],
                        slope_ghz_per_deg=slope, vg_over_c=slope / 0.0718,
                        classification=klass, beta_over_k0=beta,
                        beta_status=beta_status)
            nearest_phase, nearest_f = min(((80, low["f_ghz"]), (90, middle["f_ghz"]),
                                            (100, high["f_ghz"])), key=lambda pair: abs(pair[1] - TARGET_F))
            item["distance_to_5p8_ghz"] = abs(nearest_f - TARGET_F)
            item["nearest_phase_deg"] = nearest_phase
        out.append(item)
    return {"modes": out, "eligible_mode90_count": len(out),
            "excluded_mode90_count": sum(not eligible(m) for m in by_phase[90]["modes"])}


def analyse_geometry(item: dict, all_rows: dict) -> dict:
    name = item["name"]
    geom_rows = [all_rows[r["*"]] for r in cases_from_item(item)]
    groups = {}
    raw_mode1_90 = {}
    for state in ("OFF", "ON"):
        for direction in ("along", "across"):
            group_rows = [r for r in geom_rows if r["state"] == state and r["direction"] == direction]
            groups[f"{state}_{direction}"] = branches(group_rows)
            middle = next(r for r in group_rows if r["phase_deg"] == 90)
            raw_mode1_90[f"{state}_{direction}"] = {
                "f_ghz": middle["modes"][0]["f_ghz"], "q": middle["modes"][0]["q"],
                "eligible": eligible(middle["modes"][0])}
    distances = {}
    beta = {}
    anisotropy = {}
    for state in ("OFF", "ON"):
        candidates = [(direction, mode) for direction in ("along", "across")
                      for mode in groups[f"{state}_{direction}"]["modes"]
                      if mode["classification"] in ("RES", "HYBRID")]
        if candidates:
            direction, nearest = min(candidates, key=lambda x: x[1]["distance_to_5p8_ghz"])
            distance = nearest["distance_to_5p8_ghz"]
            distances[state] = {"nearest_ghz": distance, "status": "PASS" if distance >= 1.0 else "FAIL",
                                "direction": direction, "mode90": nearest["mode90"],
                                "nearest_phase_deg": nearest["nearest_phase_deg"],
                                "classification": nearest["classification"]}
        else:
            distances[state] = {"nearest_ghz": None, "status": "NO_RES_OR_HYBRID_MODE"}
        beta[state] = {}
        for direction in ("along", "across"):
            tracked = [mode for mode in groups[f"{state}_{direction}"]["modes"]
                       if mode["classification"] in ("RES", "HYBRID", "SW") and mode["beta_over_k0"] is not None]
            if tracked:
                chosen = min(tracked, key=lambda x: x["distance_to_5p8_ghz"])
                beta[state][direction] = {"beta_over_k0": chosen["beta_over_k0"],
                                          "status": chosen["beta_status"],
                                          "classification": chosen["classification"],
                                          "mode90": chosen["mode90"],
                                          "ambiguous_track": chosen["ambiguous_track"]}
            else:
                beta[state][direction] = {"beta_over_k0": None, "status": "NO_TRACKED_LOW_Q_BRANCH"}
        along = beta[state]["along"]["beta_over_k0"]
        across = beta[state]["across"]["beta_over_k0"]
        anisotropy[state] = abs(across - along) if along is not None and across is not None else None
    ratio = (anisotropy["ON"] / anisotropy["OFF"]
             if anisotropy["ON"] is not None and anisotropy["OFF"] not in (None, 0) else None)
    return {"geometry": name, "route": item["geometry"]["route"],
            "island_l_mm": item["geometry"]["island_l"],
            "line_w_mm": item["geometry"]["line_w"], "channel_w_mm": item["geometry"]["channel_w"],
            "groups": groups, "raw_mode1_90": raw_mode1_90,
            "res_hybrid_distance": distances, "beta_at_5p8": beta,
            "anisotropy": anisotropy, "on_over_off_anisotropy_ratio": ratio}


def cases_from_item(item: dict) -> list[dict]:
    from opt1_geometry import Geometry
    return cases(Geometry(**item["geometry"]))


def compare_to_ref(entry: dict, ref: dict) -> dict:
    out = {"reference_geometry": ref["geometry"], "line_and_channel_both_changed_for_w0p20": entry["line_w_mm"] == 0.2}
    out["distance_change_ghz"] = {}
    out["anisotropy_change"] = {}
    out["beta_change_by_direction"] = {}
    out["nearest_mode90_shift_ghz"] = {}
    out["raw_mode1_90_shift_ghz"] = {}
    out["raw_mode1_90_comparison_excluded"] = {}
    for state in ("OFF", "ON"):
        distance = entry["res_hybrid_distance"][state]["nearest_ghz"]
        reference = ref["res_hybrid_distance"][state]["nearest_ghz"]
        out["distance_change_ghz"][state] = distance - reference if distance is not None and reference is not None else None
        aniso = entry["anisotropy"][state]
        ref_aniso = ref["anisotropy"][state]
        out["anisotropy_change"][state] = aniso - ref_aniso if aniso is not None and ref_aniso is not None else None
        out["beta_change_by_direction"][state] = {}
        out["nearest_mode90_shift_ghz"][state] = {}
        out["raw_mode1_90_shift_ghz"][state] = {}
        out["raw_mode1_90_comparison_excluded"][state] = {}
        for direction in ("along", "across"):
            beta = entry["beta_at_5p8"][state][direction]["beta_over_k0"]
            ref_beta = ref["beta_at_5p8"][state][direction]["beta_over_k0"]
            out["beta_change_by_direction"][state][direction] = beta - ref_beta if beta is not None and ref_beta is not None else None
            raw = entry["raw_mode1_90"][f"{state}_{direction}"]
            raw_ref = ref["raw_mode1_90"][f"{state}_{direction}"]
            out["raw_mode1_90_shift_ghz"][state][direction] = raw["f_ghz"] - raw_ref["f_ghz"]
            out["raw_mode1_90_comparison_excluded"][state][direction] = not (raw["eligible"] and raw_ref["eligible"])
            modes = entry["groups"][f"{state}_{direction}"]["modes"]
            ref_modes = ref["groups"][f"{state}_{direction}"]["modes"]
            if modes and ref_modes:
                nearest = min(modes, key=lambda x: x["distance_to_5p8_ghz"])
                nearest_ref = min(ref_modes, key=lambda x: abs(x["f90_ghz"] - nearest["f90_ghz"]))
                out["nearest_mode90_shift_ghz"][state][direction] = nearest["f90_ghz"] - nearest_ref["f90_ghz"]
            else:
                out["nearest_mode90_shift_ghz"][state][direction] = None
    return out


def main() -> None:
    manifest = json.loads((RESULTS / "row_manifest.json").read_text())
    rows = {p.stem: json.loads(p.read_text()) for p in (RESULTS / "rows").glob("*.json")}
    if set(rows) != set(manifest):
        raise RuntimeError(f"Need all 120 rows for analysis: {len(rows)} present")
    items = all_geometries(0.5)
    results = [analyse_geometry(item, rows) for item in items]
    by_name = {x["geometry"]: x for x in results}
    for entry in results:
        if entry["route"] != "REF":
            ref_name = f"REF_L{entry['island_l_mm']:.3f}"
            entry["vs_ref"] = compare_to_ref(entry, by_name[ref_name])
    out = {"rule": {"q_excluded": "Q <= 60 or Q >= 1000", "res_slope_abs_lt": 0.005,
                    "sw_slope_abs_gt": 0.02, "hybrid": "between thresholds",
                    "vg_over_c": "signed slope / 0.0718", "target_ghz": 5.8,
                    "res_hybrid_distance_pass_ghz": 1.0,
                    "beta_branch": "nearest tracked low-Q branch at 5.8 GHz, with RES/HYBRID/SW class; bracketed or EXTRAPOLATED"},
           "geometries": results}
    (RESULTS / "analysis.json").write_text(json.dumps(out, indent=2) + "\n")
    with (RESULTS / "analysis_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["geometry", "state", "nearest_RES_HYBRID_GHz", "distance_status",
                         "beta_along", "beta_along_status", "beta_across", "beta_across_status",
                         "anisotropy", "ON_over_OFF_ratio", "delta_distance_vs_REF_GHz", "delta_anisotropy_vs_REF"])
        for entry in results:
            for state in ("OFF", "ON"):
                dist = entry["res_hybrid_distance"][state]
                comparison = entry.get("vs_ref", {})
                writer.writerow([entry["geometry"], state, dist["nearest_ghz"], dist["status"],
                                 entry["beta_at_5p8"][state]["along"]["beta_over_k0"],
                                 entry["beta_at_5p8"][state]["along"]["status"],
                                 entry["beta_at_5p8"][state]["across"]["beta_over_k0"],
                                 entry["beta_at_5p8"][state]["across"]["status"],
                                 entry["anisotropy"][state], entry["on_over_off_anisotropy_ratio"],
                                 comparison.get("distance_change_ghz", {}).get(state),
                                 comparison.get("anisotropy_change", {}).get(state)])
    with (RESULTS / "line_vs_ref.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["geometry", "reference", "state", "direction", "raw_mode1_90_shift_GHz",
                         "raw_comparison_excluded_by_Q_rule", "nearest_lowQ_mode90_shift_GHz",
                         "beta_change", "distance_change_GHz"])
        for entry in results:
            if entry["route"] == "REF":
                continue
            comparison = entry["vs_ref"]
            for state in ("OFF", "ON"):
                for direction in ("along", "across"):
                    writer.writerow([entry["geometry"], comparison["reference_geometry"], state, direction,
                                     comparison["raw_mode1_90_shift_ghz"][state][direction],
                                     comparison["raw_mode1_90_comparison_excluded"][state][direction],
                                     comparison["nearest_mode90_shift_ghz"][state][direction],
                                     comparison["beta_change_by_direction"][state][direction],
                                     comparison["distance_change_ghz"][state]])
    fails = sum(entry["res_hybrid_distance"][state]["status"] == "FAIL"
                for entry in results for state in ("OFF", "ON"))
    missing = sum(entry["res_hybrid_distance"][state]["status"] == "NO_RES_OR_HYBRID_MODE"
                  for entry in results for state in ("OFF", "ON"))
    print(f"ANALYSIS: 10 geometries, {fails} distance failures, {missing} states with no RES/HYBRID mode")


if __name__ == "__main__":
    main()
