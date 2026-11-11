"""Report every seven-phase Option 4 mode and lid-screened 5.8 GHz metrics."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from opt1b_analyze import make_tracks, match_lids
from opt4_geometry import GEOMETRIES
from opt4_spec import LIDS, PHASES, cases, variant

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1] / "results" / "sub6_cell_param" / "lab_B_opt4"
F0 = 5.8
LID_SHIFT_PCT = 0.5


def load_tracks(g, lid, state, direction):
    selected = [r for r in cases(g, lid) if r["state"] == state and r["direction"] == direction]
    files = [ROOT / "rows" / f"{r['*']}.json" for r in selected]
    if any(not p.exists() for p in files):
        return None
    rows = {int(r["phase_deg"]): json.loads(file.read_text()) for r, file in zip(selected, files)}
    return make_tracks(rows)


def main():
    failures = json.loads((ROOT / "FAILED_ROWS.json").read_text()) if (ROOT / "FAILED_ROWS.json").exists() else []
    if failures:
        raise RuntimeError(f"FAILED rows present: {len(failures)}")
    tracks = {}
    for g in GEOMETRIES:
        for lid in LIDS:
            for state in ("OFF", "ON"):
                for direction in ("along", "across"):
                    group = load_tracks(g, lid, state, direction)
                    if group is not None:
                        tracks[(g.name, lid, state, direction)] = group
    for g in GEOMETRIES:
        for state in ("OFF", "ON"):
            for direction in ("along", "across"):
                base = tracks.get((g.name, "PEC150", state, direction))
                p300 = tracks.get((g.name, "PEC300", state, direction))
                z377 = tracks.get((g.name, "Z377_150", state, direction))
                if base is None:
                    continue
                pmatches = match_lids(base, p300) if p300 else None
                zmatches = match_lids(base, z377) if z377 else None
                for i, t in enumerate(base):
                    t["PEC300_match"] = pmatches[i] if pmatches else None
                    t["Z377_match"] = zmatches[i] if zmatches else None
                    t["box_mode"] = (pmatches[i]["max_frequency_shift_pct"] > LID_SHIFT_PCT
                                     if pmatches else None)
    out = ROOT / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    fields = ("geometry", "lid", "state", "direction", "track", "class", "slope_ghz_per_deg",
              "vg_over_c", "q_min", "q_median", "q_max", "beta_brackets", "PEC300_match",
              "Z377_match", "box_mode", "absorber_only_Q_le_60", "phase_deg", "mode", "f_ghz", "q")
    with (out / "every_mode_track.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for (geometry, lid, state, direction), group in tracks.items():
            for t in group:
                for p in t["points"]:
                    writer.writerow({"geometry": geometry, "lid": lid, "state": state,
                                     "direction": direction, "track": t["track"], "class": t["class"],
                                     "slope_ghz_per_deg": t["slope_ghz_per_deg"],
                                     "vg_over_c": t["vg_over_c"], "q_min": t["q_min"],
                                     "q_median": t["q_median"], "q_max": t["q_max"],
                                     "beta_brackets": json.dumps(t["beta_brackets"]),
                                     "PEC300_match": json.dumps(t.get("PEC300_match")),
                                     "Z377_match": json.dumps(t.get("Z377_match")),
                                     "box_mode": t.get("box_mode"),
                                     "absorber_only_Q_le_60": lid == "Z377_150" and t["q_max"] <= 60,
                                     "phase_deg": p["phase_deg"], "mode": p["mode"],
                                     "f_ghz": p["f_ghz"], "q": p["q"]})
    summary = {"rows_available": len(list((ROOT / "rows").glob("OPT4_*.json"))),
               "q_ge_1000_excluded": False, "absorber_only_q_limit": 60,
               "box_mode_pec_height_shift_limit_pct": LID_SHIFT_PCT,
               "geometries": {}}
    for g in GEOMETRIES:
        entry = {"states": {}}
        for state in ("OFF", "ON"):
            state_entry = {"axes": {}}
            flat = []
            for direction in ("along", "across"):
                group = tracks.get((g.name, "PEC150", state, direction))
                if group is None:
                    state_entry["axes"][direction] = {"status": "NOT_SOLVED"}
                    continue
                screened = group[0]["box_mode"] is not None
                candidates = [{"track": t["track"], "class": t["class"],
                               "q_median": t["q_median"], "box_mode": t["box_mode"], **b}
                              for t in group if t["box_mode"] is not True and t["class"] in ("SW", "HYBRID")
                              for b in t["beta_brackets"]]
                candidates.sort(key=lambda x: x["beta_over_k0"])
                selected = candidates[0] if screened and len(candidates) == 1 else None
                status = ("UNSCREENED" if not screened else "BRACKETED" if selected else
                          "NO_BRACKET" if not candidates else "AMBIGUOUS")
                for t in group:
                    if t["box_mode"] is not True and t["class"] in ("RES", "HYBRID"):
                        nearest = min(t["points"], key=lambda p: abs(p["f_ghz"] - F0))
                        flat.append({"direction": direction, "track": t["track"],
                                     "class": t["class"], "f_ghz": nearest["f_ghz"],
                                     "q": nearest["q"], "distance_ghz": abs(nearest["f_ghz"] - F0),
                                     "lid_screened": screened})
                state_entry["axes"][direction] = {"status": status,
                                                  "bracketed_candidates": candidates,
                                                  "selected": selected,
                                                  "lid_screened": screened,
                                                  "box_track_count": sum(t["box_mode"] is True for t in group)}
            flat.sort(key=lambda x: x["distance_ghz"])
            nearest = flat[0] if flat else None
            state_entry["nearest_RES_or_HYBRID"] = nearest
            state_entry["resonance_target_pass"] = (nearest["distance_ghz"] >= 1.0
                                                    if nearest and nearest["lid_screened"] else None)
            along = state_entry["axes"].get("along", {}).get("selected")
            across = state_entry["axes"].get("across", {}).get("selected")
            state_entry["anisotropy_across_minus_along"] = (
                across["beta_over_k0"] - along["beta_over_k0"] if along and across else None)
            entry["states"][state] = state_entry
        off = entry["states"]["OFF"]["anisotropy_across_minus_along"]
        on = entry["states"]["ON"]["anisotropy_across_minus_along"]
        entry["anisotropy_ON_over_OFF"] = on / off if on is not None and off not in (None, 0) else None
        entry["targets_pass"] = (off >= 0.05 and on <= 0.2 * off and
                                 all(entry["states"][s]["resonance_target_pass"] for s in ("OFF", "ON"))
                                 if off is not None and on is not None and
                                 all(entry["states"][s]["resonance_target_pass"] is not None
                                     for s in ("OFF", "ON")) else None)
        summary["geometries"][g.name] = entry
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
