"""Track every solved Option 1b mode; classify box modes by lid movement only."""
from __future__ import annotations

import csv
from functools import lru_cache
import json
import math
from pathlib import Path
import statistics

from opt1b_spec import GEOMETRIES, LIDS, PHASES, cases, variant

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1] / "results" / "sub6_cell_param" / "lab_B_opt1b"
TARGET_GHZ = 5.8
PERIOD_MM = 11.6
LID_SHIFT_LIMIT_PCT = 0.5


def assignment(cost):
    """Minimum-cost one-to-one mapping, without a scipy dependency."""
    n = len(cost)

    @lru_cache(None)
    def solve(i, used):
        if i == n:
            return (0.0, ())
        options = []
        for j in range(n):
            if not used & (1 << j):
                rest, order = solve(i + 1, used | (1 << j))
                options.append((cost[i][j] + rest, (j,) + order))
        return min(options)

    return solve(0, 0)[1]


def slope(points):
    xs = PHASES
    ys = [p["f_ghz"] for p in points]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def beta_at_target(points):
    choices = []
    for a, b in zip(points, points[1:]):
        fa, fb = a["f_ghz"], b["f_ghz"]
        if min(fa, fb) <= TARGET_GHZ <= max(fa, fb) and fa != fb:
            phase = a["phase_deg"] + (b["phase_deg"] - a["phase_deg"]) * (TARGET_GHZ - fa) / (fb - fa)
            beta = phase * 299792458 / (360 * PERIOD_MM * 1e-3 * TARGET_GHZ * 1e9)
            choices.append({"beta_over_k0": beta, "phase_deg": phase,
                            "bracket_phase_deg": [a["phase_deg"], b["phase_deg"]],
                            "bracket_f_ghz": [fa, fb]})
    return choices


def make_tracks(rows):
    first = rows[PHASES[0]]["modes"]
    tracks = [[dict(m, phase_deg=PHASES[0])] for m in first]
    previous = first
    for phase in PHASES[1:]:
        current = rows[phase]["modes"]
        if len(current) != len(tracks):
            raise RuntimeError("Mode count changed across phase")
        cost = []
        for t in tracks:
            predicted = t[-1]["f_ghz"]
            if len(t) > 1:
                predicted += t[-1]["f_ghz"] - t[-2]["f_ghz"]
            cost.append([abs(m["f_ghz"] - predicted) / 0.1 +
                         0.08 * abs(math.log(m["q"] / t[-1]["q"])) for m in current])
        order = assignment(cost)
        for i, j in enumerate(order):
            tracks[i].append(dict(current[j], phase_deg=phase))
        previous = current
    result = []
    for i, points in enumerate(tracks, 1):
        s = slope(points)
        klass = "RES" if abs(s) < 0.005 else "SW" if abs(s) > 0.02 else "HYBRID"
        result.append({"track": i, "points": points, "slope_ghz_per_deg": s,
                       "vg_over_c": s / 0.0718, "class": klass,
                       "q_min": min(p["q"] for p in points),
                       "q_median": statistics.median(p["q"] for p in points),
                       "q_max": max(p["q"] for p in points),
                       "beta_brackets": beta_at_target(points)})
    return result


def match_lids(base, other):
    cost = [[statistics.mean(abs(a["f_ghz"] - b["f_ghz"]) / a["f_ghz"] for a, b in
                             zip(t["points"], u["points"])) +
             0.01 * abs(math.log(t["q_median"] / u["q_median"])) for u in other] for t in base]
    order = assignment(cost)
    return [{"track": order[i] + 1,
             "max_frequency_shift_pct": 100 * max(abs(a["f_ghz"] - b["f_ghz"]) / a["f_ghz"] for a, b in
                                                 zip(base[i]["points"], other[order[i]]["points"])),
             "mean_frequency_shift_pct": 100 * statistics.mean(abs(a["f_ghz"] - b["f_ghz"]) / a["f_ghz"] for a, b in
                                                               zip(base[i]["points"], other[order[i]]["points"]))}
            for i in range(len(base))]


def main():
    expected = [r["*"] for g in GEOMETRIES for lid, _, _ in LIDS for r in cases(g, lid)]
    missing = [x for x in expected if not (ROOT / "rows" / f"{x}.json").exists()]
    if missing:
        raise RuntimeError(f"Missing {len(missing)} rows; analysis requires complete sweep")
    failure = json.loads((ROOT / "FAILED_ROWS.json").read_text())
    if failure:
        raise RuntimeError(f"FAILED rows present: {len(failure)}")
    tracks = {}
    for g in GEOMETRIES:
        for lid, _, _ in LIDS:
            for state in ("OFF", "ON"):
                for direction in ("along", "across"):
                    selected = [r for r in cases(g, lid) if r["state"] == state and r["direction"] == direction]
                    rows = {int(r["phase_deg"]): json.loads((ROOT / "rows" / f"{r['*']}.json").read_text()) for r in selected}
                    key = (g.name, lid, state, direction)
                    tracks[key] = make_tracks(rows)
    # A mode is excluded only if its frequency moves by >0.5% between the two PEC heights.
    # The impedance comparison is reported separately because loss changes Q by construction.
    for g in GEOMETRIES:
        for state in ("OFF", "ON"):
            for direction in ("along", "across"):
                key = (g.name, "PEC150", state, direction)
                base = tracks[key]
                p300 = match_lids(base, tracks[(g.name, "PEC300", state, direction)])
                z377 = match_lids(base, tracks[(g.name, "Z377_150", state, direction)])
                for i, t in enumerate(base):
                    t["PEC300_match"] = p300[i]
                    t["Z377_match"] = z377[i]
                    t["box_mode"] = p300[i]["max_frequency_shift_pct"] > LID_SHIFT_LIMIT_PCT
                    t["lid_stable"] = not t["box_mode"]
    # Export all solved tracks and all points, including excluded box tracks and all Q values.
    fields = ["geometry", "lid", "state", "direction", "track", "class", "slope_ghz_per_deg",
              "vg_over_c", "q_min", "q_median", "q_max", "beta_brackets",
              "PEC300_match", "Z377_match", "box_mode", "phase_deg", "mode", "f_ghz", "q"]
    out = ROOT / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "every_mode_track.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for (geometry, lid, state, direction), group in tracks.items():
            for t in group:
                for p in t["points"]:
                    writer.writerow({"geometry": geometry, "lid": lid, "state": state,
                                     "direction": direction, "track": t["track"], "class": t["class"],
                                     "slope_ghz_per_deg": t["slope_ghz_per_deg"], "vg_over_c": t["vg_over_c"],
                                     "q_min": t["q_min"], "q_median": t["q_median"], "q_max": t["q_max"],
                                     "beta_brackets": json.dumps(t["beta_brackets"]),
                                     "PEC300_match": json.dumps(t.get("PEC300_match")),
                                     "Z377_match": json.dumps(t.get("Z377_match")),
                                     "box_mode": t.get("box_mode"), "phase_deg": p["phase_deg"],
                                     "mode": p["mode"], "f_ghz": p["f_ghz"], "q": p["q"]})
    summary = {"lid_shift_limit_pct": LID_SHIFT_LIMIT_PCT, "q_exclusion": None,
               "rows": len(expected), "by_geometry": {}}
    for g in GEOMETRIES:
        entry = {"states": {}}
        for state in ("OFF", "ON"):
            axes = {}
            flat = []
            for direction in ("along", "across"):
                group = tracks[(g.name, "PEC150", state, direction)]
                candidates = [{"track": t["track"], **b, "q_median": t["q_median"]} for t in group
                              if t["lid_stable"] and t["class"] in ("SW", "HYBRID")
                              for b in t["beta_brackets"]]
                candidates.sort(key=lambda x: x["beta_over_k0"])
                selected = candidates[0] if len(candidates) == 1 else None
                for t in group:
                    if t["lid_stable"] and t["class"] == "RES":
                        nearest = min(t["points"], key=lambda p: abs(p["f_ghz"] - TARGET_GHZ))
                        flat.append({"direction": direction, "track": t["track"],
                                     "f_ghz": nearest["f_ghz"], "q": nearest["q"],
                                     "distance_ghz": abs(nearest["f_ghz"] - TARGET_GHZ)})
                axes[direction] = {"bracketed_candidates": candidates, "selected": selected,
                                   "box_track_count": sum(t["box_mode"] for t in group)}
            flat.sort(key=lambda x: x["distance_ghz"])
            entry["states"][state] = {"axes": axes, "nearest_flat_track": flat[0] if flat else None}
            a, b = axes["along"]["selected"], axes["across"]["selected"]
            entry["states"][state]["anisotropy_across_minus_along"] = (
                b["beta_over_k0"] - a["beta_over_k0"] if a and b else None)
        off = entry["states"]["OFF"]["anisotropy_across_minus_along"]
        on = entry["states"]["ON"]["anisotropy_across_minus_along"]
        entry["anisotropy_ON_over_OFF"] = on / off if on is not None and off not in (None, 0) else None
        summary["by_geometry"][g.name] = entry
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
