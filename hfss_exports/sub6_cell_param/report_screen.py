"""Merge SCREEN run folders, re-score everything, and write the decision-session report.

python report_screen.py --runs screen_bar90_nostub_run screen_bar90_trace_run screen_bar90_full_run --out screen_bar90_report [--diode bar90]

Outputs (in --out): analysis.json (full), SCREEN_classified_modes.csv (one row per geometry/state/direction/mode),
SCREEN_geometry_summary.csv (pass/fail, anisotropy), SCREEN_stub_spread.csv (per mode rank, spread across stub 0/0.3/1/2 + trace),
SCREEN_REPORT.md (compact tables). Rows without a matching AEDT read-back are REFUSED by the analysis.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cell_param_pyaedt as cp  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--diode", default="bar90")
    a = ap.parse_args()
    cp.set_profile(a.diode)
    a.out.mkdir(parents=True, exist_ok=True)
    merged = a.out / "_merged_rows"
    if merged.exists():
        shutil.rmtree(merged)
    merged.mkdir()
    n = 0
    for run in a.runs:
        run = run if run.is_absolute() else HERE / run
        for js in run.glob("*.json"):
            try:
                if "modes" in json.loads(js.read_text()):
                    shutil.copy(js, merged / js.name); n += 1
            except Exception:
                continue
    rep = cp.analyse(merged)
    shutil.move(str(merged / "analysis.json"), str(a.out / "analysis.json"))
    nominal = {k: v for k, v in cp.nominal().items() if k not in {"$ph_x", "$ph_y", "Rs", "Xs"}}

    def label(g):
        g = dict(g); var = g.pop("_project_variant")
        diffs = [f"{k}={g[k]:g}" for k in sorted(g) if k in nominal and abs(g[k]-nominal[k]) > 1e-9
                 and k not in ("island_via_x", "island_via_y")]
        return var, (" ".join(diffs) or "NOMINAL")

    # classified modes
    with (a.out / "SCREEN_classified_modes.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["variant", "geometry", "state_direction", "mode90", "f80", "f90", "f100", "q90", "slope_GHz_per_deg", "class",
                    "ambiguous_track", "dist_from_5p8_GHz", "counts_for_criterion", "violates_1GHz", "beta_over_k0_5p8", "beta_status"])
        for e in rep["ranking"]:
            var, lab = label(e["geometry"])
            for sd, c in sorted(e["classification"].items()):
                for m in c.get("modes", []):
                    w.writerow([var, lab, sd, m["mode90"], m.get("f80_ghz"), m["f90_ghz"], m.get("f100_ghz"), m["q90"], m.get("slope_ghz_per_deg"),
                                m["klass"], m.get("ambiguous_track"), m["dist_from_5p8_ghz"], m["counts_for_1ghz_criterion"], m["violates_1ghz"],
                                m.get("beta_over_k0_5p8"), m.get("beta_status")])
    # geometry summary
    with (a.out / "SCREEN_geometry_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["variant", "geometry", "criterion_1GHz", "awaiting_80_100", "violations", "OFF_aniso", "ON_aniso", "ON_over_OFF", "beta_OFF", "beta_ON"])
        for e in rep["ranking"]:
            var, lab = label(e["geometry"]); an = e["anisotropy"]
            w.writerow([var, lab, e["criterion_1ghz"], ";".join(e["awaiting_80_100"]),
                        "; ".join(f"{v['case']} {v['f90_ghz']}GHz {v['klass']}" for v in e["violations_1ghz"]),
                        an.get("OFF", {}).get("anisotropy"), an.get("ON", {}).get("anisotropy"), an.get("ON_over_OFF_anisotropy_ratio"),
                        json.dumps(an.get("OFF", {}).get("beta_over_k0_by_direction")), json.dumps(an.get("ON", {}).get("beta_over_k0_by_direction"))])
    # stub / trace spread
    with (a.out / "SCREEN_stub_spread.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["geometry", "state_direction", "mode_rank", "spread_GHz", "FLAG_gt_0p3GHz", "f_by_label"])
        for g in rep["bias_path_stub_spread"]:
            lab = " ".join(f"{k}={v:g}" for k, v in g["geometry"].items() if k in nominal and abs(v-nominal[k]) > 1e-9
                           and k not in ("island_via_x", "island_via_y")) or "NOMINAL"
            for sd, per in sorted(g["by_state"].items()):
                for rk, v in sorted(per.items()):
                    w.writerow([lab, sd, rk, v["spread_ghz"], v["FLAG_spread_gt_0p3_ghz"], json.dumps(v["f_ghz_by_label"])])
    # compact markdown
    L = [f"# SCREEN report (merged {n} solved rows; criterion = >=1 GHz from 5.8 GHz for RESONANCE and MIXED modes)", ""]
    L.append("## Geometry pass / fail")
    L.append("| variant | geometry | criterion | violations | OFF aniso | ON aniso |")
    L.append("|---|---|---|---|---|---|")
    for e in rep["ranking"]:
        var, lab = label(e["geometry"]); an = e["anisotropy"]
        L.append(f"| {var} | {lab} | {e['criterion_1ghz']} | " + "; ".join(f"{v['case']} {v['f90_ghz']} {v['klass']}" for v in e["violations_1ghz"]) +
                 f" | {an.get('OFF', {}).get('anisotropy')} | {an.get('ON', {}).get('anisotropy')} |")
    L += ["", "## Stub 0 / 0.3 / 1 / 2 mm + trace spread (key result; FLAG if > 0.3 GHz)", "| geometry | state/dir | mode rank | spread GHz | flag | f by label |", "|---|---|---|---|---|---|"]
    for g in rep["bias_path_stub_spread"]:
        lab = " ".join(f"{k}={v:g}" for k, v in g["geometry"].items() if k in nominal and abs(v-nominal[k]) > 1e-9
                       and k not in ("island_via_x", "island_via_y")) or "NOMINAL"
        for sd, per in sorted(g["by_state"].items()):
            for rk, v in sorted(per.items()):
                L.append(f"| {lab} | {sd} | {rk} | {v['spread_ghz']} | {'**FLAG**' if v['FLAG_spread_gt_0p3_ghz'] else ''} | {json.dumps(v['f_ghz_by_label'])} |")
    (a.out / "SCREEN_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"report written to {a.out} ({n} rows, {len(rep['ranking'])} geometry groups, {len(rep['bias_path_stub_spread'])} spread sets)")


if __name__ == "__main__":
    main()
