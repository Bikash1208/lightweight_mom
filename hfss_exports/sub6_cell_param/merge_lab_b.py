"""Merge QUEUE_B results from the lab machine (results/sub6_cell_param/lab_B/) into this machine's run folders, through the SAME gates.

A lab row is ACCEPTED only if ALL hold (otherwise it is rejected with a reason in lab_B_merge_report.json):
  * the lab's VERIFY passed the main-server comparison (recomputed here against verify_reference/, f within 0.5 %)
  * run_id is a row of SCREEN_2.csv / SCREEN_6_short_nostub.csv / SCREEN_6_short_trace.csv (QUEUE_B), and its project_variant matches that file
  * the stored variables are IDENTICAL to the CSV row (geometry, phases, Rs/Xs)
  * the AEDT read-back stored with the row matches the intended variables (cell_param_pyaedt.readback_ok) and the floor depth is 30 mm (air_below, AIR bottom -34.01 mm)
  * the row has modes (a solved row), not a FAILED stub
Accepted rows are copied into screen_bar90_full_run / screen_bar90_nostub_run / screen_bar90_trace_run; FAILED lists and parallel reports are merged so the
analysis and the FAILED-row report see them. Then re-run report_screen.py.

python merge_lab_b.py --lab-dir <path>/results/sub6_cell_param/lab_B [--dry-run]
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cell_param_pyaedt as cp  # noqa: E402

QUEUE_B = {  # csv stem -> (expected project_variant, destination run dir)
    "SCREEN_2": ("full", HERE / "screen_bar90_full_run"),
    "SCREEN_6_short_nostub": ("nostub", HERE / "screen_bar90_nostub_run"),
    "SCREEN_6_short_trace": ("trace", HERE / "screen_bar90_trace_run"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lab-dir", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    cp.set_profile("bar90")
    rep = {"accepted": 0, "rejected": [], "per_file": {}}
    # 1) the lab VERIFY must pass against OUR reference (recomputed, not trusted)
    lv = a.lab_dir / "verify" / "verify.json"
    rc = subprocess.run([sys.executable, "-X", "utf8", str(HERE / "check_verify_vs_reference.py"), "--lab", str(lv)], capture_output=True).returncode if lv.exists() else 1
    if rc != 0:
        print("LAB VERIFY missing or NOT within 0.5 % of the reference -> nothing merged"); sys.exit(1)
    # 2) expected rows
    expected = {}
    for stem in QUEUE_B:
        for run_id, v in cp.read_rows(HERE / f"{stem}.csv"):
            expected[run_id] = (stem, v)
    # 3) rows
    for js in sorted((a.lab_dir / "rows").glob("*.json")):
        try:
            r = json.loads(js.read_text())
        except Exception as exc:
            rep["rejected"].append({"file": js.name, "why": f"unreadable: {exc!r}"}); continue
        rid = r.get("run_id")
        why = None
        if rid not in expected:
            why = "not a QUEUE_B row"
        else:
            stem, v = expected[rid]
            variant, dest = QUEUE_B[stem]
            if "modes" not in r:
                why = "no modes"
            elif r.get("project_variant") != variant:
                why = f"project_variant {r.get('project_variant')} != {variant}"
            elif r.get("variables") != v:
                why = "variables differ from the CSV row"
            elif abs(r["variables"].get("air_below", -1) - 30.0) > 1e-9:
                why = "air_below != 30 mm"
            else:
                ok, reason = cp.readback_ok(r)
                if not ok:
                    why = "read-back mismatch: " + reason
                elif abs(r["readback_si"].get("air_bottom", 0) * 1e3 + 34.01) > 1e-6:
                    why = "AIR bottom read-back != -34.01 mm"
        if why:
            rep["rejected"].append({"file": js.name, "why": why}); continue
        rep["accepted"] += 1; rep["per_file"][stem] = rep["per_file"].get(stem, 0) + 1
        if not a.dry_run:
            dest.mkdir(exist_ok=True); shutil.copy(js, dest / js.name)
    # 4) FAILED lists and parallel reports
    if not a.dry_run:
        for stem, (variant, dest) in QUEUE_B.items():
            ff = a.lab_dir / f"FAILED_ROWS_{stem}.json"
            if ff.exists():
                dest.mkdir(exist_ok=True)
                tgt = dest / f"FAILED_ROWS_{stem}.json"
                prev = json.loads(tgt.read_text()) if tgt.exists() else []
                have = {x["run_id"] for x in prev}
                tgt.write_text(json.dumps(prev + [x for x in json.loads(ff.read_text()) if x["run_id"] not in have], indent=2))
        for name, (stem) in (("full", "SCREEN_2"), ("nostub", "SCREEN_6_short_nostub"), ("trace", "SCREEN_6_short_trace")):
            pr = a.lab_dir / f"parallel_report_{name}.json"
            if pr.exists():
                dest = QUEUE_B[stem][1]; dest.mkdir(exist_ok=True)
                tgt = dest / "parallel_report.json"
                cur = json.loads(tgt.read_text()) if tgt.exists() else []
                cur += [dict(x, source="lab_B") for x in json.loads(pr.read_text()) if x.get("file") == f"{stem}.csv"]
                tgt.write_text(json.dumps(cur, indent=2))
    (HERE / "lab_B_merge_report.json").write_text(json.dumps(rep, indent=2))
    print(f"accepted {rep['accepted']} rows ({rep['per_file']}); rejected {len(rep['rejected'])}")
    for x in rep["rejected"][:10]:
        print("  REJECTED", x["file"], "-", x["why"])
    sys.exit(0 if not rep["rejected"] else 2)


if __name__ == "__main__":
    main()
