"""Build Option 4 projects sequentially with new local AEDT sessions."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from opt4_geometry import GEOMETRIES, check
from opt4_spec import LIDS, variant

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "lab_B_opt4_work"
RESULTS = REPO / "results" / "sub6_cell_param" / "lab_B_opt4"


def build_one(g, lid):
    dc = check(g)
    if dc["status"] != "PASS" or not dc["anode_reaches_ground_only_through_diode"]:
        raise RuntimeError(f"DC gate failed {g.name}: {dc['conflicts']}")
    height, kind = LIDS[lid]
    name = variant(g, lid)
    project = WORK / "projects" / f"{name}.aedt"
    out = WORK / "build" / name
    summary_file = out / "build_summary.json"
    if project.exists() and summary_file.exists() and json.loads(summary_file.read_text()).get("status") == "BUILT":
        print(f"BUILD {name}: already BUILT", flush=True)
        return
    if project.exists():
        raise RuntimeError(f"Refusing to overwrite incomplete project {name}")
    out.mkdir(parents=True, exist_ok=True)
    project.parent.mkdir(parents=True, exist_ok=True)
    launch = WORK / "launch" / name
    launch.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(HERE / "opt4_cell.py"), "--geometry", g.name,
           "--lid-height-mm", str(height), "--top-boundary", kind,
           "--build", "--project", str(project), "--output", str(out)]
    with (out / "build.log").open("w", encoding="utf-8") as log:
        rc = subprocess.run(cmd, cwd=launch, stdout=log, stderr=subprocess.STDOUT).returncode
    summary = json.loads(summary_file.read_text()) if summary_file.exists() else {}
    print(f"BUILD {name}: exit {rc}, status {summary.get('status', 'NO_SUMMARY')}", flush=True)
    if rc or summary.get("status") != "BUILT" or summary.get("aedt") != "2025.2":
        raise RuntimeError(f"BUILD {name} failed; inspect own build.log")
    if summary["mesh_operations"] or summary["floor_boundary"] != "SOLID_PEC_GROUND":
        raise RuntimeError(f"BUILD {name} wrong mesh/floor")
    if not {"DISC", "ISLAND"}.issubset(summary["objects"]):
        raise RuntimeError(f"BUILD {name} required copper components missing")
    if summary["top_boundary"]["type"] != kind or summary["top_boundary"]["height_mm"] != height:
        raise RuntimeError(f"BUILD {name} wrong top boundary")
    if kind == "IMPEDANCE_377" and (summary["top_boundary"]["Resistance"],
                                      summary["top_boundary"]["Reactance"]) != ("377ohm", "0ohm"):
        raise RuntimeError(f"BUILD {name} impedance readback mismatch")
    if not summary["geometry"]["anode_reaches_ground_only_through_diode"]:
        raise RuntimeError(f"BUILD {name} DC failure")
    dest = RESULTS / "builds" / f"{name}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(summary_file, dest)


def main():
    if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
        raise RuntimeError("Explicit AEDT 2025.2 selection required")
    ap = argparse.ArgumentParser()
    ap.add_argument("--lid", choices=LIDS, default="PEC150")
    ap.add_argument("--geometry", choices=[g.name for g in GEOMETRIES])
    args = ap.parse_args()
    chosen = [g for g in GEOMETRIES if args.geometry is None or g.name == args.geometry]
    for g in chosen:
        build_one(g, args.lid)


if __name__ == "__main__":
    main()
