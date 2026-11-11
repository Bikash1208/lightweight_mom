"""Build Cell 2 projects sequentially with new local AEDT sessions."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from cell2_geometry import GEOMETRIES, check
from cell2_spec import variant

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "lab_B_cell2_work"
RESULTS = REPO / "results" / "sub6_cell_param" / "lab_B_cell2"


def build_one(g):
    dc = check(g)
    if dc["status"] != "PASS" or not dc["only_diode_connects_target_nets"]:
        raise RuntimeError(f"DC gate failed {g.name}: {dc['conflicts']}")
    name = variant(g)
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
    cmd = [sys.executable, str(HERE / "cell2_cell.py"), "--geometry", g.name,
           "--lid-height-mm", "150", "--top-boundary", "PEC",
           "--build", "--project", str(project), "--output", str(out)]
    with (out / "build.log").open("w", encoding="utf-8") as log:
        rc = subprocess.run(cmd, cwd=launch, stdout=log, stderr=subprocess.STDOUT).returncode
    summary = json.loads(summary_file.read_text()) if summary_file.exists() else {}
    print(f"BUILD {name}: exit {rc}, status {summary.get('status', 'NO_SUMMARY')}", flush=True)
    if rc or summary.get("status") != "BUILT" or summary.get("aedt") != "2025.2":
        raise RuntimeError(f"BUILD {name} failed; inspect own build.log")
    if summary["mesh_operations"] or summary["floor_boundary"] != "SOLID_PEC_GROUND":
        raise RuntimeError(f"BUILD {name} wrong mesh/floor")
    required = {"LEFT", "RIGHT"} if g.family == "SD" else {"DISC"}
    if not required.issubset(summary["objects"]):
        raise RuntimeError(f"BUILD {name} required copper components missing")
    if summary["top_boundary"]["type"] != "PEC" or summary["top_boundary"]["height_mm"] != 150:
        raise RuntimeError(f"BUILD {name} wrong top boundary")
    if not summary["geometry"]["only_diode_connects_target_nets"] or summary["stack"]["vias"]:
        raise RuntimeError(f"BUILD {name} DC failure")
    dest = RESULTS / "builds" / f"{name}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(summary_file, dest)


def main():
    if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
        raise RuntimeError("Explicit AEDT 2025.2 selection required")
    stop = RESULTS / "status.json"
    if stop.exists() and json.loads(stop.read_text()).get("status", "").startswith("STOPPED"):
        raise RuntimeError("Cell 2 hard stop recorded; refusing additional AEDT builds")
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", choices=[g.name for g in GEOMETRIES])
    args = ap.parse_args()
    chosen = [g for g in GEOMETRIES if args.geometry is None or g.name == args.geometry]
    for g in chosen:
        build_one(g)


if __name__ == "__main__":
    main()
