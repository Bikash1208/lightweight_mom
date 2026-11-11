"""Build the ten Option 1 projects sequentially after geometry and VERIFY gates."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from opt1_geometry import all_geometries

here = Path(__file__).resolve().parent
repo = here.parents[1]
work = here / "lab_B_opt1_work"
results = repo / "results" / "sub6_cell_param" / "lab_B_opt1"
if json.loads((results / "geometry_precheck.json").read_text())["status"] != "PASS":
    raise RuntimeError("Geometry gate missing or failed")
if json.loads((results / "verify" / "determinism.json").read_text())["status"] != "PASS":
    raise RuntimeError("Determinism gate missing or failed")
if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
    raise RuntimeError("Explicit AEDT 2025.2 selection required")

for item in all_geometries(0.5):
    name = item["name"]
    project = work / "projects" / f"{name}.aedt"
    out = work / "build" / name
    summary_file = out / "build_summary.json"
    if project.exists() and summary_file.exists() and json.loads(summary_file.read_text()).get("status") == "BUILT":
        print(f"BUILD {name}: already BUILT", flush=True)
        continue
    if project.exists():
        raise RuntimeError(f"Refusing to overwrite incomplete project {name}")
    out.mkdir(parents=True, exist_ok=True)
    project.parent.mkdir(parents=True, exist_ok=True)
    launch = work / "launch" / name
    launch.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(here / "opt1_cell.py"), "--geometry", name,
           "--wide-channel", "0.5", "--build", "--project", str(project), "--output", str(out)]
    with (out / "build.log").open("w", encoding="utf-8") as log:
        rc = subprocess.run(cmd, cwd=launch, stdout=log, stderr=subprocess.STDOUT).returncode
    summary = json.loads(summary_file.read_text()) if summary_file.exists() else {}
    print(f"BUILD {name}: exit {rc}, status {summary.get('status', 'NO_SUMMARY')}", flush=True)
    if rc or summary.get("status") != "BUILT" or summary.get("aedt") != "2025.2":
        raise RuntimeError(f"BUILD {name} failed; inspect own build.log")
    if summary["mesh_operations"] or summary["floor_boundary"] != "SOLID_PEC_GROUND":
        raise RuntimeError(f"BUILD {name} has wrong mesh/floor")
    if summary["stack"]["vias"] or summary["stack"]["fr4"] or summary["stack"]["air_below_ground"]:
        raise RuntimeError(f"BUILD {name} has wrong stack")
    dest = results / "builds" / f"{name}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(summary_file, dest)
