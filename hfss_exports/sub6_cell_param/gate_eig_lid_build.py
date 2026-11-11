"""Build private 40/300 mm PEC-lid SD projects after the 150 mm queue."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from gate_eig_lid_spec import GEOMETRY, variant
from cell2_geometry import check

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "lab_B_gate_eig_work"
RESULTS = REPO / "results" / "sub6_cell_param" / "lab_B_gate_eig"


def build(lid):
    name = variant(lid)
    project = WORK / "projects" / f"{name}.aedt"
    output = WORK / "build" / name
    summary_file = output / "build_summary.json"
    if project.exists() and summary_file.exists() and json.loads(summary_file.read_text()).get("status") == "BUILT":
        print(f"BUILD {name}: already BUILT", flush=True)
        return
    if project.exists():
        raise RuntimeError(f"Refusing to overwrite incomplete project {name}")
    output.mkdir(parents=True, exist_ok=True)
    project.parent.mkdir(parents=True, exist_ok=True)
    launch = WORK / "launch" / name
    launch.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(HERE / "gate_eig_lid_cell.py"), "--geometry", GEOMETRY.name,
           "--project", str(project), "--output", str(output),
           "--lid-height-mm", str(lid), "--build"]
    with (output / "build.log").open("w", encoding="utf-8") as log:
        rc = subprocess.run(cmd, cwd=launch, stdout=log, stderr=subprocess.STDOUT).returncode
    summary = json.loads(summary_file.read_text()) if summary_file.exists() else {}
    print(f"BUILD {name}: exit {rc}, status {summary.get('status', 'NO_SUMMARY')}", flush=True)
    if rc or summary.get("status") != "BUILT" or summary.get("aedt") != "2025.2":
        raise RuntimeError(f"BUILD {name} failed; inspect own build.log")
    if summary["geometry"]["status"] != "PASS" or summary["stack"]["vias"]:
        raise RuntimeError(f"BUILD {name} geometry/stack gate failed")
    if not {"LEFT", "RIGHT"}.issubset(summary["objects"]):
        raise RuntimeError(f"BUILD {name} split copper missing")
    if summary["top_boundary"] != {"type": "PEC", "height_mm": float(lid)}:
        raise RuntimeError(f"BUILD {name} lid readback mismatch")
    dest = RESULTS / "builds" / f"{name}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(summary_file, dest)


if __name__ == "__main__":
    if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
        raise RuntimeError("Explicit AEDT 2025.2 selection required")
    status = RESULTS / "status.json"
    if status.exists() and json.loads(status.read_text()).get("status") == "ABORTED_BY_USER":
        raise RuntimeError("Gate is closed; refusing additional AEDT builds")
    if check(GEOMETRY)["status"] != "PASS":
        raise RuntimeError("SD geometry gate failed")
    for height in (40, 300):
        build(height)
