"""Build the nine Option 1b projects sequentially in new local AEDT sessions."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from opt1b_spec import GEOMETRIES, LIDS, variant

here = Path(__file__).resolve().parent
repo = here.parents[1]
work = here / "lab_B_opt1b_work"
results = repo / "results" / "sub6_cell_param" / "lab_B_opt1b"
if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
    raise RuntimeError("Explicit AEDT 2025.2 selection required")

for g in GEOMETRIES:
    for lid, height, kind in LIDS:
        name = variant(g, lid)
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
        cmd = [sys.executable, str(here / "opt1_cell.py"), "--geometry", g.name,
               "--wide-channel", "0.5", "--lid-height-mm", str(height),
               "--top-boundary", kind, "--build", "--project", str(project),
               "--output", str(out)]
        with (out / "build.log").open("w", encoding="utf-8") as log:
            rc = subprocess.run(cmd, cwd=launch, stdout=log, stderr=subprocess.STDOUT).returncode
        summary = json.loads(summary_file.read_text()) if summary_file.exists() else {}
        print(f"BUILD {name}: exit {rc}, status {summary.get('status', 'NO_SUMMARY')}", flush=True)
        if rc or summary.get("status") != "BUILT" or summary.get("aedt") != "2025.2":
            raise RuntimeError(f"BUILD {name} failed; inspect own build.log")
        if summary["mesh_operations"] or summary["floor_boundary"] != "SOLID_PEC_GROUND":
            raise RuntimeError(f"BUILD {name} wrong mesh/floor")
        if not {"DISC", "ISLAND"}.issubset(summary["objects"]):
            raise RuntimeError(f"BUILD {name} local copper components not separate")
        if summary["top_boundary"]["type"] != kind or summary["top_boundary"]["height_mm"] != height:
            raise RuntimeError(f"BUILD {name} wrong top boundary")
        if kind == "IMPEDANCE_377" and (summary["top_boundary"]["Resistance"], summary["top_boundary"]["Reactance"]) != ("377ohm", "0ohm"):
            raise RuntimeError(f"BUILD {name} impedance readback mismatch")
        dest = results / "builds" / f"{name}.json"
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(summary_file, dest)
