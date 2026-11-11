"""Copy only the two validated SD yline projects into this job's private workspace."""
import json
import os
from pathlib import Path
import shutil

from gate_eig_spec import GEOMETRIES, STATES, PHASES_HIGH, PHASES_C_LOW, cases, variant

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE_WORK = HERE / "lab_B_cell2_work" / "projects"
SOURCE_RESULTS = REPO / "results" / "sub6_cell_param" / "lab_B_cell2" / "builds"
WORK = HERE / "lab_B_gate_eig_work"
RESULTS = REPO / "results" / "sub6_cell_param" / "lab_B_gate_eig"

if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
    raise RuntimeError("Explicit AEDT 2025.2 selection required")

sources = []
for g in GEOMETRIES:
    name = variant(g)
    summary = json.loads((SOURCE_RESULTS / f"{name}.json").read_text())
    if (summary["status"] != "BUILT" or summary["aedt"] != "2025.2" or
            summary["geometry"]["status"] != "PASS" or
            summary["top_boundary"] != {"type": "PEC", "height_mm": 150.0} or
            summary["stack"]["vias"] or not {"LEFT", "RIGHT"}.issubset(summary["objects"])):
        raise RuntimeError(f"Source build gate failed: {name}")
    source = SOURCE_WORK / f"{name}.aedt"
    target = WORK / "projects" / f"{name}.aedt"
    if not source.exists() or target.exists():
        raise RuntimeError(f"Source missing or private target already exists: {name}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    dest_summary = RESULTS / "builds" / f"{name}.json"
    dest_summary.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE_RESULTS / f"{name}.json", dest_summary)
    sources.append({"geometry": g.name, "private_project": str(target.relative_to(REPO)),
                    "source_build": str((SOURCE_RESULTS / f"{name}.json").relative_to(REPO)),
                    "rows": len(cases(g))})
RESULTS.mkdir(parents=True, exist_ok=True)
(RESULTS / "matrix.json").write_text(json.dumps({"source_builds": sources, "states": STATES,
    "high_phases_deg": PHASES_HIGH, "extra_C_phases_deg": PHASES_C_LOW,
    "unique_rows": sum(len(cases(g)) for g in GEOMETRIES), "direction": "across",
    "phase_variable": "$ph_x", "other_phase_variable": "$ph_y=0",
    "C_states_ohm_per_square": {"OFF_C013": [29.56, -390.9], "OFF_C018": [15.45, -275.0]}},
    indent=2) + "\n")
print(f"Prepared {len(sources)} private project copies and {sum(x['rows'] for x in sources)} unique rows")
