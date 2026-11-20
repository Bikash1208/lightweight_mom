"""Read-only probe of stored driven fields in an owned v3.1 full-board project."""
import json
import os
import sys
import traceback
from pathlib import Path

os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

project = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
out.mkdir(parents=True, exist_ok=True)
temp = out / "aedt_temp"
temp.mkdir(exist_ok=True)
os.environ.update(TEMP=str(temp), TMP=str(temp))
report = {"status": "START", "project": str(project)}
h = None
try:
    h = Hfss(project=str(project), design=f"GATE_{sys.argv[3]}", solution_type="Modal",
             version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    report["own_aedt_pid"] = h.desktop_class.aedt_process_id
    h.desktop_class.temp_directory = temp
    report["sweeps"] = h.existing_analysis_sweeps
    report["setup"] = {k: str(h.setups[0].props.get(k)) for k in ("Frequency", "MaximumPasses", "SaveFields")}
    points = [[40, 0, 2.035], [60, 0, 2.035]]
    for quantity in ("Mag_E", "E"):
        f = out / f"probe_{quantity}.fld"
        report[quantity] = bool(h.post.export_field_file(quantity=quantity, solution="Setup1 : LastAdaptive",
                                  output_file=str(f), sample_points=points, export_in_si_system=False,
                                  intrinsics={"Freq": "5.68GHz", "Phase": "0deg"}))
        report[quantity + "_sample"] = f.read_text(errors="replace")[:700] if f.exists() else None
    report["status"] = "COMPLETE"
except BaseException:
    report["status"] = "FAILED"
    report["traceback"] = traceback.format_exc()
finally:
    (out / "probe.json").write_text(json.dumps(report, indent=2))
    if h is not None:
        h.release_desktop(close_projects=True, close_desktop=True)
print(report["status"], flush=True)
