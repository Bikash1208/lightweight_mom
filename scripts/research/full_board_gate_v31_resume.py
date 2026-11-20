"""Continue an owned first-pass full-board project to adaptive convergence."""
import argparse
import json
import os
import re
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import psutil

os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss


def stamp():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    path.write_text(json.dumps(value, indent=1, default=str), encoding="utf-8")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--state", required=True)
    args = p.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    temp = out / "aedt_temp"
    temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp))
    state = {"phase": "START", "state": args.state, "started_utc": stamp(), "controller_pid": os.getpid()}
    status_path = out / "status.json"
    save(status_path, state)
    h = None
    try:
        project = args.project.resolve()
        assert project.exists()
        h = Hfss(project=str(project), design=f"GATE_{args.state}", solution_type="Modal",
                 version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
        state["own_aedt_pid"] = h.desktop_class.aedt_process_id
        state["aedt_version"] = str(h.desktop_class.aedt_version_id)
        assert state["aedt_version"] == "2026.1"
        h.desktop_class.temp_directory = temp
        assert Path(str(h.desktop_class.temp_directory)).resolve() == temp.resolve()
        names = list(h.modeler.object_names)
        metals = [name for name in names if name.startswith("METAL_PEC")]
        state["object_counts"] = {"all": len(names), "metal": len(metals)}
        assert len(metals) == 518 and "SUBSTRATE" in names and "AIR" in names and "ABSORBER" in names
        assert str(h.modeler["SUBSTRATE"].material_name).lower() == "rt5880_gate"
        assert str(h.modeler["ABSORBER"].material_name).lower() == "absorber_15_10"
        assert len(h.setups) == 1 and h.setups[0].name == "Setup1"
        setup = h.setups[0]
        assert str(setup.props["Frequency"]) == "5.68GHz"
        assert int(setup.props["MaximumPasses"]) == 1
        setup.props["MaximumPasses"] = 12
        setup.props["MinimumPasses"] = 6
        setup.props["MaxDeltaS"] = 0.02
        assert setup.update()
        state["setup_readback"] = {key: str(setup.props.get(key)) for key in
                                   ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS", "DrivenSolverType")}
        assert int(setup.props["MaximumPasses"]) == 12
        assert h.validate_simple(str(out / "validation.log"))
        assert h.save_project()
        state["phase"] = "SOLVING"
        state["solve_started_utc"] = stamp()
        state["cores_requested"] = max(1, psutil.cpu_count(logical=True) - 2)
        save(status_path, state)
        start = time.monotonic()
        assert h.analyze_setup("Setup1", cores=state["cores_requested"], use_auto_settings=False, blocking=True)
        state["wall_s"] = round(time.monotonic() - start, 3)
        state["solve_ended_utc"] = stamp()
        f = Path(h.export_profile("Setup1", output_file=str(out / "solve.prof")))
        profile = f.read_text(errors="replace")
        state["profile_path"] = str(f)
        state["hfss_max_memory_gb"] = max(map(float, re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", profile)))
        state["max_tets"] = max(map(int, re.findall(r"Tetrahedra:\s*(\d+)", profile)))
        state["adaptive_passes"] = len(re.findall(r"^\s*Adaptive Pass \d+\s+Frequency:", profile, flags=re.M))
        state["normal_completion"] = "Status: Normal Completion" in profile
        # export_profile concatenates the earlier one-pass run and this continuation.
        # The final solution summary, not an earlier summary, is authoritative.
        state["converged"] = profile.rfind("Adaptive Passes converged") > profile.rfind("Adaptive Passes did not converge")
        assert state["normal_completion"] and state["converged"], "SOLVE_DID_NOT_CONVERGE"
        assert h.save_project()
        state["phase"] = "COMPLETE"
        save(status_path, state)
        print("COMPLETE", state["adaptive_passes"], state["hfss_max_memory_gb"], flush=True)
    except BaseException:
        state["phase"] = "FAILED"
        state["traceback"] = traceback.format_exc()
        save(status_path, state)
        print(state["traceback"], file=sys.stderr, flush=True)
        raise
    finally:
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
