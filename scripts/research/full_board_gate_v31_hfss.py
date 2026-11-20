"""Build and solve one v3.1 R150 full-board gate state in a fresh HFSS 2026 R1 process."""
import argparse
import json
import os
import re
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import psutil

os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss


def now():
    return datetime.now(timezone.utc).isoformat()


def write(path, obj):
    path.write_text(json.dumps(obj, indent=1, default=str), encoding="utf-8")


def need(value, label):
    if not value:
        raise RuntimeError(label)
    return value


def profile(path):
    text = path.read_text(errors="replace")
    tets = [int(v) for v in re.findall(r"Tetrahedra:\s*(\d+)", text)]
    memory = [float(v) for v in re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", text)]
    delta = [float(v) for v in re.findall(r"Max Mag\. Delta S:\s*([\d.]+)", text)]
    return {"tetrahedra": tets, "hfss_max_memory_gb": max(memory) if memory else None,
            "converged": text.rfind("Adaptive Passes converged") > text.rfind("Adaptive Passes did not converge"),
            "normal_completion": "Status: Normal Completion" in text,
            "delta_s": delta, "profile_excerpt": [line.strip() for line in text.splitlines()
            if "Max memory/process" in line or "Tetrahedra:" in line or "Delta S:" in line][-30:]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--state", required=True)
    parser.add_argument("--first-pass", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    temp = out / "aedt_temp"
    temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp))
    status_path = out / "status.json"
    status = {"state": args.state, "phase": "START", "start_utc": now(),
              "first_pass_only": args.first_pass, "controller_pid": os.getpid(),
              "aedt_pid": None, "aedt_version": None}
    write(status_path, status)
    done = threading.Event()

    def heartbeat():
        while not done.wait(60):
            status["heartbeat_utc"] = now()
            try:
                if status["aedt_pid"]:
                    proc = psutil.Process(status["aedt_pid"])
                    tree = [proc] + proc.children(recursive=True)
                    status["own_process_tree_gib"] = round(sum(p.memory_info().rss for p in tree if p.is_running()) / 2**30, 3)
                status["free_memory_gib"] = round(psutil.virtual_memory().available / 2**30, 3)
            except psutil.Error:
                pass
            write(status_path, status)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    h = None
    try:
        project = out / f"BOARD_R150_GATE_V31_{args.state}.aedt"
        need(not project.exists(), "PROJECT_ALREADY_EXISTS_NEW_TAG_REQUIRED")
        h = Hfss(project=str(project), design=f"GATE_{args.state}", solution_type="Modal",
                 version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
        status["aedt_pid"] = h.desktop_class.aedt_process_id
        status["aedt_version"] = str(h.desktop_class.aedt_version_id)
        need(status["aedt_version"] == "2026.1", "WRONG_AEDT_VERSION")
        h.desktop_class.temp_directory = temp
        status["aedt_temp_readback"] = str(h.desktop_class.temp_directory)
        need(Path(status["aedt_temp_readback"]).resolve() == temp.resolve(), "TEMP_READBACK")
        status["phase"] = "IMPORT"
        write(status_path, status)
        m = h.modeler
        m.model_units = "mm"
        rt = h.materials.add_material("RT5880_GATE")
        rt.permittivity = 2.2
        rt.dielectric_loss_tangent = 0.0009
        absorber = h.materials.add_material("ABSORBER_15_10")
        absorber.permittivity = 1.5
        absorber.dielectric_loss_tangent = 1 / 1.5
        absorber.permeability = 1.5
        absorber.magnetic_loss_tangent = 1 / 1.5
        need(m.import_3d_cad(str(args.step.resolve())), "STEP_IMPORT_FAILED")
        names = list(m.object_names)
        metals = [n for n in names if n.startswith("METAL_PEC")]
        need(len(metals) > 0 and "SUBSTRATE" in names and "BIAS_RT5880" in names
             and "AIR" in names and "PORT" in names, f"OBJECT_READBACK_FAILED:{len(metals)}")
        for name in metals:
            m[name].material_name = "pec"
        m["SUBSTRATE"].material_name = "RT5880_GATE"
        m["BIAS_RT5880"].material_name = "RT5880_GATE"
        m["AIR"].material_name = "vacuum"
        status["object_counts"] = {"all": len(names), "metal": len(metals)}
        status["materials_readback"] = {name: str(m[name].material_name) for name in
                                        (metals[0], "SUBSTRATE", "BIAS_RT5880", "AIR")}
        need(h.lumped_port(assignment="PORT", integration_line=[[0.635, 0, -3.175], [2.05, 0, -3.175]],
                           impedance=50, name="Port1", renormalize=True, deembed=False), "PORT_FAILED")
        need(h.assign_radiation_boundary_to_objects("AIR"), "RADIATION_FAILED")
        outer = need(m.create_cylinder("Z", [0, 0, 0], 186, 30.035,
                                       name="ABS_OUT", material="ABSORBER_15_10"), "ABS_OUT_FAILED")
        inner = need(m.create_cylinder("Z", [0, 0, -1], 156, 32,
                                       name="ABS_IN", material="vacuum"), "ABS_IN_FAILED")
        need(m.subtract(outer, inner, keep_originals=False), "ABS_SUBTRACT_FAILED")
        m[outer.name].name = "ABSORBER"
        m["ABSORBER"].material_name = "ABSORBER_15_10"
        status["absorber_readback"] = str(m["ABSORBER"].material_name)
        setup = need(h.create_setup("Setup1", Frequency="5.68GHz", MaximumPasses=1 if args.first_pass else 12,
                                    MaxDeltaS=0.02, MinimumPasses=1 if args.first_pass else 6), "SETUP_FAILED")
        setup.props["MaximumPasses"] = 1 if args.first_pass else 12
        setup.props["MinimumPasses"] = 1 if args.first_pass else 6
        setup.props["MaxDeltaS"] = 0.02
        setup.props["BasisOrder"] = 1
        # Try the licensed iterative solver on the measured first-pass gate.
        setup.props["DrivenSolverType"] = "Iterative"
        status["iterative_update_accepted"] = bool(setup.update())
        if not status["iterative_update_accepted"]:
            setup.props["DrivenSolverType"] = "Direct"
            need(setup.update(), "DIRECT_SETUP_UPDATE_FAILED")
        status["solver_type_readback"] = str(setup.props.get("DrivenSolverType"))
        need(h.insert_infinite_sphere(theta_start=0, theta_stop=180, theta_step=2,
                                      phi_start=0, phi_stop=360, phi_step=2, name="Global"), "SPHERE_FAILED")
        status["setup_readback"] = {key: str(setup.props.get(key)) for key in
                                    ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS", "DrivenSolverType")}
        status["boundaries_readback"] = [b.name for b in h.boundaries]
        status["spheres_readback"] = [v.name for v in h.field_setups]
        valid = h.validate_simple(str(out / "validation.log"))
        status["validation"] = valid
        need(valid == 1 or valid is True, "VALIDATION_FAILED")
        need(h.save_project(), "SAVE_FAILED")
        status["phase"] = "SOLVING_FIRST_PASS" if args.first_pass else "SOLVING"
        status["solve_start_utc"] = now()
        write(status_path, status)
        t0 = time.monotonic()
        cores = max(1, psutil.cpu_count(logical=True) - 2)
        status["cores_requested"] = cores
        need(h.analyze_setup("Setup1", cores=cores, use_auto_settings=False, blocking=True), "SOLVE_FAILED")
        status["wall_s"] = round(time.monotonic() - t0, 3)
        status["solve_end_utc"] = now()
        prof = Path(h.export_profile("Setup1", output_file=str(out / "solve.prof")))
        status["profile"] = profile(prof)
        if not args.first_pass:
            need(status["profile"]["converged"] and status["profile"]["normal_completion"],
                 "SOLVE_DID_NOT_CONVERGE")
        status["profile_path"] = str(prof)
        status["free_memory_gib_after"] = round(psutil.virtual_memory().available / 2**30, 3)
        need(h.save_project(), "POST_SOLVE_SAVE_FAILED")
        status["phase"] = "COMPLETE_FIRST_PASS" if args.first_pass else "COMPLETE"
        write(status_path, status)
        print(status["phase"], status["profile"], flush=True)
    except BaseException:
        status["phase"] = "FAILED"
        status["traceback"] = traceback.format_exc()
        write(status_path, status)
        print(status["traceback"], file=sys.stderr, flush=True)
        raise
    finally:
        done.set()
        thread.join(timeout=5)
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
