"""One fresh 2026.1 rotated anchor-strip supercell; own AEDT PID only."""

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone

import psutil

os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

from sub6_rotated_supercell_hfss import CU, P, RT, BIAS, modes, need

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sub6_meander_model" / "diode_supercells"
WORK = ROOT / "hfss_work" / "diode_supercells"
R = 5.45
SLOT_L = 8.72
SLOT_W = 2.91
STRIP_Y = 0.5
GAP_X = 0.3
VALUES = {"ON_A": (2.0, -139.6), "ON_B": (22.23, -335.4),
          "OFF_D": (1.6, 29.2)}


def utc():
    return datetime.now(timezone.utc).isoformat()


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=str) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def build(h, state, phase=30, element="RLC"):
    m = h.modeler
    m.model_units = "mm"
    mat = h.materials.add_material("RT5880_SUPERCELL")
    mat.permittivity = 2.2
    mat.dielectric_loss_tangent = 0.0009
    h["$ph_s"] = f"{phase}deg"
    h["$ph_y"] = "0deg"
    for key, value in dict(period=P, cu_t=CU, rt_t=RT, bias_t=BIAS,
                           air_top=40.035, air_bottom=-43.499).items():
        h[key] = f"{value:.12g}mm"

    def box(name, pos, size, material="pec"):
        return need(m.create_box(pos, size, name=name, material=material), "CREATE_" + name)

    diode = state in VALUES
    rotated = state not in ("SHORT_UNROT", "OFF_D_UNROT")
    copper_audit = []
    boundaries = []
    for n in range(8):
        angle = n * 22.5 if rotated else 0
        dx = (n - 3.5) * P
        disc = need(m.create_cylinder("Z", ["0mm", "0mm", "0mm"],
                                      f"{R}mm", "cu_t", name=f"DISC_{n}",
                                      material="pec"), "DISC_CREATE")
        # Two closed 4.11-mm slots leave a 0.5-mm copper bridge at y=0.
        cuts = [box(f"SLOT_BOT_{n}", [-SLOT_W/2, -SLOT_L/2, "-cu_t"],
                    [SLOT_W, (SLOT_L-STRIP_Y)/2, "3*cu_t"], "vacuum"),
                box(f"SLOT_TOP_{n}", [-SLOT_W/2, STRIP_Y/2, "-cu_t"],
                    [SLOT_W, (SLOT_L-STRIP_Y)/2, "3*cu_t"], "vacuum")]
        if diode:
            cuts.append(box(f"DIODE_GAP_{n}", [-GAP_X/2, -STRIP_Y/2, "-cu_t"],
                            [GAP_X, STRIP_Y, "3*cu_t"], "vacuum"))
        need(m.subtract(disc.name, [c.name for c in cuts], keep_originals=False),
             f"COPPER_SUBTRACT_{n}")
        expected_area = math.pi * R**2 - 2*4.11*SLOT_W - (GAP_X*STRIP_Y if diode else 0)
        area = float(disc.volume) / CU
        need(abs(area-expected_area) < 0.08, f"COPPER_AREA_{n}: {area} vs {expected_area}")
        need(disc.rotate("Z", angle), f"ROTATE_DISC_{n}")
        need(disc.move([dx, 0, 0]), f"MOVE_DISC_{n}")
        copper_audit.append({"cell": n, "rotation_deg": angle, "disc_area_mm2": area,
                             "disc_volume_mm3": float(disc.volume), "centre_x_mm": dx})
        if diode:
            # Sheet's x-directed current line bridges the full 0.3-mm cut.
            sheet = need(m.create_rectangle("XY", [-GAP_X/2, -STRIP_Y/2, CU],
                                            [GAP_X, STRIP_Y], name=f"DIODE_{n}"),
                         f"DIODE_SHEET_{n}")
            need(sheet.rotate("Z", angle), f"ROTATE_DIODE_{n}")
            need(sheet.move([dx, 0, 0]), f"MOVE_DIODE_{n}")
            ohms, x_ohms = VALUES["OFF_D" if state=="OFF_D_UNROT" else state]
            if element == "RLC":
                # Series R + C/L chosen to reproduce the stated X at 5.68 GHz.
                cap = -1 / (2*math.pi*5.68e9*x_ohms) if x_ohms < 0 else None
                ind = x_ohms / (2*math.pi*5.68e9) if x_ohms > 0 else None
                start = [dx-GAP_X/2*math.cos(math.radians(angle)),
                         -GAP_X/2*math.sin(math.radians(angle)), CU]
                end = [dx+GAP_X/2*math.cos(math.radians(angle)),
                       GAP_X/2*math.sin(math.radians(angle)), CU]
                bound = need(h.assign_lumped_rlc_to_sheet(sheet.name, [start, end],
                                                          name=f"RLC_{n}", rlc_type="Serial",
                                                          resistance=ohms, capacitance=cap,
                                                          inductance=ind), f"RLC_ASSIGN_{n}")
            else:
                scale = STRIP_Y/GAP_X
                bound = need(h.assign_impedance_to_sheet(sheet.name, name=f"IMP_{n}",
                                                          resistance=ohms*scale,
                                                          reactance=x_ohms*scale),
                             f"IMP_ASSIGN_{n}")
            boundaries.append({"name":bound.name,"props":dict(bound.props),
                               "sheet":sheet.name,"expected_R_ohm":ohms,
                               "expected_X_ohm":x_ohms,
                               "expected_Zs_ohm_sq":[ohms*STRIP_Y/GAP_X,
                                                       x_ohms*STRIP_Y/GAP_X]})

    box("GROUND", ["-4*period", "-period/2", "-rt_t-cu_t"],
        ["8*period", "period", "cu_t"])
    substrate = box("SUBSTRATE", ["-4*period", "-period/2", "-rt_t"],
                    ["8*period", "period", "rt_t"], "RT5880_SUPERCELL")
    box("BIAS_RT5880", ["-4*period", "-period/2", "-rt_t-cu_t-bias_t"],
        ["8*period", "period", "bias_t"], "RT5880_SUPERCELL")
    box("AIR", ["-4*period", "-period/2", "air_bottom"],
        ["8*period", "period", "air_top-air_bottom"], "vacuum")
    face = lambda xyz: need(m.get_faceid_from_position(xyz, assignment="AIR"), "AIR_FACE")
    zm = (40.035 - 43.499)/2
    need(h.assign_lattice_pair([face([-4*P,0,zm]),face([4*P,0,zm])],
                               phase_delay="InputPhaseDelay", phase_delay_param1="$ph_s",
                               name="Lattice_X"), "LATTICE_X")
    need(h.assign_lattice_pair([face([0,-P/2,zm]),face([0,P/2,zm])],
                               phase_delay="InputPhaseDelay", phase_delay_param1="$ph_y",
                               name="Lattice_Y"), "LATTICE_Y")
    h.oboundary.CreatePML(["NAME:PMLCreationSettings", "UserDrawnGroup:=", False,
                           "PMLFaces:=", [face([0,0,40.035])],
                           "CreateJoiningObjs:=", False, "Thickness:=", "10mm",
                           "RadDist:=", "20mm", "UseFreq:=", True, "MinFreq:=", "5.8GHz"])
    need(h.assign_perfect_e([face([0,0,-43.499])], name="PEC_BOTTOM"), "PEC_BOTTOM")
    setup = need(h.create_setup("EigenCell"), "SETUP")
    setup.props.update(MinimumFrequency="5.3GHz", NumModes=20, MaxDeltaFreq=.2,
                       MinimumPasses=4, MinimumConvergedPasses=2,
                       MaximumPasses=20, ConvergeOnRealFreq=True, SaveFields=True)
    need(setup.update(), "SETUP_UPDATE")
    validation_file = OUT / "validation" / f"{state}_{element}.log"
    validation_file.parent.mkdir(parents=True, exist_ok=True)
    need(h.validate_simple(validation_file) == 1, "DESIGN_VALIDATION")
    bnames = [b.name for b in h.boundaries]
    need(all(x in bnames for x in ("Lattice_X", "Lattice_Y", "PEC_BOTTOM")),
         "BOUNDARY_READBACK")
    need(any("pml" in x.lower() for x in m.object_names), "PML_READBACK")
    need(len(boundaries) == (8 if diode else 0), "DIODE_BOUNDARY_COUNT")
    need(substrate.material_name.lower() == "rt5880_supercell", "SUBSTRATE_READBACK")
    need(abs(float(substrate.bounding_box[5])-float(substrate.bounding_box[2])-RT)<1e-6,
         "THICKNESS_READBACK")
    need(h.save_project(), "SAVE_PROJECT")
    return {"state":state,"element":element if diode else "NONE", "phase_deg":phase,
            "copper":copper_audit,"diode_boundaries":boundaries,
            "aedt_boundaries":{b.name:dict(b.props) for b in h.boundaries},
            "setup":dict(setup.props),"aedt_version":h.desktop_class.aedt_version_id,
            "substrate":"RT5880_SUPERCELL 3.175mm","bias":"RT5880 0.254mm",
            "PML_top_mm":40.035,"PML_bottom_PEC_mm":-43.499,
            "validation":validation_file.read_text(errors="replace")}


def run(state, element, cores, attempt):
    OUT.mkdir(parents=True, exist_ok=True)
    label=f"{state}_{element}_attempt{attempt}"
    project=WORK/f"{label}.aedt"
    need(not project.exists(), f"REFUSE_EXISTING_PROJECT {project}")
    temp=WORK/"temp"/label
    temp.mkdir(parents=True,exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    status_path=OUT/f"{label}_status.json"
    s={"state":state,"element":element,"attempt":attempt,"status":"STARTING",
       "controller_pid":os.getpid(),"started_utc":utc(),"solve_timeout_s":3600}
    save(status_path,s)
    stop=threading.Event()
    def beat():
        while not stop.wait(60):
            s["heartbeat_utc"]=utc();save(status_path,s)
    threading.Thread(target=beat,daemon=True).start()
    h=None
    try:
        h=Hfss(project=str(project),design="EigenCell",solution_type="Eigenmode",
               version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
        h.desktop_class.temp_directory=temp
        need(str(h.desktop_class.aedt_version_id)=="2026.1", "AEDT_VERSION_READBACK")
        need(Path(h.desktop_class.temp_directory).resolve()==temp.resolve(), "PRIVATE_TEMP_READBACK")
        s.update(status="BUILDING",aedt_pid=h.desktop_class.aedt_process_id)
        save(status_path,s)
        built=build(h,state,element=element)
        save(OUT/f"{label}_build.json",built)
        s.update(status="SOLVING",solve_start_utc=utc(),cores_requested=cores)
        save(status_path,s)
        watchdog_log=(OUT/f"{label}_watchdog.log").open("w",encoding="utf-8")
        watcher=subprocess.Popen([sys.executable,str(Path(__file__).with_name("sub6_external_watchdog.py")),
                                  str(status_path)],cwd=str(ROOT),stdout=watchdog_log,
                                 stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
        s["external_watchdog_pid"]=watcher.pid;save(status_path,s)
        need(h.analyze_setup("EigenCell",cores=cores,use_auto_settings=True,blocking=True),
             "SOLVER_FAILURE")
        need(s["status"]=="SOLVING", "SOLVE_INTERRUPTED")
        result={"state":state,"element":element,"attempt":attempt,"status":"COMPLETE",
                "aedt_version":"2026.1","modes":modes(h),
                "solve_end_utc":utc(),"aedt_pid":h.desktop_class.aedt_process_id,
                "project":str(project)}
        result["profile_path"]=h.export_profile("EigenCell",output_file=str(OUT/f"{label}_profile.prof"))
        result["convergence_path"]=h.export_convergence("EigenCell",output_file=str(OUT/f"{label}_convergence.conv"))
        need(h.save_project(),"POSTSOLVE_SAVE_PROJECT")
        save(OUT/f"{label}_modes.json",result)
        s.update(status="COMPLETE",solve_end_utc=utc(),modes=len(result["modes"]))
        save(status_path,s)
        print(json.dumps({"state":state,"modes":len(result["modes"]),
                          "highest_ghz":result["modes"][-1]["f_ghz"]}),flush=True)
    except BaseException as exc:
        s.update(status="FAILED",error=f"{type(exc).__name__}: {exc}",traceback=traceback.format_exc())
        save(status_path,s)
        traceback.print_exc()
        raise
    finally:
        stop.set()
        if h is not None:
            own_pid=int(h.desktop_class.aedt_process_id)
            try:
                own_process=psutil.Process(own_pid)
                own_created=own_process.create_time()
            except psutil.NoSuchProcess:
                own_process=None
            try:
                h.release_desktop(close_projects=True,close_desktop=True)
            finally:
                if own_process is not None:
                    try:
                        own_process.wait(timeout=10)
                    except psutil.TimeoutExpired:
                        if own_process.create_time()==own_created:
                            s["terminated_own_aedt_after_release"]=own_pid
                            save(status_path,s)
                            own_process.terminate()
                    except psutil.NoSuchProcess:
                        pass


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--state",choices=("ON_A","ON_B","SHORT","SHORT_UNROT","OFF_D","OFF_D_UNROT"),required=True)
    ap.add_argument("--element",choices=("RLC","IMP"),default="RLC")
    ap.add_argument("--cores",type=int,default=46)
    ap.add_argument("--attempt",type=int,default=1)
    a=ap.parse_args()
    run(a.state,a.element,a.cores,a.attempt)
