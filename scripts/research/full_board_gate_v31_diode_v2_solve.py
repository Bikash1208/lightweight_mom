"""One attached HFSS 2026 R1 solve of corrected DIODE V2; own AEDT PID only."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime,timezone

import psutil

os.environ["ANSYSEM_ROOT261"]=r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"hfss_work/full_board_gate_v31_diode_v2/BOARD_R150_V31_B_ALL_DIODE_V2.aedt"
OUT=ROOT/"results/sub6_meander_model/full_board_gate/v31_diode_v2"
WORK=ROOT/"hfss_work/full_board_gate_v31_diode_v2"
VALUES={"ON_A":("3.333ohm","-232.667ohm"),"OFF_D":("2.667ohm","48.667ohm")}


def utc():return datetime.now(timezone.utc).isoformat()


def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(obj,indent=2,default=str)+"\n",encoding="utf-8")
    os.replace(tmp,path)


def need(ok,label):
    if not ok:raise RuntimeError(label)
    return ok


def parse_profile(path):
    t=path.read_text(errors="replace")
    deltas=[float(x) for x in re.findall(r"Max Mag\. Delta S:\s*([\d.]+)",t)]
    passes=len(re.findall(r"^\s*Adaptive Pass \d+\s+Frequency:",t,re.M))
    memory=[float(x) for x in re.findall(r"Max memory/process:\s*([\d.]+)\s*GB",t)]
    tets=[int(x) for x in re.findall(r"Tetrahedra:\s*(\d+)",t)]
    cores=[int(x) for x in re.findall(r"Cores:\s*(\d+)",t)]
    return {"adaptive_passes":passes,"delta_s_per_pass":[None]+deltas,
            "final_delta_s":deltas[-1] if deltas else None,
            "peak_memory_gb":max(memory) if memory else None,
            "max_tets":max(tets) if tets else None,
            "cores_granted":max(cores) if cores else None,
            "converged":t.rfind("Adaptive Passes converged")>t.rfind("Adaptive Passes did not converge"),
            "normal_completion":"Status: Normal Completion" in t}


def close(h):
    if h is None:return
    pid=int(h.desktop_class.aedt_process_id)
    try:p=psutil.Process(pid);created=p.create_time()
    except psutil.NoSuchProcess:p=None
    h.release_desktop(close_projects=True,close_desktop=True)
    if p is not None:
        try:p.wait(timeout=10)
        except psutil.TimeoutExpired:
            if p.create_time()==created:p.terminate()
        except psutil.NoSuchProcess:pass


def run(state):
    need(state in VALUES,"STATE")
    need(SOURCE.is_file(),"V2_PROJECT_MISSING")
    solve_dir=WORK/f"solve_{state}_001"
    project=solve_dir/f"BOARD_R150_V31_B_ALL_DIODE_V2_{state}.aedt"
    need(not project.exists(),"REFUSE_EXISTING_SOLVE_PROJECT")
    solve_dir.mkdir(parents=True,exist_ok=True)
    shutil.copy2(SOURCE,project)
    temp=solve_dir/"aedt_temp";temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    status_path=OUT/f"{state}_solve_status.json"
    status={"status":"STARTING","state":state,"controller_pid":os.getpid(),
            "started_utc":utc(),"project":str(project),"solve_timeout_s":14400}
    save(status_path,status)
    stop=threading.Event()
    def beat():
        while not stop.wait(60):
            status["heartbeat_utc"]=utc()
            try:
                p=psutil.Process(int(status["aedt_pid"]))
                tree=[p]+p.children(recursive=True)
                status["own_process_tree_gib"]=sum(q.memory_info().rss for q in tree if q.is_running())/2**30
                status["free_memory_gib"]=psutil.virtual_memory().available/2**30
            except (KeyError,psutil.Error):pass
            save(status_path,status)
    threading.Thread(target=beat,daemon=True).start()
    h=None
    try:
        h=Hfss(project=str(project),design="GATE_B_ALL_DIODE_V2",solution_type="Modal",
               version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
        status["aedt_pid"]=h.desktop_class.aedt_process_id
        status["aedt_version"]=str(h.desktop_class.aedt_version_id)
        need(status["aedt_version"]=="2026.1","AEDT_VERSION")
        h.desktop_class.temp_directory=temp
        need(Path(h.desktop_class.temp_directory).resolve()==temp.resolve(),"PRIVATE_TEMP")
        need(len(h.setups)==1 and h.setups[0].name=="Setup1","SETUP_COUNT")
        setup=h.get_setup("Setup1")
        settings={k:str(setup.props.get(k)) for k in
                  ("Frequency","MaximumPasses","MinimumPasses","MaxDeltaS","BasisOrder","DrivenSolverType")}
        need(settings["Frequency"]=="5.68GHz" and settings["MaximumPasses"]=="12"
             and settings["MinimumPasses"]=="6" and settings["MaxDeltaS"]=="0.02"
             and "direct" in settings["DrivenSolverType"].lower(),"SETUP_READBACK")
        diode=next((b for b in h.boundaries if b.name=="DIODE"),None)
        need(diode is not None and len(diode.props.get("Objects",[]))==250,"DIODE_BOUNDARY")
        expected=VALUES[state]
        if state=="OFF_D":
            h["Zr"],h["Zx"]=expected
        need(str(h["Zr"]).lower()==expected[0].lower() and
             str(h["Zx"]).lower()==expected[1].lower(),"DIODE_VARIABLES")
        need(len([n for n in h.modeler.object_names if n.startswith("METAL_PEC")])==518,
             "METAL_COUNT")
        status.update(status="VALIDATING",setup=settings,variables={"Zr":str(h["Zr"]),"Zx":str(h["Zx"])},
                      diode_sheets=250)
        save(status_path,status)
        valid=h.validate_simple(str(OUT/f"{state}_validation.log"))
        need(valid==1 or valid is True,"VALIDATION")
        need(h.save_project(),"SAVE_PRESOLVE")
        busy=psutil.cpu_percent(interval=2)
        logical=psutil.cpu_count(logical=True) or 48
        cores=max(1,min(logical-2,round(logical*(100-busy)/100)-2))
        status.update(status="SOLVING",solve_start_utc=utc(),cores_requested=cores,
                      free_memory_gib_before=psutil.virtual_memory().available/2**30)
        save(status_path,status)
        watchlog=(OUT/f"{state}_external_watchdog.log").open("w",encoding="utf-8")
        watcher=subprocess.Popen([sys.executable,str(Path(__file__).with_name("sub6_external_watchdog.py")),
                                  str(status_path)],cwd=str(ROOT),stdout=watchlog,
                                 stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
        status["external_watchdog_pid"]=watcher.pid;save(status_path,status)
        t0=time.monotonic()
        need(h.analyze_setup("Setup1",cores=cores,use_auto_settings=False,blocking=True),
             "SOLVER_FAILURE")
        status["solve_wall_s"]=time.monotonic()-t0
        prof=Path(h.export_profile("Setup1",output_file=str(OUT/f"{state}_solve.prof")))
        evidence=parse_profile(prof)
        status.update(profile=evidence,solve_end_utc=utc())
        need(evidence["normal_completion"] and evidence["converged"],"SOLVE_NOT_CONVERGED")
        need(6<=evidence["adaptive_passes"]<=12 and evidence["final_delta_s"]<=.02,
             "CONVERGENCE_GATE")
        need(h.save_project(),"POSTSOLVE_SAVE")
        status["status"]="COMPLETE";save(status_path,status)
        print("SOLVED",state,evidence,flush=True)
    except BaseException:
        status.update(status="FAILED",traceback=traceback.format_exc(),failed_utc=utc())
        save(status_path,status)
        traceback.print_exc()
        raise
    finally:
        stop.set();close(h)


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--state",choices=VALUES,required=True)
    run(ap.parse_args().state)
