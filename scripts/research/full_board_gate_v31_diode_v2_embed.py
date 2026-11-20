"""Bake corrected DIODE V2 CAD geometry into a results-free archive."""

import json
import os
from pathlib import Path
import threading
import traceback
import zipfile

import psutil

os.environ["ANSYSEM_ROOT261"]=r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/"hfss_work/full_board_gate_v31_diode_v2"
PROJECT=WORK/"BOARD_R150_V31_B_ALL_DIODE_V2.aedt"
OUT=ROOT/"results/sub6_meander_model/full_board_gate/v31_diode_v2"
ARCHIVE=OUT/"BOARD_R150_V31_B_ALL_DIODE_V2.aedtz"
STATUS=OUT/"embed_status.json"


def save(data):
    STATUS.write_text(json.dumps(data,indent=2,default=str)+"\n",encoding="utf-8")


def need(value,label):
    if not value: raise RuntimeError(label)
    return value


def main():
    need(PROJECT.is_file() and not ARCHIVE.exists(),"PROJECT_OR_ARCHIVE_STATE")
    temp=WORK/"embed_temp";temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    state={"phase":"START","controller_pid":os.getpid()};save(state)
    h=None
    try:
        h=Hfss(project=str(PROJECT),design="GATE_B_ALL_DIODE_V2",solution_type="Modal",
               version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
        state["aedt_pid"]=h.desktop_class.aedt_process_id
        need(str(h.desktop_class.aedt_version_id)=="2026.1","VERSION")
        h.desktop_class.temp_directory=temp
        need(Path(h.desktop_class.temp_directory).resolve()==temp.resolve(),"TEMP")
        m=h.modeler
        metal=[n for n in m.object_names if n.startswith("METAL_PEC")]
        sheets=[n for n in m.object_names if n.startswith("DIODE_")]
        need(len(metal)==518 and len(sheets)==250,"OBJECT_COUNTS_BEFORE")
        audit_path=OUT/"audit.json"
        audit=json.loads(audit_path.read_text(encoding="utf-8"))
        need(audit.get("all_gap_probes_pass") and len(audit.get("cuts",[]))==250,
             "PROBE_AUDIT_INCOMPLETE")
        boundary=next(b for b in h.boundaries if b.name=="DIODE")
        need(len(boundary.props.get("Objects",[]))==250,"BOUNDARY_BEFORE")
        setup=h.get_setup("Setup1")
        setup_before={k:str(setup.props.get(k)) for k in
                      ("Frequency","MaximumPasses","MinimumPasses","MaxDeltaS","BasisOrder","DrivenSolverType")}
        b_before={b.name:dict(b.props) for b in h.boundaries}
        sphere_before=[v.name for v in h.field_setups]
        volume_before=sum(float(m[n].volume) for n in metal)
        imported=metal+["SUBSTRATE","BIAS_RT5880","AIR","PORT"]
        need(all(n in m.object_names for n in imported),"IMPORTED_OBJECT_MISSING")
        state.update(phase="PURGING_CAD_HISTORY",imported_object_count=len(imported));save(state)
        need(m.purge_history(imported),"PURGE_HISTORY")
        need(len([n for n in m.object_names if n.startswith("METAL_PEC")])==518,"COPPER_COUNT_AFTER")
        need(len([n for n in m.object_names if n.startswith("DIODE_")])==250,"SHEET_COUNT_AFTER")
        volume_after=sum(float(m[n].volume) for n in metal)
        need(abs(volume_after-volume_before)<1e-3,"COPPER_VOLUME_CHANGED")
        need({k:str(h.get_setup("Setup1").props.get(k)) for k in setup_before}==setup_before,"SETUP_CHANGED")
        need({b.name:dict(b.props) for b in h.boundaries}==b_before,"BOUNDARIES_CHANGED")
        need([v.name for v in h.field_setups]==sphere_before,"SPHERE_CHANGED")
        valid=h.validate_simple(str(OUT/"validation_embedded.log"))
        need(valid==1 or valid is True,"VALIDATION")
        need(h.save_project(),"SAVE")
        step_tokens=PROJECT.read_bytes().lower().count(b".step")
        need(step_tokens==0,f"EXTERNAL_STEP_REFERENCES_REMAIN_{step_tokens}")
        need(h.archive_project(ARCHIVE,include_external_files=False,include_results_file=False),"ARCHIVE")
        with zipfile.ZipFile(ARCHIVE) as z:
            entries=[(x.filename,x.file_size) for x in z.infolist()]
            aedt=next(x.filename for x in z.infolist() if x.filename.lower().endswith(".aedt"))
            archived_step_tokens=z.read(aedt).lower().count(b".step")
        need(archived_step_tokens==0,"ARCHIVED_STEP_REFERENCES")
        need(not any(n.lower().endswith(".step") or "aedtresults" in n.lower()
                     for n,_ in entries),"ARCHIVE_EXTERNAL_OR_RESULTS")
        audit.update(project_bytes=PROJECT.stat().st_size,archive_bytes=ARCHIVE.stat().st_size,
                     archive_members=[n for n,_ in entries],
                     cad_history_purged_objects=len(imported),
                     external_step_references=step_tokens,
                     archived_step_references=archived_step_tokens,
                     final_metal_objects=518,diode_sheet_count=250,
                     volume_before_after_mm3=[volume_before,volume_after],
                     validation_embedded=valid)
        audit_path.write_text(json.dumps(audit,indent=2)+"\n",encoding="utf-8")
        state.update(phase="COMPLETE",project_bytes=PROJECT.stat().st_size,
                     archive_bytes=ARCHIVE.stat().st_size,step_tokens=step_tokens,
                     validation=valid)
        save(state)
        (OUT/"embed.done").write_text("COMPLETE\n",encoding="utf-8")
    except BaseException:
        state.update(phase="FAILED",traceback=traceback.format_exc());save(state)
        (OUT/"embed.done").write_text("FAILED\n",encoding="utf-8")
        raise
    finally:
        if h is not None:
            pid=int(h.desktop_class.aedt_process_id)
            try: own=psutil.Process(pid);created=own.create_time()
            except psutil.NoSuchProcess: own=None
            h.release_desktop(close_projects=True,close_desktop=True)
            if own is not None:
                try: own.wait(timeout=10)
                except psutil.TimeoutExpired:
                    if own.create_time()==created: own.terminate()
                except psutil.NoSuchProcess: pass


if __name__=="__main__": main()
