"""Unattended, one-at-a-time C4X Steps 0-2 queue; own AEDT PIDs only."""

import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone

from sub6_c4x_supercell1_analyze import R, field, tag as solved_tag, run as analyze

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=Path(__file__).resolve().parent
STATUS=R/"C4X_SUPERCELL1_queue_status.json"
FINAL=R/"C4X_SUPERCELL1.done"


def stamp(): return datetime.now(timezone.utc).isoformat()


def save(value):
    temp=STATUS.with_suffix(".json.tmp")
    temp.write_text(json.dumps(value,indent=2,default=str)+"\n",encoding="utf-8")
    os.replace(temp,STATUS)


def done(label,detail):
    (R/f"C4X_{label}.done").write_text(json.dumps({"utc":stamp(),**detail},indent=2)+"\n",encoding="utf-8")


def row(cell,phase,state):
    tag="c4xrun" if state=="ON" and phase==15 else "c4xauto1"
    return R/f"SUPER_C4X_{state}_ph{phase}_PML40_{tag}_modes.json"


def command(args,label,state):
    logfile=R/f"C4X_queue_{label}.log"
    with logfile.open("w",encoding="utf-8") as log:
        proc=subprocess.Popen([sys.executable,*map(str,args)],cwd=str(ROOT),stdout=log,stderr=subprocess.STDOUT,
                              creationflags=subprocess.CREATE_NO_WINDOW)
        state["current"]={"label":label,"controller_pid":proc.pid,"start_utc":stamp(),"log":str(logfile)}
        save(state)
        while proc.poll() is None:
            # AEDT PID is recorded only from this row's own status/readback file.
            if label.startswith("SUPER_"):
                parts=label.split("_")
                status=R/f"SUPER_C4X_{parts[1]}_ph{parts[2]}_PML40_{parts[3]}_status.json"
                if status.exists():
                    try:
                        own=json.loads(status.read_text())
                        state["current"]["aedt_pid"]=own.get("aedt_pid")
                    except (OSError,ValueError): pass
            state["heartbeat_utc"]=stamp(); save(state)
            time.sleep(60)
        state["current"]["exit_code"]=proc.returncode
        state["current"]["end_utc"]=stamp();save(state)
        return proc.returncode


def solve(state,kind,phase):
    label=f"SUPER_{kind}_{phase}_c4xauto1"
    if row(kind,phase,kind).exists():
        done(label,{"status":"COMPLETE_EXISTING"}); return True
    for attempt in (1,2):
        tag=f"c4xauto{attempt}"
        args=[SCRIPTS/"sub6_supercell_phase.py","--cell",f"C4X_{kind}","--phase",str(phase),
              "--cores","48","--tag",tag,"--min-freq-ghz","5.3","--num-modes","20",
              "--timeout-s","3600"]
        actual=f"SUPER_{kind}_{phase}_{tag}"
        code=command(args,actual,state)
        status_path=R/f"SUPER_C4X_{kind}_ph{phase}_PML40_{tag}_status.json"
        status=json.loads(status_path.read_text()) if status_path.exists() else {}
        state["rows"].append({"label":actual,"exit":code,"status":status.get("status"),
                              "aedt_pid":status.get("aedt_pid"),"controller_pid":status.get("controller_pid")})
        save(state)
        if code==0 and status.get("status")=="COMPLETE":
            done(actual,{"status":"COMPLETE","aedt_pid":status.get("aedt_pid")});return True
        # Only a timeout or launch failure gets the one retry.
        launch_fail=not status.get("solve_start_utc")
        timeout=status.get("status")=="TIMEOUT" or "TIMEOUT" in str(status.get("error",""))
        if not (timeout or launch_fail):
            raise RuntimeError(f"SOLVER_OR_READBACK_FAILURE {actual}: {status.get('error')}")
        done(actual,{"status":"RETRY" if attempt==1 else "ABORTED","reason":status.get("error")})
    return False


def export(state,kind,phase):
    tag=solved_tag(kind,phase)
    label=f"FIELDS_{kind}_{phase}"
    args=[SCRIPTS/"sub6_export_supercell_fields.py","--cell",f"C4X_{kind}","--phase",str(phase),"--tag",tag]
    code=command(args,label,state)
    meta=R/f"SUPER_C4X_{kind}_ph{phase}_PML40_{tag}_field_exports.json"
    info=json.loads(meta.read_text()) if meta.exists() else {}
    state["rows"].append({"label":label,"exit":code,"aedt_pid":info.get("aedt_pid")});save(state)
    if code: raise RuntimeError(f"FIELD_EXPORT_FAILURE {label}")
    done(label,{"status":"COMPLETE","aedt_pid":info.get("aedt_pid")})


def preflight(state):
    refs={("ON","across"):5.730,("ON","along"):5.578,
          ("OFF","across"):5.609,("OFF","along"):5.610}
    spur=[];cal=[]
    for (kind,direction),reference in refs.items():
        prefix=f"C4X_{kind}"
        pec=json.loads((R/f"{prefix}_PEC_{direction}_85_air40.035_c4xcal_modes.json").read_text())["modes"][0]
        pct=100*(pec["f_ghz"]-reference)/reference
        cal.append({"state":kind,"direction":direction,"f_ghz":pec["f_ghz"],"reference_ghz":reference,
                    "error_pct":pct,"Q":pec["Q"]})
        if abs(pct)>.5: raise RuntimeError(f"CALIBRATION_MISS {kind} {direction} {pct}")
        pml=[]
        for phase in (85,95):
            pml.append(json.loads((R/f"{prefix}_PML_{direction}_{phase}_air40.035_c4xbase_modes.json").read_text())["modes"][0])
        vg=abs((pml[1]["f_ghz"]-pml[0]["f_ghz"])*1e9/10)*360*.0116
        alpha=math.pi*pml[0]["f_ghz"]*1e9*(1/pml[0]["Q"]-1/pec["Q"])/vg
        spur.append({"state":kind,"direction":direction,"alpha_spur_np_m":alpha,
                     "Q_PEC":pec["Q"],"Q_PML":pml[0]["Q"],"vg_m_s":vg})
        if alpha>.1: raise RuntimeError(f"PML_SPURIOUS_LOSS {kind} {direction} {alpha}")
        done(f"CAL_{kind}_{direction}",{"status":"PASS","frequency_error_pct":pct})
        done(f"BASE_{kind}_{direction}",{"status":"PASS","alpha_spur_np_m":alpha})
    state["calibration"]=cal;state["spurious_alpha"]=spur;save(state)


def write_report(state,analysis=None):
    rows=[]
    path=R/"C4X_SUPERCELL1_rows.csv"
    if path.exists():
        with path.open(encoding="utf-8") as fh: rows=list(csv.DictReader(fh))
    candidates=(analysis or {}).get("on_gate_candidates",[])
    lines=["# C4X SUPERCELL 1", "",f"STATUS: {state['status']}","",
           "GATES: Step 0 calibration and Step 1 PML spurious loss PASS." if state.get("calibration") else
           "GATES: calibration or PML gate incomplete.",
           f"ON single-mode alpha gate: {'PASS' if candidates else 'NOT DEMONSTRATED'}; OFF common-frequency gate: {state.get('off_gate','NOT RUN')}.","",
           "DATA: HFSS 2025.2, 40 mm air/PML, 20 modes from 5.3 GHz, 3600 s row timeout.",
           f"Mode inventory: {len(rows)} rows in `C4X_SUPERCELL1_rows.csv`."]
    for x in state.get("calibration",[]):
        lines.append(f"- PEC {x['state']} {x['direction']} 85°: {x['f_ghz']:.6f} GHz, Q {x['Q']:.1f}, reference error {x['error_pct']:+.4f}%.")
    for x in state.get("spurious_alpha",[]):
        lines.append(f"- PML spurious α {x['state']} {x['direction']}: {x['alpha_spur_np_m']:.5f} Np/m.")
    for x in candidates:
        lines.append(f"- ON candidate φ={x['phase_deg']}° mode {x['mode']}: f={x['f_ghz']:.6f} GHz, Q={x['Q']:.1f}, α={x['alpha_np_m']:.4f} Np/m, field cosine={x['field_cosine']:.3f}.")
    if state.get("off_match"):
        x=state["off_match"]
        lines.append(f"- Common-frequency match: ON {x['on']['f_ghz']:.6f} GHz, α={x['on']['alpha_np_m']:.4f}; OFF track {x['off']['track']} α={x['off']['alpha_np_m']:.4f} Np/m, field cosine={x['off']['field_cosine_to_ON_at_15']}.")
    for x in rows:
        if x["cell"]=="OFF" and (abs(float(x["f_ghz"])-5.61)<.08 or abs(float(x["f_ghz"])-6.06)<.08):
            lines.append(f"- OFF near symmetric resonance: φ={x['phase_deg']}° mode {x['mode']}, f={float(x['f_ghz']):.6f} GHz, Q={float(x['Q']):.1f}, α diagnostic={x['alpha_np_m']} Np/m.")
    lines += ["","FLAGS: Field matching uses 16 sampled |E| points and frequency continuity, not full complex overlap. Both Bloch folds in the CSV are kinematic for low-Q and box modes. A supercell mode's leakage does not establish pin-launch coupling.",
              "Only own fresh local AEDT projects and exact own PIDs were used; private project files were not pushed.",
              "", "FILES: `C4X_SUPERCELL1_rows.csv`, `C4X_SUPERCELL1_analysis.json`, raw row JSON, convergence/profile files, queue status and logs.",
              "",f"HANDOFF: {state['status']} | C4X ideal ON/OFF | ON α gate {'PASS' if candidates else 'NOT DEMONSTRATED'} | OFF {state.get('off_gate','NOT RUN')} | HFSS 2025.2"]
    (R/"C4X_SUPERCELL1.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


def main():
    R.mkdir(parents=True,exist_ok=True)
    state={"status":"RUNNING","chain_pid":os.getpid(),"started_utc":stamp(),"rows":[]}
    save(state)
    try:
        preflight(state)
        done("SUPER_ON_15_c4xrun",{"status":"COMPLETE_EXISTING"})
        incomplete=[]
        for phase in (30,45):
            if not solve(state,"ON",phase): incomplete.append(f"ON_{phase}")
        if incomplete:
            state["status"]="ABORTED";state["reason"]="ON rows incomplete: "+",".join(incomplete)
            return
        for phase in (15,30,45): export(state,"ON",phase)
        analysis=analyze(("ON",))
        state["on_gate_candidates"]=len(analysis["on_gate_candidates"]);save(state)
        if not analysis["on_pass"]:
            state["status"]="STOPPED_ON_GATE";return
        incomplete=[]
        for phase in (15,30):
            if not solve(state,"OFF",phase): incomplete.append(f"OFF_{phase}")
        if incomplete:
            state["status"]="ABORTED";state["reason"]="OFF rows incomplete: "+",".join(incomplete)
            return
        for phase in (15,30): export(state,"OFF",phase)
        analysis=analyze(("ON","OFF"))
        # A common-frequency OFF match needs an OFF track bracketing the chosen ON frequency.
        best=max(analysis["on_gate_candidates"],key=lambda x:x["alpha_np_m"])
        offrows=list(csv.DictReader((R/"C4X_SUPERCELL1_rows.csv").open(encoding="utf-8")))
        trackids={r["track_id"] for r in offrows if r["cell"]=="OFF" and r["track_id"]}
        matches=[]
        for trackid in trackids:
            pair=sorted([r for r in offrows if r["cell"]=="OFF" and r["track_id"]==trackid],key=lambda x:float(x["phase_deg"]))
            if len(pair)!=2: continue
            fa,fb=[float(r["f_ghz"]) for r in pair]
            if min(fa,fb)<=best["f_ghz"]<=max(fa,fb):
                t=(best["f_ghz"]-fa)/(fb-fa)
                aa,ab=[float(r["alpha_np_m"]) for r in pair]
                matches.append({"track":trackid,"alpha_np_m":aa+t*(ab-aa),
                                "field_cosine_to_ON_at_15":None})
        if matches:
            # Same-phase field profile chooses the OFF continuation closest to ON.
            compare_phase=best["phase_deg"] if best["phase_deg"] in (15,30) else 30
            onfield=field("ON",compare_phase,best["mode"])
            for item in matches:
                same=next(r for r in offrows if r["cell"]=="OFF" and r["track_id"]==item["track"] and r["phase_deg"]==str(compare_phase))
                offield=field("OFF",compare_phase,int(same["mode"]))
                item["field_cosine_to_ON_at_15"]=float(onfield@offield) if onfield is not None and offield is not None else None
            chosen=max(matches,key=lambda x:x["field_cosine_to_ON_at_15"] or -1)
            state["off_gate"]="PASS" if chosen["alpha_np_m"]<=best["alpha_np_m"]/12 else "FAIL"
            state["off_match"]={"on":best,"off":chosen,"alternatives":matches}
        else:
            state["off_gate"]="UNRESOLVED_NO_BRACKET"
        state["status"]="COMPLETE";save(state)
    except BaseException:
        state["status"]="STOPPED_BLOCKER";state["traceback"]=traceback.format_exc()
    finally:
        try:
            analysis=json.loads((R/"C4X_SUPERCELL1_analysis.json").read_text()) if (R/"C4X_SUPERCELL1_analysis.json").exists() else None
            write_report(state,analysis)
        except BaseException:
            state["report_traceback"]=traceback.format_exc()
        state["ended_utc"]=stamp();save(state)
        FINAL.write_text(json.dumps({"status":state["status"],"utc":stamp(),"chain_pid":os.getpid()})+"\n",encoding="utf-8")


if __name__=="__main__": main()
