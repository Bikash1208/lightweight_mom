"""One-at-a-time unattended diode supercell queue with own-PID records."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone

import psutil

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"results"/"sub6_meander_model"/"diode_supercells"
DRIVER=Path(__file__).with_name("sub6_diode_supercell.py")
ORDER=("ON_A","ON_B","SHORT","SHORT_UNROT","OFF_D","OFF_D_UNROT")
STATUS=OUT/"queue_status.json"
FINAL=OUT/"queue.done"


def utc(): return datetime.now(timezone.utc).isoformat()


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(data,indent=2,default=str)+"\n",encoding="utf-8")
    os.replace(temp,path)


def row(state,element,attempt,queue):
    tag=f"{state}_{element}_attempt{attempt}"
    log=OUT/f"{tag}_controller.log"
    sample=psutil.cpu_percent(interval=2)
    logical=psutil.cpu_count(logical=True) or 48
    cores=max(1,min(logical-2,round(logical*(100-sample)/100)-2))
    with log.open("w",encoding="utf-8") as out:
        proc=subprocess.Popen([sys.executable,str(DRIVER),"--state",state,
                               "--element",element,"--cores",str(cores),
                               "--attempt",str(attempt)],cwd=str(ROOT),
                              stdout=out,stderr=subprocess.STDOUT,
                              creationflags=subprocess.CREATE_NO_WINDOW)
        queue["current"]={"state":state,"element":element,"attempt":attempt,
                          "controller_pid":proc.pid,"cores_requested":cores,
                          "start_utc":utc(),"log":str(log)}
        write(STATUS,queue)
        while proc.poll() is None:
            source=OUT/f"{tag}_status.json"
            if source.exists():
                try:
                    own=json.loads(source.read_text(encoding="utf-8"))
                    queue["current"]["aedt_pid"]=own.get("aedt_pid")
                    queue["current"]["phase"]=own.get("status")
                except (ValueError,OSError): pass
            queue["heartbeat_utc"]=utc();write(STATUS,queue)
            time.sleep(60)
        queue["current"].update(exit_code=proc.returncode,end_utc=utc())
        source=OUT/f"{tag}_status.json"
        result=json.loads(source.read_text(encoding="utf-8")) if source.exists() else {}
        queue["current"]["status"]=result.get("status","NO_STATUS")
        write(STATUS,queue)
    return proc.returncode,result,tag


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    queue={"status":"RUNNING","controller_pid":os.getpid(),"started_utc":utc(),"rows":[]}
    write(STATUS,queue)
    for state in ORDER:
        success=False
        # RLC assignment failed before solve on ON_A. The second element type
        # is now the only allowed ON_A retry; later rows start with it.
        attempts=((2,"IMP"),) if state=="ON_A" else ((1,"IMP"),(2,"RLC"))
        for attempt,element in attempts:
            code,result,tag=row(state,element,attempt,queue)
            queue["rows"].append(dict(queue["current"]))
            write(STATUS,queue)
            if code==0 and result.get("status")=="COMPLETE":
                success=True
                (OUT/f"{state}.done").write_text(tag+"\n",encoding="utf-8")
                break
            if not result.get("solve_start_utc"):
                queue.update(status="STOP_BUILD_OR_VALIDATION",stop_state=state,
                             end_utc=utc(),reason=result.get("error"))
                write(STATUS,queue)
                write(FINAL,{"status":queue["status"],"utc":utc()})
                return
            # A failed AEDT solve must retain its own logs and retry once with the other element.
            ownproject=ROOT/"hfss_work"/"diode_supercells"/f"{tag}.aedt"
            logfiles=[str(p) for p in ownproject.parent.glob(f"{tag}*.log") if p.is_file()]
            queue["rows"][-1]["own_message_logs"]=logfiles
            write(STATUS,queue)
        if not success:
            queue.update(status="STOP_SOLVER_FAILURE",stop_state=state,end_utc=utc())
            write(STATUS,queue);break
        if state=="ON_B":
            # Stop after both ON cases only if neither has a Q<150 radiating mode
            # anywhere in 5.4-6.0 GHz. Q<10 modes are known PML/box modes.
            on=[]
            for name in ("ON_A","ON_B"):
                tag=(OUT/f"{name}.done").read_text(encoding="utf-8").strip()
                modes=json.loads((OUT/f"{tag}_modes.json").read_text(encoding="utf-8"))["modes"]
                on += [m for m in modes if 5.4<=m["f_ghz"]<=6.0 and 10<=m["Q"]<150]
            if not on:
                queue.update(status="STOP_BOTH_ON_FAIL",end_utc=utc())
                write(STATUS,queue);break
    else:
        queue.update(status="COMPLETE",end_utc=utc());write(STATUS,queue)
    write(FINAL,{"status":queue["status"],"completed":[p.stem for p in OUT.glob("*.done")],
                 "utc":utc()})


if __name__=="__main__":
    try: run()
    except BaseException:
        OUT.mkdir(parents=True,exist_ok=True)
        write(OUT/"queue_failure.json",{"utc":utc(),"traceback":traceback.format_exc()})
        write(FINAL,{"status":"CONTROLLER_FAILURE","utc":utc()})
        raise
