"""Stop only Stage 1 AEDT PIDs recorded by our failed worker controllers."""
import json
from datetime import datetime, timezone
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/sub6_planA/stage1"
RUNS = [f"ANCHOR_{v}_{s}" for v in ("V", "T8", "T5", "T2", "T1")
        for s in ("OFF", "ON")]
reports = []
for run in RUNS:
    path = BASE / run / f"{run}.json"
    d = json.loads(path.read_text())
    if d.get("run_id") != run or d.get("status") != "SOLVING":
        raise RuntimeError(f"UNEXPECTED_STATUS {run}: {d.get('status')}")
    if psutil.pid_exists(int(d["own_python_pid"])):
        raise RuntimeError(f"CONTROLLER_STILL_ALIVE {run}")
    pid = int(d["own_aedt_pid"])
    if not psutil.pid_exists(pid):
        raise RuntimeError(f"AEDT_ALREADY_GONE {run}")
    proc = psutil.Process(pid)
    if proc.name().lower() != "ansysedt.exe":
        raise RuntimeError(f"PID_NOT_AEDT {run}")
    reports.append((run, path, d, proc))

for run, path, d, proc in reports:
    children = proc.children(recursive=True)
    for q in reversed(children):
        q.terminate()
    _, alive = psutil.wait_procs(children, timeout=20)
    for q in alive:
        q.kill()
    proc.terminate()
    try:
        proc.wait(timeout=20)
    except psutil.TimeoutExpired:
        proc.kill()
    d["status"] = "ABORTED_CONTROLLER_LOST"
    d["abort_reason"] = "Own Python worker exited code 1 without final report; AEDT continued orphaned"
    d["aborted_utc"] = datetime.now(timezone.utc).isoformat()
    d["own_aedt_descendant_pids_stopped"] = [q.pid for q in children]
    path.write_text(json.dumps(d, indent=2))
    print(run, "stopped", flush=True)
