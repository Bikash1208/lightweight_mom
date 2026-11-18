"""Stop one explicitly recorded local AEDT run by exact owned PIDs."""

import json
from pathlib import Path
import sys
from datetime import datetime, timezone

import psutil


status_path = Path(sys.argv[1]).resolve()
root = Path(__file__).resolve().parents[2]
if root not in status_path.parents:
    raise RuntimeError("STATUS_OUTSIDE_OWN_WORKTREE")
status = json.loads(status_path.read_text(encoding="utf-8"))
if status["status"] not in ("SOLVING", "BUILDING"):
    raise RuntimeError(f"NOT_ACTIVE {status['status']}")
aedt_pid = int(status["aedt_pid"])
controller_pid = int(status["controller_pid"])
watchdog_pid = int(status.get("external_watchdog_pid", 0))
expected = str(Path(r"C:\Program Files\ANSYS Inc\v252\AnsysEM\ansysedt.exe")).lower()
record = {"reason": "USER_SERVER_RESTART", "utc": datetime.now(timezone.utc).isoformat(),
          "status_path": str(status_path), "own_aedt_pid": aedt_pid,
          "own_controller_pid": controller_pid, "own_watchdog_pid": watchdog_pid,
          "terminated_exact_pids": []}
if psutil.pid_exists(aedt_pid):
    aedt = psutil.Process(aedt_pid)
    if aedt.exe().lower() != expected:
        raise RuntimeError("AEDT_PID_EXE_MISMATCH")
    targets = aedt.children(recursive=True) + [aedt]
    record["terminated_exact_pids"].extend(p.pid for p in targets)
    for proc in targets:
        try:
            if proc.is_running():
                proc.terminate()
        except psutil.NoSuchProcess:
            pass
if psutil.pid_exists(controller_pid):
    ctrl = psutil.Process(controller_pid)
    if Path(ctrl.exe()).name.lower() not in ("python.exe", "pythonw.exe"):
        raise RuntimeError("CONTROLLER_PID_EXE_MISMATCH")
    record["terminated_exact_pids"].append(controller_pid)
    ctrl.terminate()
if watchdog_pid and psutil.pid_exists(watchdog_pid):
    watch = psutil.Process(watchdog_pid)
    if Path(watch.exe()).name.lower() not in ("python.exe", "pythonw.exe"):
        raise RuntimeError("WATCHDOG_PID_EXE_MISMATCH")
    record["terminated_exact_pids"].append(watchdog_pid)
    watch.terminate()
out = status_path.with_name(status_path.stem + "_user_stop.json")
out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
print(json.dumps(record))
