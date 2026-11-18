"""Independent heartbeat and exact-PID timeout for one locally launched AEDT solve."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time
import traceback

import psutil


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def save(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def main(status):
    heartbeat = status.with_name(status.stem + "_external_heartbeat.json")
    initial = read(status)
    aedt_pid = int(initial["aedt_pid"])
    controller_pid = int(initial["controller_pid"])
    aedt_created = psutil.Process(aedt_pid).create_time()
    controller_created = psutil.Process(controller_pid).create_time()
    started = datetime.fromisoformat(initial["solve_start_utc"])
    while True:
        state = read(status)
        aedt = psutil.Process(aedt_pid) if psutil.pid_exists(aedt_pid) else None
        ctrl = psutil.Process(controller_pid) if psutil.pid_exists(controller_pid) else None
        own_aedt_alive = bool(aedt and aedt.create_time() == aedt_created and aedt.is_running())
        own_controller_alive = bool(ctrl and ctrl.create_time() == controller_created and ctrl.is_running())
        elapsed = (datetime.now(timezone.utc) - started).total_seconds()
        record = {"utc": datetime.now(timezone.utc).isoformat(),
                  "source_status": str(status), "status": state["status"],
                  "own_aedt_pid": aedt_pid, "own_aedt_alive": own_aedt_alive,
                  "own_controller_pid": controller_pid,
                  "own_controller_alive": own_controller_alive,
                  "solve_elapsed_s": elapsed}
        save(heartbeat, record)
        if state["status"] in ("COMPLETE", "FAILED", "TIMEOUT"):
            return
        if state["status"] == "SOLVING" and (elapsed >= 1500 or not own_controller_alive):
            reason = "TIMEOUT_1500S" if elapsed >= 1500 else "CONTROLLER_DEATH"
            record["reason"] = reason
            if own_aedt_alive:
                targets = aedt.children(recursive=True) + [aedt]
                record["terminated_own_pids"] = [p.pid for p in targets]
                save(heartbeat, record)
                for proc in targets:
                    if proc.is_running():
                        proc.terminate()
            save(heartbeat, record)
            return
        time.sleep(60)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("status", type=Path)
    args = ap.parse_args()
    try:
        main(args.status)
    except BaseException:
        args.status.with_name(args.status.stem + "_watchdog_error.log").write_text(
            traceback.format_exc(), encoding="utf-8")
        raise
