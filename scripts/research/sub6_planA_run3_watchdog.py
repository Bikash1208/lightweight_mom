"""Stop only this controller's recorded AEDT PID if its controller dies."""
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import psutil


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    parent_pid = int(sys.argv[1])
    parent_created = float(sys.argv[2])
    status_file = Path(sys.argv[3])
    log_file = Path(sys.argv[4])
    while True:
        time.sleep(10)
        try:
            status = json.loads(status_file.read_text())
            if status.get("status") in {"COMPLETE", "ABORTED"}:
                return
            parent = psutil.Process(parent_pid)
            alive = abs(parent.create_time() - parent_created) < 1 and parent.is_running()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            alive = False
        if alive:
            continue
        status = json.loads(status_file.read_text())
        pid = status.get("own_aedt_pid")
        created = status.get("own_aedt_created")
        stopped = []
        if pid and created:
            try:
                proc = psutil.Process(pid)
                if abs(proc.create_time() - created) < 1:
                    own = proc.children(recursive=True) + [proc]
                    for item in own:
                        try:
                            item.terminate()
                            stopped.append(item.pid)
                        except psutil.Error:
                            pass
            except psutil.Error:
                pass
        status["status"] = "ABORTED"
        status["error"] = "CONTROLLER_DIED_WATCHDOG_STOPPED_OWN_AEDT"
        status["watchdog_stop_pids"] = stopped
        status["ended_utc"] = utc()
        status_file.write_text(json.dumps(status, indent=2))
        with log_file.open("a", encoding="utf-8") as f:
            f.write("\nCONTROLLER_DEATH_TRACEBACK\n")
            f.write("".join(traceback.format_stack()))
            f.write(f"\nrecorded AEDT PID={pid}; stopped={stopped}\n")
        return


if __name__ == "__main__":
    try:
        main()
    except BaseException:
        with Path(sys.argv[4]).open("a", encoding="utf-8") as f:
            f.write(traceback.format_exc())
        raise
