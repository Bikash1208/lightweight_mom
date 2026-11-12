"""Independent heartbeat and RAM sampler for an attached Run 3 solve."""
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


def write(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2))
    os.replace(tmp, path)


def main():
    controller_pid = int(sys.argv[1])
    controller_created = float(sys.argv[2])
    status_path = Path(sys.argv[3])
    out = Path(sys.argv[4])
    peaks = {}
    last_heartbeat = 0.0
    while True:
        try:
            status = json.loads(status_path.read_text())
            pid = status.get("own_aedt_pid")
            created = status.get("own_aedt_created")
            run = status.get("active_run")
            if run and pid and created:
                try:
                    proc = psutil.Process(pid)
                    if abs(proc.create_time() - created) < 1:
                        tree = [proc] + proc.children(recursive=True)
                        rss = sum(p.memory_info().rss for p in tree if p.is_running())
                        peaks[run] = max(peaks.get(run, 0), rss)
                except psutil.Error:
                    pass
            now = time.monotonic()
            if now - last_heartbeat >= 60:
                write(out, {"heartbeat_utc": utc(), "controller_pid": controller_pid,
                            "active_run": run, "controller_status": status.get("status"),
                            "peak_ram_gb_by_run": {k: round(v / 2**30, 3) for k, v in peaks.items()}})
                last_heartbeat = now
            try:
                parent = psutil.Process(controller_pid)
                alive = parent.is_running() and abs(parent.create_time() - controller_created) < 1
            except psutil.Error:
                alive = False
            if not alive or status.get("status") in {"COMPLETE", "ABORTED"}:
                write(out, {"heartbeat_utc": utc(), "controller_pid": controller_pid,
                            "active_run": run, "controller_status": status.get("status"),
                            "controller_alive": alive,
                            "peak_ram_gb_by_run": {k: round(v / 2**30, 3) for k, v in peaks.items()}})
                return
            time.sleep(5)
        except BaseException:
            with out.with_suffix(".log").open("a", encoding="utf-8") as f:
                f.write(traceback.format_exc())
            raise


if __name__ == "__main__":
    main()
