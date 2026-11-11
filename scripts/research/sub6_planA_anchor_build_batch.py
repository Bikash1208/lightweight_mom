"""Build the remaining Stage 1 variants in private AEDT instances."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "stage1_work"
BASE = ROOT / "results/sub6_planA/stage1/preflight/anchor_import_preflight.aedt"
REPORT = ROOT / "results/sub6_planA/stage1/build_batch.json"
SCRIPT = ROOT / "scripts/research/sub6_planA_anchor_build.py"
PY = Path(sys.executable)


def save(d):
    REPORT.write_text(json.dumps(d, indent=2))


def main():
    jobs = [(v, s) for v in ("V", "T8", "T5", "T2", "T1")
            for s in ("OFF", "ON") if not (v == "T5" and s == "OFF")]
    free_cores = round(os.cpu_count() * (1 - psutil.cpu_percent(interval=2) / 100))
    free_gb = psutil.virtual_memory().available / 2**30
    if free_cores < len(jobs) or free_gb < 9 * len(jobs):
        raise RuntimeError(f"CPU_OR_RAM_NOT_FREE cores={free_cores} GB={free_gb:.1f}")
    d = {"status": "RUNNING", "free_cores_at_start": free_cores,
         "free_ram_gb_at_start": round(free_gb, 1), "jobs": {}}
    save(d)
    procs = {}
    handles = {}
    for v, s in jobs:
        key = f"ANCHOR_{v}_{s}"
        f = open(WORK / f"{key}_build_console.log", "w", encoding="utf-8")
        handles[key] = f
        cmd = [str(PY), str(SCRIPT), "--base", str(BASE), "--out", str(WORK),
               "--variant", v, "--state", s]
        p = subprocess.Popen(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT,
                             creationflags=subprocess.CREATE_NO_WINDOW)
        procs[key] = p
        d["jobs"][key] = {"python_pid": p.pid, "status": "BUILDING"}
        save(d)
    while procs:
        for key, p in list(procs.items()):
            code = p.poll()
            if code is None:
                continue
            handles[key].close()
            rpath = WORK / f"{key}_build.json"
            r = json.loads(rpath.read_text()) if rpath.exists() else {}
            d["jobs"][key].update({"exit_code": code, "status": r.get("status", "NO_REPORT"),
                                    "aedt_pid": r.get("own_aedt_pid"), "error": r.get("error")})
            del procs[key]
            save(d)
            if code != 0 or r.get("status") != "BUILT_UNSOLVED":
                d["status"] = "STOPPED_BUILD_FAILURE"
                save(d)
                for other, q in procs.items():
                    q.terminate()  # exact own Popen PID only
                for other, q in procs.items():
                    q.wait(timeout=30)
                    handles[other].close()
                    own = WORK / f"{other}_build.json"
                    if own.exists():
                        rr = json.loads(own.read_text())
                        apid = rr.get("own_aedt_pid")
                        if apid and psutil.pid_exists(apid):
                            psutil.Process(apid).terminate()  # exact AEDT PID recorded by own script
                return
        if procs:
            time.sleep(10)
    d["status"] = "ALL_BUILT_UNSOLVED"
    save(d)


if __name__ == "__main__":
    main()
