r"""Memory-guarded launcher for the jobs/r5-supercell-gradients tasks (server 1). One process per task via the job's own script (`--only NAME`, unchanged); this file only decides WHEN to start.

Start rule (user, 2026-10-24): start a task only if free RAM >= 30 GB + its expected peak. Matrix order N = 1510 * nx * ny (measured on 1x1, 2x1, 3x1, 1x3: exact), expected peak =
1.3 * 16 * N^2 bytes + 1 GB (the factor 1.3 comes from the measured random-job worker: N = 24160 -> 16 N^2 = 9.34 GB, measured peak 11.92 GB = 1.28x). Threads: 3 BLAS threads per task;
running tasks * 3 + 9 (random job workers) <= 48. No new task starts after START_CUTOFF (21:42 local); tasks already running may finish until HARD_STOP (23:59 local), then they are
terminated (their task JSON is not written -> they rerun from scratch). Each task process is logged to C:\lw_jobs\gradients_logs\<name>.log. Resumable: tasks whose JSON exists are skipped.
"""

import ctypes
import ctypes.wintypes as wt
import datetime as dt
import json
import os
import subprocess
import sys
import time
from pathlib import Path

TREE = Path(r"C:\lw_jobs\gradients")
OUT = TREE / "results" / "phase20_r5" / "supercell_gradients"
PY = r"C:\lightweight_mom\.venv\Scripts\python.exe"
LOGDIR = Path(r"C:\lw_jobs\gradients_logs")
EVENTS = Path(r"C:\lw_jobs\gradients_orchestrator.log")
THREADS = 3
RANDOM_THREADS = 9
THREAD_CAP = 48
MIN_FREE_GB = 30.0
N_PER_CELL_CELL = 1510  # unknowns per cell at 48 cells per period
FACTOR, MARGIN_GB = 1.3, 1.0


def today_at(h, m):
    n = dt.datetime.now()
    return n.replace(hour=h, minute=m, second=0, microsecond=0)


START_CUTOFF = today_at(21, 42)
HARD_STOP = today_at(23, 59)

k32 = ctypes.WinDLL("kernel32", use_last_error=True)


class MEMSTAT(ctypes.Structure):
    _fields_ = [("dwLength", wt.DWORD), ("dwMemoryLoad", wt.DWORD), ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong), ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong), ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong), ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]


def free_gb():
    m = MEMSTAT()
    m.dwLength = ctypes.sizeof(m)
    k32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 1e9


class PMC(ctypes.Structure):
    _fields_ = [("cb", wt.DWORD), ("PageFaultCount", wt.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t), ("a", ctypes.c_size_t), ("b", ctypes.c_size_t),
                ("c", ctypes.c_size_t), ("d", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t)]


class PE32(ctypes.Structure):
    _fields_ = [("dwSize", wt.DWORD), ("cntUsage", wt.DWORD), ("th32ProcessID", wt.DWORD), ("th32DefaultHeapID", ctypes.c_size_t),
                ("th32ModuleID", wt.DWORD), ("cntThreads", wt.DWORD), ("th32ParentProcessID", wt.DWORD), ("pcPriClassBase", wt.LONG),
                ("dwFlags", wt.DWORD), ("szExeFile", ctypes.c_char * 260)]


def descendants(pid):
    snap = k32.CreateToolhelp32Snapshot(2, 0)
    pe = PE32()
    pe.dwSize = ctypes.sizeof(pe)
    kids = {}
    ok = k32.Process32First(snap, ctypes.byref(pe))
    while ok:
        kids.setdefault(pe.th32ParentProcessID, []).append(pe.th32ProcessID)
        ok = k32.Process32Next(snap, ctypes.byref(pe))
    k32.CloseHandle(snap)
    out, stack = [], [pid]
    while stack:
        p = stack.pop()
        for c in kids.get(p, []):
            out.append(c)
            stack.append(c)
    return out


def peak_ws_gb(pid):
    k32.OpenProcess.restype = wt.HANDLE
    h = k32.OpenProcess(0x1000 | 0x0010, False, pid)  # QUERY_LIMITED_INFORMATION | VM_READ
    if not h:
        return 0.0, 0.0
    p = PMC()
    p.cb = ctypes.sizeof(p)
    k32.K32GetProcessMemoryInfo.argtypes = [wt.HANDLE, ctypes.POINTER(PMC), wt.DWORD]
    k32.K32GetProcessMemoryInfo(h, ctypes.byref(p), p.cb)
    k32.CloseHandle(h)
    return p.PeakWorkingSetSize / 1e9, p.WorkingSetSize / 1e9


def log(msg):
    line = f"{dt.datetime.now().isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with open(EVENTS, "a") as fh:
        fh.write(line + "\n")


def task_list():
    sys.path.insert(0, str(TREE))
    sys.path.insert(0, str(TREE / "scripts" / "research"))
    import phase20_r5_supercell as s

    from src.ris.calibrated_library import CalibratedLibrary

    out = []
    for t in s.gradient_tasks(CalibratedLibrary.load()):
        n = N_PER_CELL_CELL * t["nx"] * t["ny"]
        out.append({"name": t["name"], "nx": t["nx"], "ny": t["ny"], "N": n, "est_gb": FACTOR * 16 * n * n / 1e9 + MARGIN_GB})
    return sorted(out, key=lambda t: -t["N"])  # largest first (the slowest tasks start first)


def main():
    LOGDIR.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = task_list()
    log(f"orchestrator start; start cutoff {START_CUTOFF:%H:%M}, hard stop {HARD_STOP:%H:%M}; tasks: " + ", ".join(f"{t['name']}(N={t['N']}, est {t['est_gb']:.1f} GB)" for t in tasks))
    pending = [t for t in tasks if not (OUT / f"{t['name']}.json").exists()]
    running = {}
    env = dict(os.environ, OPENBLAS_NUM_THREADS=str(THREADS), OMP_NUM_THREADS=str(THREADS), MKL_NUM_THREADS=str(THREADS))
    last_start = 0.0
    while pending or running:
        now = dt.datetime.now()
        # finished tasks
        for name, r in list(running.items()):
            rc = r["proc"].poll()
            pk = max([peak_ws_gb(p)[0] for p in descendants(r["proc"].pid)] + [peak_ws_gb(r["proc"].pid)[0], r["peak"]])
            r["peak"] = pk
            if rc is not None:
                ok = (OUT / f"{name}.json").exists()
                log(f"FINISHED {name}: exit {rc}, json {'present' if ok else 'MISSING'}, wall {(time.time() - r['t0']) / 3600:.2f} h, peak working set {pk:.2f} GB (N={r['N']}, estimated {r['est']:.1f} GB)")
                del running[name]
        if now >= HARD_STOP:
            for name, r in running.items():
                log(f"HARD STOP: terminating {name} after {(time.time() - r['t0']) / 3600:.2f} h (peak {r['peak']:.2f} GB); it reruns from scratch")
                subprocess.run(["taskkill", "/F", "/T", "/PID", str(r["proc"].pid)], capture_output=True)
            break
        if now >= START_CUTOFF and pending:
            log(f"start cutoff reached: {len(pending)} tasks not started: " + ", ".join(t["name"] for t in pending))
            pending = []
        # start the next task if the guards allow it
        if pending and time.time() - last_start > 20:
            t = pending[0]
            fg = free_gb()
            threads_after = (len(running) + 1) * THREADS + RANDOM_THREADS
            if fg >= MIN_FREE_GB + t["est_gb"] and threads_after <= THREAD_CAP:
                lf = open(LOGDIR / f"{t['name']}.log", "a")
                p = subprocess.Popen([PY, "-u", "scripts/research/phase20_r5_supercell.py", "gradients", str(OUT.relative_to(TREE)), "--workers", "1", "--only", t["name"]],
                                     cwd=str(TREE), env=env, stdout=lf, stderr=subprocess.STDOUT)
                running[t["name"]] = {"proc": p, "t0": time.time(), "N": t["N"], "est": t["est_gb"], "peak": 0.0}
                last_start = time.time()
                pending.pop(0)
                log(f"STARTED {t['name']}: N={t['N']}, expected peak {t['est_gb']:.1f} GB, free RAM {fg:.1f} GB (needs {MIN_FREE_GB + t['est_gb']:.1f}), threads after start {threads_after}/{THREAD_CAP}")
            else:
                log(f"WAIT {t['name']}: free RAM {fg:.1f} GB < {MIN_FREE_GB + t['est_gb']:.1f} needed or threads {threads_after} > {THREAD_CAP}")
                last_start = time.time() + 40
        time.sleep(15)
    log("orchestrator end")


if __name__ == "__main__":
    main()
