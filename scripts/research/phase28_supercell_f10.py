"""Phase 28: the Phase 26 SOLVER_SUPERCELL table (beta and alpha of a slot-angle-rotating row vs modulation period) rebuilt at another frequency (10.0 GHz, the PI's
operating frequency, 2026-10-26). SAME METHOD, SAME CODE: `phase26_supercell_table.run_point` (periodic MoM, analytic-continuation complex root) is called unchanged; only the
module constants `F` and `K0` of the two Phase 26 scripts are set at run time, so no Phase 26 script or solver module is modified. `validate` reproduces a committed 10.54 GHz
table point with the overridden machinery before anything is trusted.

    python scripts/research/phase28_supercell_f10.py validate
    python scripts/research/phase28_supercell_f10.py task   <outdir> <L> <W> <N>
    python scripts/research/phase28_supercell_f10.py run    <outdir> [--workers 32] [--mem-gb 190] [--deadline ISO8601] [--tier 1|2|all]
    python scripts/research/phase28_supercell_f10.py assemble <outdir>

`run` is the orchestrator: one subprocess per (slot, N) task (1 BLAS thread each), largest first, a task starts only if the free physical memory covers its estimated peak and
(optionally) it is estimated to finish before the deadline; every finished task is a checkpoint file `task_<L>x<W>_N<N>.json` (a rerun skips finished tasks, so the job is resumable);
`manifest.json` records every task (status, seconds, peak working set, sha256 of its file) and the provenance.
Labels: SOLVER_SUPERCELL, MODEL_ONLY, NOT_HFSS_VALIDATED, GAMMA_0_COLLINEAR_ONLY.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import ctypes  # noqa: E402
import datetime  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

FREQ_HZ = 10.0e9
SLOTS = [(4.4, 1.2), (4.6, 1.6), (4.8, 1.6), (4.8, 2.0)]  # the four slot sizes of the committed 10.54 GHz table
DIAM, MESH_N = 6.0, 24
# modulation periods (cells): K/k0 = lambda0 / (N P). Tier 1 reproduces the committed table's five periods (beams up to ~30 deg: K in 0.56..1.56 k0);
# tier 2 extends the range (N = 2 for K above, N = 10, 12, 16 for K below: beams 40..50 deg need K down to ~0.3 k0)
TIER1 = [3, 4, 5, 6, 8]
TIER2 = [2, 10, 12, 16]
# cost model, calibrated from the first timings (single thread): seconds ~ 145 (N/3)**2.8, peak memory ~ 0.1 N**2 GB (checked against the manifest)
COST_S_AT_N3, COST_EXP = 145.0, 2.8  # measured: N = 3 145 s, N = 5 569 s (exponent 2.67; 2.8 used to be conservative)
MEM_GB_PER_N2 = 0.1  # a guess; the peak working set of every task is recorded in the manifest


def set_frequency(freq_hz):
    """Point the two Phase 26 scripts' module constants at `freq_hz` (they are read at call time inside the functions)."""
    import phase26_ramp_locality as rl
    import phase26_supercell_table as st

    from src.frequency import free_space_wavenumber

    k0 = free_space_wavenumber(freq_hz)
    rl.F, rl.K0 = freq_hz, k0
    st.F, st.K0, st.LAMBDA0 = freq_hz, k0, 299792458.0 / freq_hz
    return st


def peak_working_set_gb():
    class PMC(ctypes.Structure):
        _fields_ = [
            ("cb", ctypes.c_ulong),
            ("PageFaultCount", ctypes.c_ulong),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    from ctypes import wintypes

    k32, psapi = ctypes.WinDLL("kernel32"), ctypes.WinDLL("psapi")
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    pmc = PMC()
    pmc.cb = ctypes.sizeof(PMC)
    if not psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
        return float("nan")
    return pmc.PeakWorkingSetSize / 2**30


def free_memory_gb():
    class MS(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    m = MS()
    m.dwLength = ctypes.sizeof(MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 2**30


def task_name(slot, n):
    return f"task_{slot[0]}x{slot[1]}_N{n}"


def run_task(outdir, slot, n, freq_hz=FREQ_HZ):
    st = set_frequency(freq_hz)
    t0 = time.time()
    rec = st.run_point(n, slot, DIAM, MESH_N)
    rec["frequency_ghz"] = freq_hz / 1e9
    rec["task_seconds"] = time.time() - t0
    rec["peak_working_set_gb"] = peak_working_set_gb()
    out = Path(outdir) / (task_name(slot, n) + ".json")
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec, indent=1), encoding="utf-8", newline="\n")
    tmp.replace(out)  # atomic: a checkpoint is either complete or absent
    return rec


def validate():
    """Reproduce a committed 10.54 GHz table point with the overridden machinery (slot 4.4 x 1.2, N = 3, mesh 24)."""
    ref = json.loads((ROOT / "results/phase26/supercell_table_4p4x1p2.json").read_text(encoding="utf-8"))
    want = next(p for p in ref["points"] if p["n_cells"] == 3)
    st = set_frequency(10.54e9)
    got = st.run_point(3, (4.4, 1.2), DIAM, MESH_N)
    rows = []
    for k in ("K_over_k0", "beta_real_root_over_k0", "beta_over_k0", "alpha_np_per_m", "alpha_leak_np_per_m"):
        rows.append(
            {"quantity": k, "this_run": got[k], "committed": want[k], "rel_diff": abs(got[k] - want[k]) / abs(want[k])}
        )
    out = {
        "label": "VALIDATION_10P54_REPRODUCTION",
        "rows": rows,
        "max_rel_diff": max(r["rel_diff"] for r in rows),
        "ok": max(r["rel_diff"] for r in rows) < 1e-6,
    }
    print(json.dumps(out, indent=1))
    return out


def est_seconds(n):
    return COST_S_AT_N3 * (n / 3.0) ** COST_EXP


def est_mem_gb(n):
    return MEM_GB_PER_N2 * n**2


def tasks_for(tier):
    ns = {"1": TIER1, "2": TIER2, "all": TIER1 + TIER2}[tier]
    return [(s, n) for n in ns for s in SLOTS]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def running_task_processes(outdir):
    """`{(slot, n): pid}` of task subprocesses of this script already running (so a restarted orchestrator adopts them instead of starting duplicates)."""
    cmd = 'Get-CimInstance Win32_Process -Filter "Name=\'python.exe\'" | ForEach-Object { "$($_.ProcessId)|$($_.CommandLine)" }'
    res = subprocess.run(["powershell", "-NoProfile", "-Command", cmd], capture_output=True, text=True, timeout=60)
    out = {}
    for line in res.stdout.splitlines():
        if "phase28_supercell_f10.py" not in line or " task " not in line:
            continue
        pid, _, cl = line.partition("|")
        parts = cl.replace('"', "").split(" task ", 1)[1].split()
        if len(parts) >= 4 and Path(parts[0]).resolve() == Path(outdir).resolve():
            out[((float(parts[1]), float(parts[2])), int(parts[3]))] = int(pid)
    return out


def orchestrate(outdir, workers, mem_gb, deadline, tier):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pending = sorted(tasks_for(tier), key=lambda t: -t[1])  # largest first
    pending = [t for t in pending if not (outdir / (task_name(*t) + ".json")).exists()]
    running = {}  # name -> (Popen, start, (slot, n), logfile)
    adopted = running_task_processes(outdir)
    pending = [t for t in pending if t not in adopted]
    if adopted:
        print("adopting running tasks:", sorted((s, n) for s, n in adopted), flush=True)
    manifest_path = outdir / "manifest.json"
    skipped = []
    print(f"{len(pending)} tasks pending, {workers} workers, {mem_gb} GB budget, deadline {deadline}", flush=True)
    while pending or running or adopted:
        for key in list(adopted):
            if not (Path(outdir) / (task_name(*key) + ".json")).exists() and _alive(adopted[key]):
                continue
            del adopted[key]
        used_adopted = sum(est_mem_gb(n) for (_, n) in adopted)
        # reap
        for name in list(running):
            p, t0, (slot, n), lf = running[name]
            if p.poll() is not None:
                lf.close()
                print(f"  finished {name} rc={p.returncode} in {time.time() - t0:.0f} s", flush=True)
                del running[name]
        # start
        used = used_adopted + sum(est_mem_gb(n) for (_, _, (_, n), _) in running.values())
        i = 0
        while i < len(pending) and len(running) < workers:
            slot, n = pending[i]
            need = est_mem_gb(n)
            finish = time.time() + est_seconds(n) * 1.5
            if deadline is not None and finish > deadline:
                skipped.append((slot, n))
                pending.pop(i)
                continue
            if used + need > mem_gb or free_memory_gb() < need + 8:
                i += 1  # try a smaller one that fits
                continue
            name = task_name(slot, n)
            lf = open(outdir / (name + ".log"), "w", encoding="utf-8")  # noqa: SIM115 - handle lives as long as the subprocess
            p = subprocess.Popen(
                [sys.executable, str(Path(__file__)), "task", str(outdir), str(slot[0]), str(slot[1]), str(n)],
                stdout=lf,
                stderr=subprocess.STDOUT,
            )
            running[name] = (p, time.time(), (slot, n), lf)
            used += need
            pending.pop(i)
            print(
                f"  started {name} (est {est_seconds(n):.0f} s, {need:.1f} GB); {len(running)} running, {len(pending)} pending",
                flush=True,
            )
        write_manifest(outdir, manifest_path, tier, skipped)
        time.sleep(10)
    write_manifest(outdir, manifest_path, tier, skipped)
    print("done; skipped for the deadline:", skipped, flush=True)


def _alive(pid):
    res = subprocess.run(
        ["powershell", "-NoProfile", "-Command", f"(Get-Process -Id {pid} -ErrorAction SilentlyContinue) -ne $null"],
        capture_output=True,
        text=True,
        timeout=60,
    )
    return res.stdout.strip() == "True"


def write_manifest(outdir, path, tier, skipped):
    import provenance as prov

    rows = []
    for slot, n in tasks_for("all"):
        f = Path(outdir) / (task_name(slot, n) + ".json")
        row = {
            "slot_L_W_mm": list(slot),
            "n_cells": n,
            "file": f.name,
            "status": "DONE" if f.exists() else ("SKIPPED_FOR_DEADLINE" if (slot, n) in skipped else "PENDING"),
        }
        if f.exists():
            r = json.loads(f.read_text(encoding="utf-8"))
            row.update(
                {
                    "sha256": sha256(f),
                    "seconds": r["task_seconds"],
                    "peak_working_set_gb": r["peak_working_set_gb"],
                    "alpha_leak_np_per_m": r["alpha_leak_np_per_m"],
                    "K_over_k0": r["K_over_k0"],
                }
            )
        rows.append(row)
    rec = prov.provenance(
        {
            "tool": "phase28_supercell_f10.py",
            "frequency_ghz": FREQ_HZ / 1e9,
            "slots": SLOTS,
            "tier1": TIER1,
            "tier2": TIER2,
            "mesh_per_cell": MESH_N,
            "diameter_mm": DIAM,
        }
    )
    doc = {
        "label": "PHASE28_SUPERCELL_TABLE_10GHZ_MANIFEST",
        "generated": datetime.datetime.now().astimezone().isoformat(),
        "provenance": rec,
        "tasks": rows,
    }
    tmp = Path(path).with_suffix(".tmp")
    tmp.write_text(json.dumps(doc, indent=1), encoding="utf-8", newline="\n")
    tmp.replace(path)


def assemble(outdir):
    """Combine the finished task files into one table file per slot, in the committed Phase 26 format (`SupercellTable.from_files` reads them)."""
    outdir = Path(outdir)
    n_done = 0
    for slot in SLOTS:
        pts = []
        for n in sorted(TIER1 + TIER2, reverse=True):
            f = outdir / (task_name(slot, n) + ".json")
            if f.exists():
                pts.append(json.loads(f.read_text(encoding="utf-8")))
        if not pts:
            continue
        pts.sort(key=lambda p: -p["n_cells"])
        doc = {
            "labels": ["SOLVER_SUPERCELL", "MODEL_ONLY", "NOT_HFSS_VALIDATED", "GAMMA_0_COLLINEAR_ONLY"],
            "mesh_per_cell": MESH_N,
            "frequency_ghz": FREQ_HZ / 1e9,
            "diameter_mm": DIAM,
            "pitch_mm": 6.4,
            "points": [{k: v for k, v in p.items() if k not in ("task_seconds", "peak_working_set_gb")} for p in pts],
        }
        name = f"supercell_table_{str(slot[0]).replace('.', 'p')}x{str(slot[1]).replace('.', 'p')}.json"
        (outdir / name).write_text(json.dumps(doc, indent=1), encoding="utf-8", newline="\n")
        n_done += len(pts)
    print("assembled", n_done, "points")


if __name__ == "__main__":
    cmd, a = sys.argv[1], sys.argv[2:]
    if cmd == "validate":
        validate()
    elif cmd == "task":
        run_task(a[0], (float(a[1]), float(a[2])), int(a[3]))
    elif cmd == "run":
        w = int(a[a.index("--workers") + 1]) if "--workers" in a else 32
        m = float(a[a.index("--mem-gb") + 1]) if "--mem-gb" in a else 190.0
        d = datetime.datetime.fromisoformat(a[a.index("--deadline") + 1]).timestamp() if "--deadline" in a else None
        t = a[a.index("--tier") + 1] if "--tier" in a else "all"
        orchestrate(a[0], w, m, d, t)
    elif cmd == "assemble":
        assemble(a[0])
    else:
        raise SystemExit(__doc__)
