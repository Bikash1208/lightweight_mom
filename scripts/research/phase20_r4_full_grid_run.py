"""Phase 20 R4 full approved grid (f 8-12 GHz step 0.1, theta 0-45 step 5, phi 0-345 step 15, four states, mesh 60, F = 2n) with a per-task checkpoint, a priority order and a
wall-clock stop. Server-authored, MODEL_ONLY. Uses `phase20_r4_build_library.solve_point` unchanged (same solver call as the R4 library build).

Order: (1) f = 9.5-11.5 GHz, all theta/phi, mesh 60; (2) the mesh-90 correction points (CORR_F x CORR_POINTS, R4's set); (3) the remaining frequencies, mesh 60.
No mesh-120 points (deferred). Every task is one `.npy` file (4 states x 2 x 2 complex), written atomically; a rerun skips finished tasks. Tasks that cannot finish before the
deadline are not started (mesh 60: 5 min margin, mesh 90: 40 min margin), so the run ends cleanly at the stop time.

    python scripts/research/phase20_r4_full_grid_run.py run <out_dir> [--workers 40] [--hours 9]
    python scripts/research/phase20_r4_full_grid_run.py status <out_dir>
    python scripts/research/phase20_r4_full_grid_run.py consolidate <out_dir>     # per-frequency chunks of the finished frequencies (+ manifest)
"""

import os
import sys

for _k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase20_r4_build_library as lib  # noqa: E402

MARGIN_S = {60: 300.0, 90: 2400.0}
LIB_STEPS = {"f_step": 0.1, "theta_step": 5.0, "phi_step": 15.0}


def task_name(f, th, ph, mesh):
    return f"m{mesh}_f{f:.2f}_t{th:g}_p{ph:g}.npy"


def build_tasks():
    """Ordered task list [(f, theta, phi, mesh)] and the group of each."""
    lib.configure_grid(**LIB_STEPS)
    pts = lib.grid_points()
    freqs = list(lib.GRID_F)
    g1 = [f for f in freqs if 9.5 - 1e-9 <= f <= 11.5 + 1e-9]
    g3 = [f for f in freqs if f not in g1]
    t1 = [(f, th, ph, 60) for f in g1 for th, ph in pts]
    t2 = [(f, th, ph, 90) for f in lib.CORR_F for th, ph in lib.CORR_POINTS]
    t3 = [(f, th, ph, 60) for f in g3 for th, ph in pts]
    return [(t, 1) for t in t1] + [(t, 2) for t in t2] + [(t, 3) for t in t3]


def work(args):
    task, out_dir, deadline = args
    f, th, ph, mesh = task
    path = Path(out_dir) / "tasks" / task_name(f, th, ph, mesh)
    if path.exists():
        return task, "cached", 0.0
    if time.monotonic() > deadline - MARGIN_S[mesh]:
        return task, "skipped_deadline", 0.0
    t0 = time.monotonic()
    try:
        out = lib.solve_point((f, th, ph, mesh))
        tmp = path.with_name("." + path.name + ".partial.npy")
        np.save(tmp, out)
        os.replace(tmp, path)
    except Exception as e:  # noqa: BLE001 - report and leave the task for the rerun
        return task, f"error:{type(e).__name__}:{e}", time.monotonic() - t0
    return task, "done", time.monotonic() - t0


def run(out_dir, workers, hours):
    import provenance as prov

    out = Path(out_dir)
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    tasks = build_tasks()
    todo = [t for t, _ in tasks if not (out / "tasks" / task_name(*t)).exists()]
    total = len(tasks)
    print(f"{total} tasks ({sum(g == 1 for _, g in tasks)} + {sum(g == 2 for _, g in tasks)} + {sum(g == 3 for _, g in tasks)} by group), {total - len(todo)} already done, "
          f"{len(todo)} to do, {workers} workers, stop after {hours} h", flush=True)
    rec = prov.provenance({"job": "r4-full-grid", "grid": LIB_STEPS, "hours_cap": hours, "tasks_total": total, "tasks_already_done": total - len(todo)}, workers=workers)
    (out / "provenance_start.json").write_text(json.dumps(rec, indent=1))
    t_start = time.monotonic()
    deadline = t_start + hours * 3600.0
    n_done = n_err = n_skip = 0
    busy = 0.0
    last = t_start
    with Pool(workers) as pool:
        for task, status, dt in pool.imap_unordered(work, [(t, str(out), deadline) for t in todo], chunksize=1):
            if status == "done":
                n_done += 1
                busy += dt
            elif status == "skipped_deadline":
                n_skip += 1
            elif status.startswith("error"):
                n_err += 1
                print(f"  ERROR {task}: {status}", flush=True)
            now = time.monotonic()
            if (n_done and n_done % 100 == 0) or now - last > 600:
                last = now
                el = now - t_start
                rate = n_done / el
                left = len(todo) - n_done - n_skip - n_err
                print(f"[{el / 3600:.2f} h] done {n_done}/{len(todo)} this run, errors {n_err}, skipped(deadline) {n_skip}, {rate * 3600:.0f} tasks/h, "
                      f"mean {busy / max(n_done, 1):.0f} s/task, ETA {left / max(rate, 1e-9) / 3600:.1f} h", flush=True)
    el = time.monotonic() - t_start
    print(f"finished: done {n_done}, errors {n_err}, skipped at the deadline {n_skip}, wall {el / 3600:.2f} h", flush=True)
    (out / "provenance_end.json").write_text(json.dumps(prov.finish(rec) | {"tasks_done_this_run": n_done, "errors": n_err, "skipped_deadline": n_skip, "wall_s": el}, indent=1))
    status(out_dir)


def status(out_dir):
    out = Path(out_dir)
    tasks = build_tasks()
    have = {g: [0, 0] for g in (1, 2, 3)}
    for t, g in tasks:
        have[g][1] += 1
        have[g][0] += (out / "tasks" / task_name(*t)).exists()
    for g, name in ((1, "f 9.5-11.5 GHz, mesh 60"), (2, "mesh-90 correction points"), (3, "remaining f, mesh 60")):
        print(f"group {g} ({name}): {have[g][0]}/{have[g][1]}")
    print("mesh 120: DEFERRED (not run)")
    return have


def consolidate(out_dir):
    """Per-frequency chunks (same layout as `phase20_r4_build_library.chunk_path` files: f_ghz, points, s, mesh) for every fully finished frequency, plus the manifest."""
    out = Path(out_dir)
    lib.configure_grid(**LIB_STEPS)
    pts = lib.grid_points()
    (out / "chunks").mkdir(exist_ok=True)
    rows = []
    n_f = 0
    for f in lib.GRID_F:
        files = [out / "tasks" / task_name(f, th, ph, 60) for th, ph in pts]
        if not all(p.exists() for p in files):
            continue
        s = np.array([np.load(p) for p in files])
        lib.CHUNKS = out / "chunks"
        lib.save_atomic(lib.chunk_path("grid", f), f_ghz=f, points=np.array(pts), s=s, mesh=60)
        n_f += 1
    f90 = []
    for f in lib.CORR_F:
        files = [out / "tasks" / task_name(f, th, ph, 90) for th, ph in lib.CORR_POINTS]
        if all(p.exists() for p in files):
            lib.save_atomic(lib.chunk_path("corr90", f), f_ghz=f, points=np.array(lib.CORR_POINTS), s=np.array([np.load(p) for p in files]), mesh=90)
            f90.append(f)
    h = hashlib.sha256()
    n_tasks = 0
    for p in sorted((out / "tasks").glob("m*.npy")):
        h.update(p.name.encode() + hashlib.sha256(p.read_bytes()).digest())
        n_tasks += 1
    for p in sorted((out / "chunks").glob("*.npz")) + [q for q in (out / "provenance_start.json", out / "provenance_end.json") if q.exists()]:
        b = p.read_bytes()
        rows.append({"file": str(p.relative_to(out)), "size": len(b), "sha256": hashlib.sha256(b).hexdigest()})
    commit = os.popen(f'git -C "{ROOT}" rev-parse HEAD').read().strip()
    manifest = {
        "job": "r4-full-grid",
        "label": "SERVER_AUTHORED_PENDING_MAC_REVIEW",
        "data_dir": str(out),
        "producing_commit": commit,
        "machine": os.environ.get("COMPUTERNAME", "unknown"),
        "grid": LIB_STEPS,
        "frequencies_complete_mesh60": n_f,
        "mesh90_correction_frequencies_complete": f90,
        "mesh120": "DEFERRED",
        "phi_range": "0-345 deg, full (the mirror check failed: MIRROR_NOT_USABLE)",
        "task_files": {"count": n_tasks, "combined_sha256": h.hexdigest()},
        "files": rows,
    }
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    print(f"consolidated {n_f} complete frequencies (mesh 60) and {len(f90)} mesh-90 correction chunks; {n_tasks} task files; MANIFEST.json written")


if __name__ == "__main__":
    cmd, outd = sys.argv[1], sys.argv[2]
    opts = sys.argv[3:]
    w = int(opts[opts.index("--workers") + 1]) if "--workers" in opts else 40
    hrs = float(opts[opts.index("--hours") + 1]) if "--hours" in opts else 9.0
    {"run": lambda: run(outd, w, hrs), "status": lambda: status(outd), "consolidate": lambda: consolidate(outd)}[cmd]()
