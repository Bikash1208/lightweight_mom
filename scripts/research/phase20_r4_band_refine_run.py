"""Phase 20 R4 band refinement (SERVER JOB `r4-band-refine`): the library's resonance band 10.4-10.9 GHz gets frequency nodes every 0.05 GHz. Adds the frequencies that are NOT on the 0.1 GHz grid
(10.45, 10.55, 10.65, 10.75, 10.85) for the full theta (0-45 step 5) / phi (0-345 step 15) grid at mesh 60, F = 2n, plus the correction-map points at those frequencies:
mesh 90 on CORR_POINTS (always), mesh 120 on CORR_POINTS (only with --mesh120; needed for the three-mesh Richardson correction of the library, see the job file). The mesh-60 correction points
are grid points, so no separate corr60 tasks exist. Per-task checkpoint (finished tasks are skipped on a rerun), wall-clock stop, same solver call as the library build. MODEL_ONLY.

    python scripts/research/phase20_r4_band_refine_run.py run <out_dir> [--workers 40] [--hours 6] [--mesh120]
    python scripts/research/phase20_r4_band_refine_run.py status <out_dir> [--mesh120]
    python scripts/research/phase20_r4_band_refine_run.py consolidate <out_dir>     # chunks in the library-build layout + MANIFEST.json
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
import phase20_r4_full_grid_run as fg  # noqa: E402

BAND_F = [10.45, 10.55, 10.65, 10.75, 10.85]
fg.MARGIN_S[120] = 5400.0  # a mesh-120 task takes ~1 h: do not start one that cannot finish before the stop time


def build_tasks(mesh120=False):
    """Ordered [(f, theta, phi, mesh)]: longest first (mesh 120, then 90, then 60) so the pool stays busy to the end."""
    lib.configure_grid(**fg.LIB_STEPS, f_list=BAND_F, corr_f_list=BAND_F)
    pts = lib.grid_points()
    t120 = [(f, th, ph, 120) for f in BAND_F for th, ph in lib.CORR_POINTS] if mesh120 else []
    t90 = [(f, th, ph, 90) for f in BAND_F for th, ph in lib.CORR_POINTS]
    t60 = [(f, th, ph, 60) for f in BAND_F for th, ph in pts]
    return t120 + t90 + t60


def run(out_dir, workers, hours, mesh120):
    import provenance as prov

    out = Path(out_dir)
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    tasks = build_tasks(mesh120)
    todo = [t for t in tasks if not (out / "tasks" / fg.task_name(*t)).exists()]
    print(
        f"{len(tasks)} tasks ({sum(t[3] == 120 for t in tasks)} mesh 120, {sum(t[3] == 90 for t in tasks)} mesh 90, {sum(t[3] == 60 for t in tasks)} mesh 60), "
        f"{len(tasks) - len(todo)} done, {len(todo)} to do, {workers} workers, stop after {hours} h",
        flush=True,
    )
    rec = prov.provenance(
        {"job": "r4-band-refine", "band_f": BAND_F, "mesh120": mesh120, "hours_cap": hours, "tasks_total": len(tasks)},
        workers=workers,
    )
    (out / "provenance_start.json").write_text(json.dumps(rec, indent=1))
    t_start = time.monotonic()
    deadline = t_start + hours * 3600.0
    n_done = n_err = n_skip = 0
    with Pool(workers) as pool:
        for task, outcome, _ in pool.imap_unordered(fg.work, [(t, str(out), deadline) for t in todo], chunksize=1):
            n_done += outcome == "done"
            n_skip += outcome == "skipped_deadline"
            if outcome.startswith("error"):
                n_err += 1
                print(f"  ERROR {task}: {outcome}", flush=True)
            if n_done and n_done % 100 == 0:
                print(
                    f"[{(time.monotonic() - t_start) / 3600:.2f} h] done {n_done}/{len(todo)}, errors {n_err}",
                    flush=True,
                )
    el = time.monotonic() - t_start
    print(
        f"finished: done {n_done}, errors {n_err}, skipped at the deadline {n_skip}, wall {el / 3600:.2f} h", flush=True
    )
    (out / "provenance_end.json").write_text(
        json.dumps(
            prov.finish(rec)
            | {"tasks_done_this_run": n_done, "errors": n_err, "skipped_deadline": n_skip, "wall_s": el},
            indent=1,
        )
    )


def status(out_dir, mesh120=False):
    out = Path(out_dir)
    for mesh in (120, 90, 60):
        ts = [t for t in build_tasks(mesh120) if t[3] == mesh]
        if ts:
            print(f"mesh {mesh}: {sum((out / 'tasks' / fg.task_name(*t)).exists() for t in ts)}/{len(ts)}")


def consolidate(out_dir):
    out = Path(out_dir)
    tasks = build_tasks(mesh120=True)
    pts = lib.grid_points()
    (out / "chunks").mkdir(exist_ok=True)
    lib.CHUNKS = out / "chunks"
    done = {"grid": [], "corr60": [], "corr90": [], "corr120": []}
    for kind, mesh, points in (
        ("grid", 60, pts),
        ("corr60", 60, lib.CORR_POINTS),
        ("corr90", 90, lib.CORR_POINTS),
        ("corr120", 120, lib.CORR_POINTS),
    ):
        for f in BAND_F:
            files = [out / "tasks" / fg.task_name(f, th, ph, mesh) for th, ph in points]
            if all(p.exists() for p in files):
                lib.save_atomic(
                    lib.chunk_path(kind, f),
                    f_ghz=f,
                    points=np.array(points),
                    s=np.array([np.load(p) for p in files]),
                    mesh=mesh,
                )
                done[kind].append(f)
    h = hashlib.sha256()
    n_tasks = 0
    for p in sorted((out / "tasks").glob("m*.npy")):
        h.update(p.name.encode() + hashlib.sha256(p.read_bytes()).digest())
        n_tasks += 1
    rows = []
    for p in sorted((out / "chunks").glob("*.npz")) + [
        q for q in (out / "provenance_start.json", out / "provenance_end.json") if q.exists()
    ]:
        b = p.read_bytes()
        rows.append({"file": str(p.relative_to(out)), "size": len(b), "sha256": hashlib.sha256(b).hexdigest()})
    manifest = {
        "job": "r4-band-refine", "label": "SERVER_AUTHORED_PENDING_MAC_REVIEW", "data_dir": str(out), "band_f": BAND_F,
        "producing_commit": os.popen(f'git -C "{ROOT}" rev-parse HEAD').read().strip(), "machine": os.environ.get("COMPUTERNAME", "unknown"),
        "chunks_complete": done, "task_files": {"count": n_tasks, "combined_sha256": h.hexdigest()}, "tasks_defined": len(tasks), "files": rows,
    }  # fmt: skip
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    print({k: len(v) for k, v in done.items()}, n_tasks, "task files; MANIFEST.json written")


if __name__ == "__main__":
    cmd, outd, opts = sys.argv[1], sys.argv[2], sys.argv[3:]
    w = int(opts[opts.index("--workers") + 1]) if "--workers" in opts else 40
    hrs = float(opts[opts.index("--hours") + 1]) if "--hours" in opts else 6.0
    m120 = "--mesh120" in opts
    {
        "run": lambda: run(outd, w, hrs, m120),
        "status": lambda: status(outd, m120),
        "consolidate": lambda: consolidate(outd),
    }[cmd]()
