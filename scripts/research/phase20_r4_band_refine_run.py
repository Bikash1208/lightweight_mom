"""Phase 20 R4 band refinement (SERVER JOB `r4-band-refine`): the library's resonance band 10.4-10.9 GHz gets nodes every 0.05 GHz, in the table AND in the mesh-correction map.

Stages (run in this order; every stage is checkpointed per task, a rerun resumes):
  grid     the table nodes NOT on the 0.1 GHz grid: 10.45, 10.55, 10.65, 10.75, 10.85 GHz, theta 0-45 step 5, phi 0-345 step 15, four states, mesh 60, F = 2n (1085 tasks)
  corr90   correction nodes at every 0.05 GHz step in 10.4-10.9 GHz that has none yet (10.40, 10.45, 10.55, 10.60, 10.65, 10.70, 10.75, 10.80, 10.85, 10.90), the 5 R4 correction
           points each (theta, phi) = (0,0), (30,0), (30,90), (45,0), (45,90), mesh 90 (50 tasks)
  corr120  the same 50 points at mesh 120 (50 tasks, ~9 GB each: use <= 16 workers)
The mesh-60 correction values at those nodes are table points (the new grid tasks, or the existing v2 grid chunks for 10.40/10.60/10.70/10.80/10.90): no separate tasks.
The wall-clock cap is for the WHOLE job: the first `run` writes <out_dir>/job_start.txt and every later `run` stops at job_start + hours. Same solver call as the library build. MODEL_ONLY.

    python scripts/research/phase20_r4_band_refine_run.py run <out_dir> --stage grid|corr90|corr120 [--workers 40] [--hours 8]
    python scripts/research/phase20_r4_band_refine_run.py status <out_dir>
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
CORR_NODES = [10.40, 10.45, 10.55, 10.60, 10.65, 10.70, 10.75, 10.80, 10.85, 10.90]
STAGES = ("grid", "corr90", "corr120")
fg.MARGIN_S[120] = 5400.0  # a mesh-120 task takes ~1 h: do not start one that cannot finish before the stop time


def build_tasks(stage=None):
    """[(f, theta, phi, mesh)] of one stage (or all stages in order grid, corr90, corr120)."""
    lib.configure_grid(**fg.LIB_STEPS, f_list=BAND_F, corr_f_list=CORR_NODES)
    pts = lib.grid_points()
    by = {
        "grid": [(f, th, ph, 60) for f in BAND_F for th, ph in pts],
        "corr90": [(f, th, ph, 90) for f in CORR_NODES for th, ph in lib.CORR_POINTS],
        "corr120": [(f, th, ph, 120) for f in CORR_NODES for th, ph in lib.CORR_POINTS],
    }
    return by[stage] if stage else sum((by[k] for k in STAGES), [])


def job_deadline(out, hours):
    """Monotonic deadline for the whole job (first run writes job_start.txt as epoch seconds)."""
    f = out / "job_start.txt"
    if not f.exists():
        f.write_text(str(time.time()))
    return time.monotonic() + (float(f.read_text()) + hours * 3600.0 - time.time())


def run(out_dir, workers, hours, stage):
    import provenance as prov

    out = Path(out_dir)
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    tasks = build_tasks(stage)
    todo = [t for t in tasks if not (out / "tasks" / fg.task_name(*t)).exists()]
    deadline = job_deadline(out, hours)
    print(
        f"stage {stage}: {len(tasks)} tasks, {len(tasks) - len(todo)} done, {len(todo)} to do, {workers} workers, {max(deadline - time.monotonic(), 0) / 3600:.2f} h left of the {hours} h cap",
        flush=True,
    )
    rec = prov.provenance(
        {"job": "r4-band-refine", "stage": stage, "hours_cap_whole_job": hours, "tasks_total": len(tasks)},
        workers=workers,
    )
    (out / f"provenance_start_{stage}.json").write_text(json.dumps(rec, indent=1))
    t_start = time.monotonic()
    n_done = n_err = n_skip = 0
    with Pool(workers) as pool:
        for task, outcome, _ in pool.imap_unordered(fg.work, [(t, str(out), deadline) for t in todo], chunksize=1):
            n_done += outcome == "done"
            n_skip += outcome == "skipped_deadline"
            if outcome.startswith("error"):
                n_err += 1
                print(f"  ERROR {task}: {outcome}", flush=True)
            if n_done and n_done % 50 == 0:
                print(
                    f"[{(time.monotonic() - t_start) / 3600:.2f} h] done {n_done}/{len(todo)}, errors {n_err}",
                    flush=True,
                )
    el = time.monotonic() - t_start
    print(
        f"stage {stage} finished: done {n_done}, errors {n_err}, skipped at the deadline {n_skip}, wall {el / 3600:.2f} h",
        flush=True,
    )
    (out / f"provenance_end_{stage}.json").write_text(
        json.dumps(
            prov.finish(rec)
            | {"tasks_done_this_run": n_done, "errors": n_err, "skipped_deadline": n_skip, "wall_s": el},
            indent=1,
        )
    )


def status(out_dir):
    out = Path(out_dir)
    for stage in STAGES:
        ts = build_tasks(stage)
        print(f"{stage}: {sum((out / 'tasks' / fg.task_name(*t)).exists() for t in ts)}/{len(ts)}")


def consolidate(out_dir):
    out = Path(out_dir)
    tasks = build_tasks()
    pts = lib.grid_points()
    (out / "chunks").mkdir(exist_ok=True)
    lib.CHUNKS = out / "chunks"
    done = {"grid": [], "corr90": [], "corr120": []}
    for kind, mesh, points, freqs in (
        ("grid", 60, pts, BAND_F),
        ("corr90", 90, lib.CORR_POINTS, CORR_NODES),
        ("corr120", 120, lib.CORR_POINTS, CORR_NODES),
    ):
        for f in freqs:
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
    for p in sorted((out / "chunks").glob("*.npz")) + sorted(out.glob("provenance_*.json")):
        b = p.read_bytes()
        rows.append({"file": str(p.relative_to(out)), "size": len(b), "sha256": hashlib.sha256(b).hexdigest()})
    manifest = {
        "job": "r4-band-refine", "label": "SERVER_AUTHORED_PENDING_MAC_REVIEW", "data_dir": str(out), "band_f": BAND_F, "corr_nodes": CORR_NODES,
        "producing_commit": os.popen(f'git -C "{ROOT}" rev-parse HEAD').read().strip(), "machine": os.environ.get("COMPUTERNAME", "unknown"),
        "chunks_complete": done, "task_files": {"count": n_tasks, "combined_sha256": h.hexdigest()}, "tasks_defined": len(tasks), "files": rows,
    }  # fmt: skip
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    print({k: len(v) for k, v in done.items()}, n_tasks, "task files; MANIFEST.json written")


if __name__ == "__main__":
    cmd, outd, opts = sys.argv[1], sys.argv[2], sys.argv[3:]
    w = int(opts[opts.index("--workers") + 1]) if "--workers" in opts else 40
    hrs = float(opts[opts.index("--hours") + 1]) if "--hours" in opts else 8.0
    stg = opts[opts.index("--stage") + 1] if "--stage" in opts else None
    if cmd == "run" and stg not in STAGES:
        sys.exit(f"--stage must be one of {STAGES}")
    {"run": lambda: run(outd, w, hrs, stg), "status": lambda: status(outd), "consolidate": lambda: consolidate(outd)}[
        cmd
    ]()
