"""Run SCREEN/SENS/FULL CSVs as several parallel AEDT processes (rows split round-robin).

Layout (processes x cores each) is re-measured before EVERY csv file:
  * free CPU  = logical CPUs - (busy fraction * logical CPUs), measured with a short Get-Counter sample
  * licences  = A/B/C/D counts that the first process's Ansys client log writes at checkout
                (B = used by us, C = used by all users, D = pool available). lmstat is unusable here.
Each process owns a private copy of the built project (AEDT locks projects) and its own output part dir.
All gates (VERIFY, artefact test, absorbing floor, diode profile) are enforced by cell_param_pyaedt.py
before AEDT starts. Nothing here caps cores below what is free.

python run_parallel.py --project BUILT.aedt --build-output screen_bar90 --verify-gate verify.json \
       --csv SCREEN_1.csv SCREEN_2.csv SCREEN_3.csv --out screen_bar90 [--cores-each 16]
"""
import argparse
import csv
import glob
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
QUEUE = HERE / "screen_queue.json"   # optional: {"override": {"<out dir name>": {"skip": true} | {"jobs": [{project, build_output, out, csv[]}]}}}
PY = sys.executable
ANSYS_LOGS = Path(os.environ.get("TEMP", "")) / ".ansys"
LOGICAL = os.cpu_count() or 1


def busy_cpus():
    cmd = ["powershell", "-NoProfile", "-Command",
           "(Get-Counter '\\Processor(_Total)\\% Processor Time' -SampleInterval 2 -MaxSamples 4)"
           ".CounterSamples.CookedValue | Measure-Object -Average | % Average"]
    try:
        pct = float(subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip())
    except Exception:
        pct = 100.0  # unknown -> assume busy, never over-subscribe silently
    return pct / 100 * LOGICAL


def licence_counts(since):
    """Latest A/B/C/D per feature from the newest Ansys client log written after `since` (epoch s).
    A=requested, B=used by us, C=used by all users, D=pool available. (lmstat is unusable on this server.)"""
    logs = [Path(p) for p in glob.glob(str(ANSYS_LOGS / "ansyscl.*.*.log")) if Path(p).stat().st_mtime >= since]
    out = {}
    for path in sorted(logs, key=lambda q: q.stat().st_mtime):
        for line in path.read_text(errors="replace").splitlines():
            m = re.search(r"(?:SPLIT_)?CHECKOUT\s+(\S+)\s+\S+\s+\(\S+\)\s+(\d+)/(\d+)/(\d+)/(\d+)", line)
            if m:
                out[m.group(1)] = tuple(int(x) for x in m.groups()[1:])
    return out


def plan(rows, cores_each, max_proc, fixed=False):
    free = LOGICAL - busy_cpus()
    n = max(1, min(rows, max_proc, int(free // cores_each)))
    # Use ALL free cores: spread the remainder over the processes (never below the requested cores each).
    return n, (cores_each if fixed else max(cores_each, int(free // n))), free


def split(csv_path, n, dest):
    rows = list(csv.DictReader(csv_path.open(newline="", encoding="utf-8")))
    parts = []
    for k in range(n):
        p = dest / f"{csv_path.stem}_part{k}.csv"
        with p.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows[k::n])
        parts.append(p)
    return parts, len(rows)


def is_result(js):
    try:
        return "modes" in json.loads(js.read_text())
    except Exception:
        return False


def launch_parts(a, csv_path, rows_csv, n, cores, work, tag):
    """Run `rows_csv` as n parallel solver processes; returns (exit codes, part dirs)."""
    parts, _ = split(rows_csv, n, work)
    procs = []; pdirs = []
    for k, part in enumerate(parts):
        pdir = work / f"{tag}part{k}"; pdir.mkdir(exist_ok=True); pdirs.append(pdir)
        shutil.copy(a.build_output / "build_summary.json", pdir / "build_summary.json")
        proj = work / f"{a.project.stem}_{tag}p{k}.aedt"
        if not proj.exists():
            shutil.copy(a.project, proj)          # geometry/setup only; results are not copied
        ld = Path(os.environ.get("SUB6_LAUNCH_ROOT", str(HERE / "aedt_launch"))) / f"aedt_launch_{csv_path.stem}_{tag}p{k}"
        ld.mkdir(parents=True, exist_ok=True)
        log = (work / f"{tag}part{k}.log").open("w")
        procs.append(subprocess.Popen(
            [PY, str(HERE / "cell_param_pyaedt.py"), "--sweep", str(part), "--diode", a.diode,
             "--project", str(proj), "--output", str(pdir), "--cores", str(cores),
             "--row-timeout-min", str(a.row_timeout_min), "--verify-gate", str(a.verify_gate)],
            cwd=str(ld), stdout=log, stderr=subprocess.STDOUT))
        time.sleep(a.stagger)   # stagger: PyAEDT validation/launch is not concurrency-safe
        if k == 0 and tag == "":
            lic = licence_counts(time.time() - a.stagger - 60)
            print("  licences (A requested / B ours / C all users / D pool):", lic, flush=True)
            print("  held by OTHER users:", {f: c - b for f, (a_, b, c, d) in lic.items()}, flush=True)
    return [p.wait() for p in procs], pdirs


def run_job(a):
    a.out.mkdir(parents=True, exist_ok=True)
    report = []
    sys.path.insert(0, str(HERE))
    import cell_param_pyaedt as cp
    cp.set_profile(a.diode)
    for csv_path in a.csv:
        csv_path = csv_path if csv_path.is_absolute() else HERE / csv_path
        all_rows = list(csv.DictReader(csv_path.open(newline="", encoding="utf-8")))
        solved = {js.stem for js in a.out.glob("*.json") if is_result(js)}
        ffp = a.out / f"FAILED_ROWS_{csv_path.stem}.json"
        logged_failed = {x["run_id"] for x in json.loads(ffp.read_text())} if ffp.exists() else set()   # failed twice already: not retried on resume
        todo = [r for r in all_rows if r["*"] not in solved and r["*"] not in logged_failed]
        print(f"[{time.strftime('%H:%M:%S')}] {csv_path.name}: {len(all_rows)} rows, {len(all_rows)-len(todo)} already solved, {len(todo)} to do", flush=True)
        if not todo:
            continue
        work = a.out / f"_{csv_path.stem}_{time.strftime('%m%d%H%M%S')}"
        work.mkdir(exist_ok=True)
        todo_csv = work / f"{csv_path.stem}_todo.csv"
        with todo_csv.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0])); w.writeheader(); w.writerows(todo)
        if a.layout:
            n, cores = (int(x) for x in a.layout.lower().split("x")); n = min(n, len(todo)); free = LOGICAL - busy_cpus()
        else:
            n, cores, free = plan(len(todo), a.cores_each, a.max_proc, a.fixed_cores)
        print(f"[{time.strftime('%H:%M:%S')}] layout {n} processes x {cores} cores (free CPUs {free:.1f})", flush=True)
        t0 = time.time()
        codes, pdirs = launch_parts(a, csv_path, todo_csv, n, cores, work, "")

        def collect(pdirs_):
            got, failed_ = [], []
            for pdir in pdirs_:
                for js in sorted(pdir.glob("*.json")):
                    if js.name != "FAILED_ROWS.json" and is_result(js):
                        r = json.loads(js.read_text()); ok, why = cp.readback_ok(r)
                        if ok:
                            shutil.copy(js, a.out / js.name); got.append(r)
                        else:
                            failed_.append({"run_id": r["run_id"], "attempts": 0, "errors": ["READBACK_REFUSED: " + why]})
                fr = pdir / "FAILED_ROWS.json"
                if fr.exists():
                    failed_ += [{k: v for k, v in x.items() if k != "variables"} for x in json.loads(fr.read_text())]
            return got, failed_

        merged, failed = collect(pdirs)
        done_ids = {r["run_id"] for r in merged} | {x["run_id"] for x in failed}
        missing = [r for r in todo if r["*"] not in done_ids]
        codes2 = []
        if missing:   # a process died before reaching these rows: one more pass, single process, then FAILED
            print(f"  {len(missing)} rows not reached (exit codes {codes}); one retry pass", flush=True)
            rcsv = work / f"{csv_path.stem}_retry.csv"
            with rcsv.open("w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(all_rows[0])); w.writeheader(); w.writerows(missing)
            codes2, pdirs2 = launch_parts(a, csv_path, rcsv, 1, max(cores, 13), work, "retry_")
            m2, f2 = collect(pdirs2); merged += m2; failed += f2
            done_ids = {r["run_id"] for r in merged} | {x["run_id"] for x in failed}
            for r in todo:
                if r["*"] not in done_ids:
                    failed.append({"run_id": r["*"], "attempts": 0, "errors": [f"NOT_SOLVED: process exits {codes} / {codes2}"]})
        allres = [json.loads(js.read_text()) for js in sorted(a.out.glob("*.json"))
                  if js.stem.startswith(csv_path.stem + "_") and is_result(js)]
        cp.export_modes(a.out / f"{csv_path.stem}_modes.csv", allres)
        ff = a.out / f"FAILED_ROWS_{csv_path.stem}.json"
        prev = json.loads(ff.read_text()) if ff.exists() else []
        ff.write_text(json.dumps([x for x in prev if x["run_id"] not in {r["run_id"] for r in merged}] + failed, indent=2))
        secs = time.time() - t0
        info = {"file": csv_path.name, "rows_total": len(all_rows), "rows_solved_this_run": len(merged), "rows_failed_this_run": len(failed),
                "processes": n, "cores_each": cores, "free_cpus_at_start": round(free, 1), "wall_seconds": round(secs),
                "mean_row_seconds": round(sum(r["seconds"] for r in merged) / max(len(merged), 1), 1), "exit_codes": codes}
        report.append(info); print("  ", json.dumps(info), flush=True)
        (a.out / "parallel_report.json").write_text(json.dumps(report, indent=2))
        # NEVER stop the queue because of a failed row: failures are listed in FAILED_ROWS_*.json and reported at the end.


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", type=Path, required=True, help="built base project (copied per process)")
    ap.add_argument("--build-output", type=Path, required=True, help="dir holding the base build_summary.json")
    ap.add_argument("--verify-gate", type=Path, required=True)
    ap.add_argument("--csv", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--diode", default="bar90")
    ap.add_argument("--cores-each", type=int, default=16)
    ap.add_argument("--max-proc", type=int, default=8)
    ap.add_argument("--row-timeout-min", type=float, default=30)
    ap.add_argument("--layout", default=None, help="force NxC (processes x cores), e.g. 8x6; default = measured from free CPU")
    ap.add_argument("--stagger", type=float, default=75, help="seconds between process launches")
    ap.add_argument("--fixed-cores", action="store_true", help="use exactly --cores-each per process (as many processes as free CPU allows)")
    a = ap.parse_args()
    a = ap.parse_args()
    ov = None
    if QUEUE.exists():
        ov = json.loads(QUEUE.read_text()).get("override", {}).get(a.out.name)
    if ov and ov.get("skip"):
        print(f"[{time.strftime('%H:%M:%S')}] screen_queue.json: skip invocation for {a.out.name} (its work is scheduled in another job)", flush=True)
        return
    jobs = [a]
    if ov:
        import argparse as _ap
        jobs = []
        for j in ov["jobs"]:
            ns = _ap.Namespace(**vars(a))
            ns.project = Path(j["project"]); ns.build_output = Path(j["build_output"]); ns.out = Path(j["out"])
            ns.csv = [Path(c) for c in j["csv"]]
            jobs.append(ns)
        print(f"[{time.strftime('%H:%M:%S')}] screen_queue.json override for {a.out.name}: "
              f"{[ (j.out.name, [c.name for c in j.csv]) for j in jobs]}", flush=True)
    for j in jobs:
        run_job(j)


if __name__ == "__main__":
    main()
