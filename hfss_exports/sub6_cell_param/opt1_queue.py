"""Resumable CPU-sized Option 1 queue; own local AEDT processes only."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import cell_param_pyaedt as cp
from opt1_cell import configure, values
from opt1_geometry import Geometry, all_geometries

configure()
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "lab_B_opt1_work"
RESULTS = REPO / "results" / "sub6_cell_param" / "lab_B_opt1"
GEOMETRIES = [Geometry(**x["geometry"]) for x in all_geometries(0.5)]
CORES_EACH = 8
MAX_PROC = 5


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def cases(g: Geometry) -> list[dict]:
    prefix = "OPT1_" + g.name.replace(".", "p")
    return [{"*": f"{prefix}_{state}_{direction}_{phase}", "state": state,
             "direction": direction, "phase_deg": str(phase)}
            for state in ("OFF", "ON") for direction in ("along", "across")
            for phase in (80, 90, 100)]


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["*", "state", "direction", "phase_deg"])
        writer.writeheader()
        writer.writerows(rows)


def free_cpus() -> float:
    cmd = ["powershell", "-NoProfile", "-Command",
           "(Get-Counter '\\Processor(_Total)\\% Processor Time' -SampleInterval 2 -MaxSamples 4)"
           ".CounterSamples.CookedValue | Measure-Object -Average | % Average"]
    try:
        pct = float(subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True).stdout.strip())
        if not 0 <= pct <= 100:
            raise ValueError(pct)
        return (os.cpu_count() or 1) * (1 - pct / 100)
    except Exception:
        return 0.0


def layout(rows: int) -> tuple[int, float]:
    while True:
        free = free_cpus()
        n = min(rows, MAX_PROC, math.floor(free / CORES_EACH))
        if n:
            return n, free
        print(f"CPU free {free:.1f}; waiting for {CORES_EACH} free logical CPUs", flush=True)
        time.sleep(60)


def validate_result(path: Path, expected: dict, g: Geometry) -> dict:
    result = json.loads(path.read_text())
    if result["run_id"] != expected["*"] or path.stem != expected["*"]:
        raise RuntimeError("Result ID mismatch")
    if result["project_variant"] != g.name or len(result["modes"]) != 8:
        raise RuntimeError(f"Variant/modes mismatch: {path.name}")
    intended = values(g, expected["state"], expected["direction"], int(expected["phase_deg"]))
    if result["variables"] != intended:
        raise RuntimeError(f"Intended variables mismatch: {path.name}")
    good, reason = cp.readback_ok(result)
    if not good:
        raise RuntimeError(f"Readback refused: {path.name}: {reason}")
    return result


def launch(g: Geometry, rows: list[dict], n: int, tag: str) -> tuple[list[int], list[Path]]:
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    root = WORK / "runs" / g.name / f"{timestamp}_{tag}"
    if root.exists():
        raise RuntimeError("Refusing to reuse launch folder")
    base = WORK / "projects" / f"{g.name}.aedt"
    procs = []
    dirs = []
    for i in range(n):
        part_rows = rows[i::n]
        part = root / f"part{i}"
        part.mkdir(parents=True)
        private_project = part / f"{g.name}.aedt"
        shutil.copy2(base, private_project)
        csv_path = part / "rows.csv"
        write_csv(csv_path, part_rows)
        launch_dir = part / "launch"
        launch_dir.mkdir()
        cmd = [sys.executable, str(HERE / "opt1_cell.py"), "--geometry", g.name,
               "--wide-channel", "0.5", "--csv", str(csv_path), "--project", str(private_project),
               "--output", str(part), "--cores", str(CORES_EACH)]
        log = (part / "worker.log").open("w", encoding="utf-8")
        try:
            procs.append(subprocess.Popen(cmd, cwd=launch_dir, stdout=log, stderr=subprocess.STDOUT))
        finally:
            log.close()
        dirs.append(part)
        if i < n - 1:
            time.sleep(40)
    return [p.wait() for p in procs], dirs


def collect(dirs: list[Path], by_id: dict, g: Geometry) -> tuple[dict, dict]:
    solved = {}
    failed = {}
    for part in dirs:
        for file in part.glob("OPT1_*.json"):
            run_id = file.stem
            if run_id not in by_id:
                raise RuntimeError(f"Unexpected row {run_id}")
            validate_result(file, by_id[run_id], g)
            solved[run_id] = file
        failure_file = part / "FAILED_ROWS.json"
        if failure_file.exists():
            for item in json.loads(failure_file.read_text()):
                failed[item["run_id"]] = item
    return solved, failed


def push_results() -> None:
    branch = subprocess.run(["git", "branch", "--show-current"], cwd=REPO,
                            capture_output=True, text=True, check=True).stdout.strip()
    if branch != "compute/sub6-cell-param-lab-b":
        raise RuntimeError(f"Wrong branch {branch!r}")
    rel = RESULTS.relative_to(REPO).as_posix()
    subprocess.run(["git", "add", rel], cwd=REPO, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    staged = subprocess.run(["git", "diff", "--cached", "--name-only"], cwd=REPO,
                            capture_output=True, text=True, check=True).stdout.splitlines()
    if any(not path.startswith(rel + "/") for path in staged):
        raise RuntimeError("Unexpected staged path")
    if staged:
        subprocess.run(["git", "commit", "-m", "SUB6 lab B Option 1 batch results"],
                       cwd=REPO, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["git", "push", "gh", branch], cwd=REPO, check=True, stdout=subprocess.DEVNULL)
        print(f"PUSHED {len(staged)} result files", flush=True)


def run_geometry(g: Geometry, failed_global: dict) -> dict:
    expected_rows = cases(g)
    by_id = {row["*"]: row for row in expected_rows}
    done = {p.stem for p in (RESULTS / "rows").glob("*.json")}
    todo = [row for row in expected_rows if row["*"] not in done and row["*"] not in failed_global]
    if not todo:
        print(f"SKIP {g.name}: 12 rows already resolved", flush=True)
        return {"geometry": g.name, "rows_to_do": 0}
    n, free = layout(len(todo))
    print(f"RUN {g.name}: {len(todo)} rows, {n} processes x {CORES_EACH} cores; free {free:.1f}", flush=True)
    start = time.time()
    codes, dirs = launch(g, todo, n, "main")
    solved, failed = collect(dirs, by_id, g)
    missing = [row for row in todo if row["*"] not in solved and row["*"] not in failed]
    retry_codes = []
    if missing:
        layout(1)  # remeasure before retry; launch one 8-core process only
        print(f"RETRY {g.name}: {len(missing)} unreached rows", flush=True)
        retry_codes, retry_dirs = launch(g, missing, 1, "retry")
        retried, retry_failed = collect(retry_dirs, by_id, g)
        solved.update(retried)
        failed.update(retry_failed)
        for row in missing:
            if row["*"] not in solved and row["*"] not in failed:
                failed[row["*"]] = {"run_id": row["*"], "attempts": 1,
                                      "errors": [f"NOT_SOLVED: process exits {codes}/{retry_codes}"],
                                      "state": row["state"], "direction": row["direction"],
                                      "phase_deg": int(row["phase_deg"]), "geometry": g.name}
    for run_id, source in solved.items():
        target = RESULTS / "rows" / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            validate_result(target, by_id[run_id], g)
        else:
            shutil.copy2(source, target)
    failed_global.update(failed)
    write_json(RESULTS / "FAILED_ROWS.json", list(failed_global.values()))
    report = {"geometry": g.name, "rows_to_do": len(todo), "solved": len(solved),
              "failed": len(failed), "processes": n, "cores_each": CORES_EACH,
              "free_cpus_at_start": round(free, 1), "wall_seconds": round(time.time() - start),
              "exit_codes": codes, "retry_exit_codes": retry_codes}
    write_json(RESULTS / "parallel_reports" / f"{g.name}.json", report)
    if sum(code != 0 for code in codes + retry_codes) > 1:
        raise RuntimeError("Repeated worker crashes: stop for review")
    return report


def main() -> None:
    if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
        raise RuntimeError("Explicit AEDT 2025.2 selection required")
    for gate in (RESULTS / "geometry_precheck.json", RESULTS / "verify" / "determinism.json"):
        if json.loads(gate.read_text())["status"] != "PASS":
            raise RuntimeError(f"Gate failed: {gate.name}")
    if any(not (WORK / "projects" / f"{g.name}.aedt").exists() for g in GEOMETRIES):
        raise RuntimeError("Missing built projects")
    expected = {row["*"] for g in GEOMETRIES for row in cases(g)}
    write_json(RESULTS / "row_manifest.json", sorted(expected))
    failed_file = RESULTS / "FAILED_ROWS.json"
    failed = {x["run_id"]: x for x in json.loads(failed_file.read_text())} if failed_file.exists() else {}
    for index, g in enumerate(GEOMETRIES, 1):
        report = run_geometry(g, failed)
        print(f"DONE {g.name}: {report.get('solved', 0)} solved, {report.get('failed', 0)} failed", flush=True)
        if index % 2 == 0 or index == len(GEOMETRIES):
            push_results()
            attempted = len(list((RESULTS / "rows").glob("*.json"))) + len(failed)
            if attempted and len(failed) / attempted > 0.05:
                raise RuntimeError("FAILED rate >5%; stop for review")
    solved = {p.stem for p in (RESULTS / "rows").glob("*.json")}
    if solved | set(failed) != expected:
        raise RuntimeError(f"Incomplete queue: {len(solved)} solved + {len(failed)} failed != 120")
    print(f"QUEUE COMPLETE: {len(solved)} solved, {len(failed)} FAILED", flush=True)


if __name__ == "__main__":
    main()
