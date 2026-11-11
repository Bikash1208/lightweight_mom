"""Resumable Cell 2 queue; every AEDT worker is a new local instance."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import cell_param_pyaedt as cp
from cell2_cell import configure, values
from opt1_queue import layout, write_csv, write_json
from cell2_spec import GEOMETRIES, cases as cell2_cases, variant as cell2_variant

configure()
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "lab_B_cell2_work"
RESULTS = REPO / "results" / "sub6_cell_param" / "lab_B_cell2"
CORES_EACH = 8
LIDS = {"PEC150": (150.0, "PEC")}


def variant(g, lid):
    if lid != "PEC150":
        raise ValueError(lid)
    return cell2_variant(g)


def cases(g, lid):
    if lid != "PEC150":
        raise ValueError(lid)
    return cell2_cases(g)


def validate_result(file: Path, expected: dict, g, lid: str, height: float) -> None:
    row = json.loads(file.read_text())
    if file.stem != expected["*"] or row["run_id"] != expected["*"]:
        raise RuntimeError("Row ID mismatch")
    if row["project_variant"] != variant(g, lid) or len(row["modes"]) != 8:
        raise RuntimeError(f"Variant/modes mismatch {file.name}")
    intended = values(g, expected["state"], expected["direction"], int(expected["phase_deg"]), height)
    good, reason = cp.readback_ok(row)
    if row["variables"] != intended or not good:
        raise RuntimeError(f"Readback refused {file.name}: {reason}")


def launch(g, lid: str, height: float, kind: str, rows: list[dict], n: int, tag: str):
    name = variant(g, lid)
    root = WORK / "runs" / name / f"{time.strftime('%Y%m%d_%H%M%S')}_{tag}"
    if root.exists():
        raise RuntimeError("Refusing to reuse launch folder")
    base = WORK / "projects" / f"{name}.aedt"
    procs, dirs = [], []
    for i in range(n):
        part = root / f"part{i}"
        part.mkdir(parents=True)
        project = part / f"{name}.aedt"
        shutil.copy2(base, project)
        row_csv = part / "rows.csv"
        write_csv(row_csv, rows[i::n])
        launch_dir = part / "launch"
        launch_dir.mkdir()
        cmd = [sys.executable, str(HERE / "cell2_cell.py"), "--geometry", g.name,
               "--lid-height-mm", str(height),
               "--top-boundary", kind, "--project-variant", name,
               "--csv", str(row_csv), "--project", str(project),
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


def collect(dirs, by_id, g, lid, height):
    solved, failed = {}, {}
    for part in dirs:
        for file in part.glob("CELL2_*.json"):
            if file.stem not in by_id:
                raise RuntimeError(f"Unexpected row {file.name}")
            validate_result(file, by_id[file.stem], g, lid, height)
            solved[file.stem] = file
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
    subprocess.run(["git", "add", rel], cwd=REPO, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    staged = subprocess.run(["git", "diff", "--cached", "--name-only"], cwd=REPO,
                            capture_output=True, text=True, check=True).stdout.splitlines()
    if any(not path.startswith(rel + "/") for path in staged):
        raise RuntimeError("Unexpected staged file")
    if staged:
        subprocess.run(["git", "commit", "-m", "SUB6 lab B Cell 2 batch results"],
                       cwd=REPO, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["git", "push", "gh", branch], cwd=REPO,
                       check=True, stdout=subprocess.DEVNULL)
        print(f"PUSHED {len(staged)} result files", flush=True)


def run_variant(g, lid: str, height: float, kind: str, failed_global: dict) -> dict:
    name = variant(g, lid)
    expected = cases(g, lid)
    by_id = {r["*"]: r for r in expected}
    done = {p.stem for p in (RESULTS / "rows").glob("*.json")}
    todo = [r for r in expected if r["*"] not in done and r["*"] not in failed_global]
    if not todo:
        print(f"SKIP {name}: 28 rows resolved", flush=True)
        return {"variant": name, "rows_to_do": 0}
    n, free = layout(len(todo))
    print(f"RUN {name}: {len(todo)} rows, {n} processes x {CORES_EACH} cores; free {free:.1f}", flush=True)
    start = time.time()
    codes, dirs = launch(g, lid, height, kind, todo, n, "main")
    solved, failed = collect(dirs, by_id, g, lid, height)
    missing = [r for r in todo if r["*"] not in solved and r["*"] not in failed]
    retry_codes = []
    if missing:
        layout(1)
        print(f"RETRY {name}: {len(missing)} unreached rows", flush=True)
        retry_codes, retry_dirs = launch(g, lid, height, kind, missing, 1, "retry")
        retried, retry_failed = collect(retry_dirs, by_id, g, lid, height)
        solved.update(retried)
        failed.update(retry_failed)
        for r in missing:
            if r["*"] not in solved and r["*"] not in failed:
                failed[r["*"]] = {"run_id": r["*"], "attempts": 1,
                                   "errors": [f"NOT_SOLVED: exits {codes}/{retry_codes}"],
                                   "state": r["state"], "direction": r["direction"],
                                   "phase_deg": int(r["phase_deg"]), "variant": name}
    for run_id, source in solved.items():
        target = RESULTS / "rows" / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            validate_result(target, by_id[run_id], g, lid, height)
        else:
            shutil.copy2(source, target)
    failed_global.update(failed)
    write_json(RESULTS / "FAILED_ROWS.json", list(failed_global.values()))
    report = {"variant": name, "rows_to_do": len(todo), "solved": len(solved),
              "failed": len(failed), "processes": n, "cores_each": CORES_EACH,
              "free_cpus_at_start": round(free, 1), "wall_seconds": round(time.time() - start),
              "exit_codes": codes, "retry_exit_codes": retry_codes}
    write_json(RESULTS / "parallel_reports" / f"{name}.json", report)
    if sum(code != 0 for code in codes + retry_codes) > 1:
        raise RuntimeError("Repeated worker crashes; stop for review")
    return report


def main() -> None:
    if os.environ.get("SUB6_AEDT_VERSION") != "2025.2":
        raise RuntimeError("Explicit AEDT 2025.2 selection required")
    stop = RESULTS / "status.json"
    if stop.exists() and json.loads(stop.read_text()).get("status", "").startswith("STOPPED"):
        raise RuntimeError("Cell 2 hard stop recorded; refusing additional AEDT solves")
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", choices=[g.name for g in GEOMETRIES])
    args = ap.parse_args()
    height, kind = LIDS["PEC150"]
    matrix = [(g, "PEC150", height, kind) for g in GEOMETRIES
              if args.geometry is None or g.name == args.geometry]
    for g, lid, height, kind in matrix:
        summary = json.loads((RESULTS / "builds" / f"{variant(g,lid)}.json").read_text())
        if summary["status"] != "BUILT" or summary["top_boundary"]["type"] != kind:
            raise RuntimeError("Build gate failed")
    expected = {r["*"] for g, lid, _, _ in matrix for r in cases(g, lid)}
    manifest = RESULTS / "row_manifest.json"
    previous = set(json.loads(manifest.read_text())) if manifest.exists() else set()
    write_json(manifest, sorted(previous | expected))
    failure_file = RESULTS / "FAILED_ROWS.json"
    failed = {x["run_id"]: x for x in json.loads(failure_file.read_text())} if failure_file.exists() else {}
    for g, lid, height, kind in matrix:
        report = run_variant(g, lid, height, kind, failed)
        print(f"DONE {variant(g,lid)}: {report.get('solved',0)} solved, {report.get('failed',0)} FAILED", flush=True)
        push_results()
        attempted = len(list((RESULTS / "rows").glob("*.json"))) + len(failed)
        if attempted and len(failed) / attempted > 0.05:
            raise RuntimeError("FAILED rate >5%; stop for review")
    solved = {p.stem for p in (RESULTS / "rows").glob("*.json")}
    if not expected.issubset(solved | set(failed)):
        raise RuntimeError(f"Incomplete queue: {len(solved)} solved + {len(failed)} FAILED")
    print(f"QUEUE COMPLETE: {len(expected)} requested rows, {len(solved)} total solved, {len(failed)} FAILED", flush=True)


if __name__ == "__main__":
    main()
