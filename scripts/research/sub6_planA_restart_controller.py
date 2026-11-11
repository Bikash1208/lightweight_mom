"""Stage 1 batch controller. Owns only newly launched AEDT processes."""
import json
import os
import shutil
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[2]
PYTHON = ROOT / ".venv_lab" / "Scripts" / "python.exe"
AEDT = Path(r"C:\Program Files\ANSYS Inc\v252\AnsysEM\ansysedtng.exe")
SOURCE = ROOT / "stage1_work"
WORK = ROOT / "stage1_restart_r1"
RESULTS = ROOT / "results" / "sub6_planA" / "stage1"
PAIRS = [("V", "OFF", "ON"), ("T2", "OFF", "ON"),
         ("T5", "OFF", "ON"), ("T1", "OFF", "ON"),
         ("T8", "OFF", "ON")]
BUILD_VALIDATE = Path(__file__).with_name("sub6_planA_restart_validate.py")
EXTRACT = Path(__file__).with_name("sub6_planA_restart_extract.py")


def utc():
    return datetime.now(timezone.utc).isoformat()


def write(path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, default=str))
    os.replace(tmp, path)


def gate_build(run, variant, state):
    source = SOURCE / f"{run}.aedt"
    report = json.loads((SOURCE / f"{run}_build.json").read_text())
    assert source.is_file() and source.stat().st_size > 10_000_000
    assert report["run_id"] == run and report["state"] == state
    assert report["status"] == "BUILT_UNSOLVED" and report["geometry_gate"] == "PASS"
    assert report["aedt_version"] == "2025.2"
    expected_t = {"V": .8, "T8": .8, "T5": .508, "T2": .254, "T1": .127}[variant]
    assert report["target_bias_thickness_mm"] == expected_t
    assert report["counts_before"] == {"metal_bodies": 1000, "diode_sheets": 240,
                                        "via_bodies": 480, "substrates": 2, "total": 1244}
    assert report["material_readback"] == {
        "SUBSTRATE": "RT5880_STAGE1", "SUBSTRATE_FR4": "vacuum" if variant == "V" else "RT5880_STAGE1",
        "AIR": "vacuum"}
    assert abs(report["bias_after"][2] - (-4.010 + .8 - expected_t)) < 1e-5
    assert len(report["impedance_readback"]["Objects"]) == 240
    rs, xs = ((3.885, -58.70) if state == "OFF" else (.28, 5.11))
    assert float(report["impedance_readback"]["Resistance"]) == rs
    assert float(report["impedance_readback"]["Reactance"]) == xs
    assert report["setup_readback"] == {"Frequency": "5.8GHz", "MaximumPasses": 4, "MaxDeltaS": .02}
    assert report["port_readback"]["Impedance"] == "50ohm"
    return source, report


def run_logged(args, log, env=None, affinity=None, peak=None):
    with log.open("w", encoding="utf-8", errors="replace") as stream:
        stream.write(f"START {utc()} CMD {args!r}\n")
        stream.flush()
        p = subprocess.Popen([str(x) for x in args], cwd=ROOT, env=env,
                             stdout=stream, stderr=subprocess.STDOUT,
                             creationflags=subprocess.CREATE_NO_WINDOW)
        if affinity:
            try:
                psutil.Process(p.pid).cpu_affinity(affinity)
            except psutil.Error as exc:
                raise RuntimeError(f"AFFINITY_SET_FAILED pid={p.pid}: {exc}")
        while p.poll() is None:
            if peak is not None:
                try:
                    tree = [psutil.Process(p.pid)] + psutil.Process(p.pid).children(recursive=True)
                    peak[0] = max(peak[0], sum(q.memory_info().rss for q in tree if q.is_running()))
                except psutil.Error:
                    pass
            time.sleep(5)
        stream.write(f"\nEXIT {utc()} CODE {p.returncode}\n")
        return p.returncode, p.pid


def push_pair(pair_index, jobs):
    def git(*args):
        p = subprocess.run(["git", *args], cwd=ROOT, text=True,
                           capture_output=True, check=False)
        if p.returncode:
            raise RuntimeError(f"GIT_{args[0]}_FAILED code={p.returncode} {p.stderr.strip()}")
        return p.stdout.strip()

    if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=ROOT).returncode:
        raise RuntimeError("UNRELATED_STAGED_CHANGES")
    paths = [str((RESULTS / run).relative_to(ROOT)) for run, *_ in jobs]
    if pair_index == 1:
        paths.extend(["scripts/research/sub6_planA_restart_validate.py",
                      "scripts/research/sub6_planA_restart_controller.py",
                      "scripts/research/sub6_planA_restart_extract.py"])
    git("add", "--", *paths)
    if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=ROOT).returncode:
        git("commit", "-m", f"Stage 1 anchor batch pair {pair_index}")
        git("push", "gh", "compute/sub6-cell-param-lab-b")
    return git("rev-parse", "HEAD")


def main():
    WORK.mkdir(exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    status_path = WORK / "controller_status.json"
    status = {"status": "STARTING", "start_utc": utc(), "controller_pid": os.getpid(),
              "max_concurrent_solves": 2, "pair_order": [v for v, *_ in PAIRS], "runs": {}}
    active = []
    stop = threading.Event()

    def heartbeat():
        while not stop.is_set():
            status["heartbeat_utc"] = utc()
            write(status_path, status)
            stop.wait(60)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        assert PYTHON.is_file() and AEDT.is_file() and EXTRACT.is_file()
        allowed = psutil.Process().cpu_affinity()
        assert len(allowed) >= 2
        halves = [allowed[:len(allowed)//2], allowed[len(allowed)//2:]]
        status["cpu_affinity_per_batch"] = [len(x) for x in halves]
        for pair_index, (variant, *states) in enumerate(PAIRS, 1):
            status["status"] = f"PAIR_{pair_index}_PREP"
            jobs = []
            for state in states:
                run = f"ANCHOR_{variant}_{state}"
                src, build = gate_build(run, variant, state)
                folder = WORK / run
                folder.mkdir(exist_ok=True)
                dst = folder / f"{run}.aedt"
                if not dst.exists():
                    shutil.copy2(src, dst)
                elif dst.stat().st_size != src.stat().st_size:
                    if (folder / "batch.log").exists():
                        raise RuntimeError(f"RESTART_PROJECT_MISMATCH {run}")
                    shutil.copy2(src, dst)
                validation = folder / "validation.json"
                code, _ = run_logged([PYTHON, BUILD_VALIDATE, dst, SOURCE / f"{run}_build.json", validation],
                                     folder / "validate.log")
                if code or json.loads(validation.read_text())["status"] != "PASS":
                    raise RuntimeError(f"VALIDATION_FAILED {run}: {code}")
                result_dir = RESULTS / run
                result_dir.mkdir(exist_ok=True)
                jobs.append((run, dst, folder, result_dir))
                status["runs"][run] = {"status": "VALIDATED", "project": str(dst),
                                       "validation": str(validation)}
            mem = psutil.virtual_memory()
            if mem.available < 64 * 2**30:
                raise RuntimeError(f"MEMORY_HEADROOM_LOW available_gb={mem.available/2**30:.1f}")
            active = []
            status["status"] = f"PAIR_{pair_index}_SOLVING"
            for i, (run, dst, folder, result_dir) in enumerate(jobs):
                log = folder / "batch.log"
                stream = log.open("w", encoding="utf-8", errors="replace")
                stream.write(f"START {utc()} AEDT 2025.2 AFFINITY_CPUS {len(halves[i])}\n")
                stream.flush()
                p = subprocess.Popen([str(AEDT), "-Local", "-monitor", "-batchsolve",
                                      "AnchorPreflight:Nominal:Setup1", str(dst)],
                                     cwd=folder, stdout=stream, stderr=subprocess.STDOUT,
                                     creationflags=subprocess.CREATE_NO_WINDOW)
                psutil.Process(p.pid).cpu_affinity(halves[i])
                active.append((run, p, stream, folder, result_dir, time.monotonic(), [0]))
                status["runs"][run].update(status="SOLVING", batch_pid=p.pid)
            while active:
                for item in list(active):
                    run, p, stream, folder, result_dir, start, peak = item
                    try:
                        tree = [psutil.Process(p.pid)] + psutil.Process(p.pid).children(recursive=True)
                        peak[0] = max(peak[0], sum(q.memory_info().rss for q in tree if q.is_running()))
                    except psutil.Error:
                        pass
                    if p.poll() is None:
                        continue
                    stream.write(f"\nEXIT {utc()} CODE {p.returncode}\n")
                    stream.close()
                    active.remove(item)
                    status["runs"][run].update(batch_exit_code=p.returncode,
                                                batch_seconds=round(time.monotonic()-start, 1),
                                                peak_ram_gb=round(peak[0]/2**30, 3))
                    if p.returncode:
                        raise RuntimeError(f"BATCH_EXIT_FAILED {run} code={p.returncode}")
                    status["runs"][run]["status"] = "BATCH_DONE"
                time.sleep(5) if active else None
            status["status"] = f"PAIR_{pair_index}_EXTRACTING"
            for run, dst, folder, result_dir in jobs:
                report = result_dir / f"{run}.json"
                code, _ = run_logged([PYTHON, EXTRACT, dst, result_dir,
                                      str(status["runs"][run]["peak_ram_gb"]),
                                      str(status["runs"][run]["batch_seconds"])],
                                     folder / "extract.log")
                if code or not report.is_file() or json.loads(report.read_text())["status"] != "COMPLETE":
                    raise RuntimeError(f"EXTRACTION_FAILED {run} code={code}")
                status["runs"][run]["status"] = "COMPLETE"
            if pair_index in (1, 2):
                status["status"] = f"PAIR_{pair_index}_PUSHING"
                status[f"pair_{pair_index}_push_head"] = push_pair(pair_index, jobs)
                status[f"pair_{pair_index}_pushed_utc"] = utc()
            status["status"] = f"PAIR_{pair_index}_COMPLETE"
        status["final_push_head"] = push_pair(5, [
            (f"ANCHOR_{variant}_{state}",) for variant in ("T5", "T1", "T8")
            for state in ("OFF", "ON")])
        status["status"] = "COMPLETE"
    except BaseException:
        status["status"] = "ABORTED"
        status["traceback"] = traceback.format_exc()
        with (WORK / "controller.log").open("a", encoding="utf-8") as log:
            log.write(status["traceback"])
        for run, p, stream, folder, result_dir, start, peak in active:
            if p.poll() is None:
                try:
                    proc = psutil.Process(p.pid)
                    for child in proc.children(recursive=True):
                        child.terminate()
                    proc.terminate()
                except psutil.Error:
                    pass
            try:
                stream.close()
            except OSError:
                pass
        raise
    finally:
        status["end_utc"] = utc()
        write(status_path, status)
        stop.set()
        thread.join(timeout=2)


if __name__ == "__main__":
    main()
