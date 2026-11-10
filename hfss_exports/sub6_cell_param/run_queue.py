"""QUEUE_A (this machine). QUEUE_B (SCREEN_2 + short islands) is run by the lab machine: see REMOTE_RUN.md.

1. finish SCREEN_4_nostub (rows lost when part 2 crashed)
2. benchmark ONE batch of 8 real rows: 8 processes x 6 cores vs 4 processes x 13 cores; keep the faster layout
3. nostub 80/100 supplement -> trace with the winning core size, all free CPU used (processes = free CPU // cores)
STOP before SCREEN_1 / SCREEN_3 (decision session); SCREEN_2 and the short islands belong to QUEUE_B. A failed row is retried once, logged FAILED and skipped; all FAILED rows are
listed at the end in screen_FAILED_ROWS.json.
"""
import csv
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PY = sys.executable
GATE = HERE / "verify_cal_nomesh_b" / "verify.json"
NOSTUB = dict(project=HERE / "work_screen" / "CellParam_BAR90_nostub.aedt", build=HERE / "screen_bar90_nostub", out=HERE / "screen_bar90_nostub_run")
TRACE = dict(project=HERE / "work_screen" / "CellParam_BAR90_trace.aedt", build=HERE / "screen_bar90_trace", out=HERE / "screen_bar90_trace_run")
FULL = dict(project=HERE / "work_screen" / "CellParam_BAR90.aedt", build=HERE / "screen_bar90", out=HERE / "screen_bar90_full_run")
LOG = HERE / "run_queue.log"


def say(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def rp(job, csvs, *extra, out=None):
    cmd = [PY, "-X", "utf8", str(HERE / "run_parallel.py"), "--project", str(job["project"]), "--build-output", str(job["build"]),
           "--verify-gate", str(GATE), "--csv", *[str(c) for c in csvs], "--out", str(out or job["out"]), "--diode", "bar90",
           "--row-timeout-min", "30", *extra]
    say("run_parallel " + " ".join(extra) + " :: " + ", ".join(Path(str(c)).name for c in csvs))
    return subprocess.run(cmd, cwd=str(HERE)).returncode


def last_wall(out, csv_name):
    rep = json.loads((out / "parallel_report.json").read_text())
    return next(r for r in reversed(rep) if r["file"] == csv_name)


def wait_for_pids(pids, timeout_min=120):
    """Earlier solver runs must be gone (run_parallel merges results before it exits)."""
    t0 = time.time()
    for pid in pids:
        while time.time() - t0 < timeout_min * 60:
            r = subprocess.run(["tasklist", "/FI", f"PID eq {pid}"], capture_output=True, text=True).stdout
            if str(pid) not in r:
                break
            time.sleep(30)
        else:
            raise SystemExit(f"process {pid} still running after {timeout_min} min")


def main():
    say("QUEUE START (local only)")
    wait_for_pids([43300, 22636, 22208])   # 22208 = the real run_parallel of the SCREEN_4 step started by the previous queue driver
    say("earlier runs finished; their results are merged")
    # 1) rows lost to the part-2 crash
    rp(NOSTUB, ["SCREEN_4_nostub.csv"], "--cores-each", "13", "--max-proc", "8", "--stagger", "40")
    # 2) benchmark on 8 real rows = first 8 rows of the nostub 80/100 supplement (nominal geometry, both states/directions)
    rows = list(csv.DictReader((HERE / "SCREEN_4_nostub_extra.csv").open(newline="", encoding="utf-8")))[:8]
    bench = HERE / "BENCH8.csv"
    with bench.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    outB = HERE / "screen_bar90_bench_4x13"
    rp(NOSTUB, [bench], "--layout", "8x6", "--stagger", "30")                       # real results are kept (rows are needed anyway)
    A = last_wall(NOSTUB["out"], "BENCH8.csv")
    rp(NOSTUB, [bench], "--layout", "4x13", "--stagger", "30", out=outB)          # same 8 rows, discarded after timing
    B = last_wall(outB, "BENCH8.csv")
    # fairness: A and B both measured on the same 8 rows, same stagger, same contention window (back to back)
    winner = "8x6" if A["wall_seconds"] <= B["wall_seconds"] else "4x13"
    bench_report = {"8x6": A, "4x13": B, "winner": winner,
                    "note": "wall time for 8 rows incl. process launch stagger 30 s and AEDT start-up; free CPU at start in each entry"}
    (HERE / "layout_benchmark.json").write_text(json.dumps(bench_report, indent=2))
    say(f"BENCHMARK 8x6 wall {A['wall_seconds']} s (free CPUs {A['free_cpus_at_start']}) vs 4x13 wall {B['wall_seconds']} s "
        f"(free CPUs {B['free_cpus_at_start']}) -> keep {winner}")
    cores = winner.split("x")[1]
    flags = ["--cores-each", cores, "--fixed-cores", "--max-proc", "16", "--stagger", "40"]
    # 3) the queue up to the stub/trace comparison; SCREEN_1 and SCREEN_3 are NOT run
    rp(NOSTUB, ["SCREEN_4_nostub_extra.csv"], *flags)
    rp(TRACE, ["SCREEN_5_trace.csv"], *flags)
    # SCREEN_2 (stub 0.3/1/2) and the short-island rows (SCREEN_6_*) are QUEUE_B: they run on the lab machine (run_queue_b.py)
    # and are merged back with merge_lab_b.py through the same read-back gate. They are NOT run here, so nothing runs twice.
    # failed rows, reported at the end (by geometry and state; flags a possible physics/setup issue)
    sys.path.insert(0, str(HERE))
    allfailed = []
    for job in (NOSTUB, TRACE, FULL):
        for ff in sorted(job["out"].glob("FAILED_ROWS_*.json")):
            allfailed += [{"file": ff.stem.replace("FAILED_ROWS_", ""), **x} for x in json.loads(ff.read_text())]
    (HERE / "screen_FAILED_ROWS.json").write_text(json.dumps(allfailed, indent=2))
    try:
        import failed_rows_report as fr
        o = fr.write_report([NOSTUB["out"], TRACE["out"], FULL["out"]], HERE / "screen_failed_report")
        say(f"FAILED ROWS: {o['rows_failed']} of {o['rows_attempted']} ({o['failure_rate_pct']} %); flags: {o['flags'] or 'none'} -> screen_failed_report/FAILED_ROWS_REPORT.md")
    except Exception as exc:
        say(f"failed-row report error: {exc!r}")
    try:
        subprocess.run([PY, "-X", "utf8", str(HERE / "report_screen.py"), "--runs", str(NOSTUB["out"]), str(TRACE["out"]), str(FULL["out"]),
                        "--out", str(HERE / "screen_bar90_report"), "--diode", "bar90"], cwd=str(HERE), check=False)
    except Exception as exc:
        say(f"report_screen error: {exc!r}")
    say(f"QUEUE DONE. Held: SCREEN_1, SCREEN_3 (decision session).")

if __name__ == "__main__":
    main()
