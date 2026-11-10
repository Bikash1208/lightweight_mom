"""QUEUE_B driver (lab machine, LOCAL only; no connection to the main server). See REMOTE_RUN.md.

QUEUE_B = SCREEN_2.csv (stub 0.3 / 1 / 2 mm + fr4_t, disc_r, via_drill, pad_dia, full-cell project) and the short-island rows
SCREEN_6_short_nostub.csv / SCREEN_6_short_trace.csv (island_l 3.5/3.0/2.5/2.0 x island_w 1.3/1.0, no-stub and trace projects),
all at 80/90/100 deg, OFF and ON, BAR90, 377-ohm impedance floor at -34.01 mm, NO mesh operations.

Steps (python run_queue_b.py --step <name>; default 'all' runs env -> verify -> build -> queue -> collect and STOPS if VERIFY fails):
  env        record machine/software facts (no host name)         -> lab_B_work/lab_B_machine.json
  verify     build + 4-row VERIFY on the BAR50 calibration cell, then compare with verify_reference/ (f within 0.5 %)
  build      build the three BAR90 projects (full, nostub, trace)
  benchmark  OPTIONAL: 8 real rows as 8x6 vs 4x13 cores; prints the faster layout (then pass --cores-each to queue)
  queue      run QUEUE_B through run_parallel.py (read-back gate, retry-once, FAILED rows logged and skipped)
  collect    copy rows / reports / FAILED lists to <repo>/results/sub6_cell_param/lab_B/ (you then git add/commit/push, see REMOTE_RUN.md)
"""
import argparse
import csv
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PY = sys.executable
WORK = HERE / "lab_B_work"
PROJ = WORK / "projects"
RESULTS = HERE.parents[1] / "results" / "sub6_cell_param" / "lab_B"
GATE = WORK / "verify" / "verify.json"
LOG = WORK / "run_queue_b.log"

# (project file stem, build variant, build output dir, run out dir, CSV files)
JOBS = [
    dict(name="full", variant="full", project=PROJ / "CellParam_BAR90.aedt", build=WORK / "build_full", out=WORK / "run_full", csv=["SCREEN_2.csv"]),
    dict(name="nostub", variant="nostub", project=PROJ / "CellParam_BAR90_nostub.aedt", build=WORK / "build_nostub", out=WORK / "run_nostub",
         csv=["SCREEN_6_short_nostub.csv"]),
    dict(name="trace", variant="trace", project=PROJ / "CellParam_BAR90_trace.aedt", build=WORK / "build_trace", out=WORK / "run_trace",
         csv=["SCREEN_6_short_trace.csv"]),
]


def say(msg):
    WORK.mkdir(parents=True, exist_ok=True)
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def run(cmd, cwd=None, env_extra=None):
    env = dict(os.environ, **(env_extra or {}))
    return subprocess.run([str(c) for c in cmd], cwd=str(cwd or HERE), env=env).returncode


def launch_env():
    d = WORK / "aedt_launch"; d.mkdir(parents=True, exist_ok=True)
    return {"SUB6_LAUNCH_ROOT": str(d)}


def step_env():
    import importlib.metadata as md
    info = {"python": sys.version.split()[0], "platform": platform.platform(), "logical_cpus": os.cpu_count(),
            "pyaedt": md.version("pyaedt"), "numpy": md.version("numpy")}
    root = os.environ.get("ANSYSEM_ROOT251")
    cands = [Path(root)] if root else []
    cands += [Path(r"C:\Program Files\ANSYS Inc\v251\AnsysEM"), Path(r"C:\Program Files\AnsysEM\v251\Win64")]
    info["aedt_v251_found"] = any((c / "ansysedt.exe").exists() for c in cands)
    (WORK / "lab_B_machine.json").write_text(json.dumps(info, indent=2))
    say(f"ENV {info}")
    if not info["aedt_v251_found"]:
        say("WARNING: AEDT 2025 R1 (v251) not found at the default locations; set ANSYSEM_ROOT251 or report the installed version (do NOT proceed silently)")
    return info


def step_verify():
    out = WORK / "verify"; out.mkdir(parents=True, exist_ok=True)
    rc = run([PY, "-X", "utf8", HERE / "cell_param_pyaedt.py", "--all", "--diode", "bar50", "--floor", "pec", "--floor-z", "-33.21", "--no-mesh-ops",
              "--cores", os.environ.get("LAB_VERIFY_CORES", "16"), "--row-timeout-min", "30", "--project", WORK / "verify_project" / "CellParam_CAL_NOMESH.aedt",
              "--output", out], env_extra=launch_env())
    say(f"VERIFY run exit {rc}")
    if not (out / "verify.json").exists():
        say("VERIFY produced no verify.json -> STOP"); return False
    rc = run([PY, "-X", "utf8", HERE / "check_verify_vs_reference.py", "--lab", out / "verify.json", "--out", out / "verify_vs_reference.json"])
    say("VERIFY vs reference: " + ("PASS (f within 0.5 %)" if rc == 0 else "FAIL -> STOP, do not run the queue; report the table"))
    return rc == 0


def step_build():
    ok = True
    for j in JOBS:
        if (j["build"] / "build_summary.json").exists() and json.loads((j["build"] / "build_summary.json").read_text()).get("status") == "BUILT":
            say(f"BUILD {j['name']}: already BUILT"); continue
        j["build"].mkdir(parents=True, exist_ok=True); j["project"].parent.mkdir(parents=True, exist_ok=True)
        cmd = [PY, "-X", "utf8", HERE / "cell_param_pyaedt.py", "--build", "--diode", "bar90", "--floor", "impedance377", "--air-below", "30", "--no-mesh-ops",
               "--variant", j["variant"], "--project", j["project"], "--output", j["build"]]
        rc = run(cmd, env_extra=launch_env())      # one build at a time (PyAEDT validation uses a shared temp file)
        bs = j["build"] / "build_summary.json"
        status = json.loads(bs.read_text()).get("status") if bs.exists() else "NO_SUMMARY"
        say(f"BUILD {j['name']}: exit {rc}, status {status}")
        ok = ok and status == "BUILT"
    return ok


def rp_cmd(job, csvs, *extra, out=None):
    return [PY, "-X", "utf8", HERE / "run_parallel.py", "--project", job["project"], "--build-output", job["build"], "--verify-gate", GATE,
            "--csv", *csvs, "--out", out or job["out"], "--diode", "bar90", "--row-timeout-min", "30", *extra]


def step_benchmark():
    """8 real rows (first 8 rows of SCREEN_2): 8 procs x 6 cores vs 4 x 13; prints the faster layout."""
    job = JOBS[0]
    rows = list(csv.DictReader((HERE / "SCREEN_2.csv").open(newline="", encoding="utf-8")))[:8]
    bench = WORK / "BENCH8.csv"
    with bench.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    outB = WORK / "run_full_bench_4x13"
    run(rp_cmd(job, [bench], "--layout", "8x6", "--stagger", "30"), env_extra=launch_env())
    A = [r for r in json.loads((job["out"] / "parallel_report.json").read_text()) if r["file"] == "BENCH8.csv"][-1]
    run(rp_cmd(job, [bench], "--layout", "4x13", "--stagger", "30", out=outB), env_extra=launch_env())
    B = [r for r in json.loads((outB / "parallel_report.json").read_text()) if r["file"] == "BENCH8.csv"][-1]
    win = "6" if A["wall_seconds"] <= B["wall_seconds"] else "13"
    (WORK / "layout_benchmark.json").write_text(json.dumps({"8x6": A, "4x13": B, "faster_cores_each": win}, indent=2))
    say(f"BENCHMARK 8x6 {A['wall_seconds']} s (free CPUs {A['free_cpus_at_start']}) vs 4x13 {B['wall_seconds']} s (free CPUs {B['free_cpus_at_start']}) -> --cores-each {win}")
    return win


def step_queue(cores_each, max_proc):
    if not GATE.exists() or json.loads(GATE.read_text()).get("status") != "PASS":
        say("VERIFY gate missing or not PASS -> STOP"); return False
    chk = WORK / "verify" / "verify_vs_reference.json"
    if not chk.exists() or json.loads(chk.read_text()).get("status") != "PASS":
        say("VERIFY vs reference not PASS -> STOP"); return False
    flags = ["--cores-each", str(cores_each), "--fixed-cores", "--max-proc", str(max_proc), "--stagger", "40"]
    for j in JOBS:     # SCREEN_2 first (stub 0.3/1/2 for the stub-spread result), then the short islands
        say(f"QUEUE_B {j['name']}: {', '.join(j['csv'])}")
        run(rp_cmd(j, j["csv"], *flags), env_extra=launch_env())   # never stops for failed rows; FAILED lists are written
    return True


def step_collect():
    sys.path.insert(0, str(HERE))
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows_dir = RESULTS / "rows"; rows_dir.mkdir(exist_ok=True)
    n = 0
    for j in JOBS:
        if not j["out"].exists():
            continue
        for js in j["out"].glob("*.json"):
            if js.name.startswith("FAILED_ROWS_"):
                shutil.copy(js, RESULTS / js.name)
            elif js.name == "parallel_report.json":
                shutil.copy(js, RESULTS / f"parallel_report_{j['name']}.json")
            else:
                try:
                    if "modes" in json.loads(js.read_text()):
                        shutil.copy(js, rows_dir / js.name); n += 1
                except Exception:
                    pass
        bs = j["build"] / "build_summary.json"
        if bs.exists():
            shutil.copy(bs, RESULTS / f"build_summary_{j['name']}.json")
    v = RESULTS / "verify"; v.mkdir(exist_ok=True)
    for name in ("verify.json", "VERIFY_modes.csv", "verify_vs_reference.json", "environment.json"):
        if (WORK / "verify" / name).exists():
            shutil.copy(WORK / "verify" / name, v / name)
    for name in ("lab_B_machine.json", "layout_benchmark.json", "run_queue_b.log"):
        if (WORK / name).exists():
            shutil.copy(WORK / name, RESULTS / name)
    try:
        import failed_rows_report as fr
        o = fr.write_report([j["out"] for j in JOBS], RESULTS / "failed_report")
        say(f"FAILED rows: {o['rows_failed']} of {o['rows_attempted']} ({o['failure_rate_pct']} %); flags: {o['flags'] or 'none'}")
    except Exception as exc:
        say(f"failed-row report error: {exc!r}")
    say(f"COLLECT: {n} row JSONs -> {RESULTS}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", default="all", choices=["all", "env", "verify", "build", "benchmark", "queue", "collect"])
    ap.add_argument("--cores-each", type=int, default=8, help="cores per solver process (size to the CPU free on THIS machine)")
    ap.add_argument("--max-proc", type=int, default=16)
    a = ap.parse_args()
    if a.step in ("all", "env"):
        step_env()
    if a.step in ("all", "verify"):
        if not step_verify():
            sys.exit(1)
    if a.step in ("all", "build"):
        if not step_build():
            say("a build failed -> STOP"); sys.exit(1)
    if a.step == "benchmark":
        step_benchmark()
    if a.step in ("all", "queue"):
        if not step_queue(a.cores_each, a.max_proc):
            sys.exit(1)
    if a.step in ("all", "collect"):
        step_collect()


if __name__ == "__main__":
    main()
