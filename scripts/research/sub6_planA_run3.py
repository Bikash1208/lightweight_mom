"""One-at-a-time attached PyAEDT Stage 1 anchor solves with a private temp folder."""
import json
import math
import os
import re
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
SOURCE = ROOT / "stage1_work"
WORK = ROOT / "stage1_run3_work"
TEMP = ROOT / "aedt_temp"
OUT = ROOT / "results" / "sub6_planA" / "stage1" / "run3"
STEP = Path(r"C:\lwmom_labB_anchor\hfss_exports\sub6_board\BOARD_N240_RF_full_win.step")
PY = Path(sys.executable)
ORDER = [(v, s) for v in ("V", "T2", "T5", "T1", "T8") for s in ("OFF", "ON")]
STATUS = WORK / "controller_status.json"
LOG = WORK / "controller.log"


def utc():
    return datetime.now(timezone.utc).isoformat()


def write(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, path)


def check_build(run, variant, state, report):
    assert report["run_id"] == run and report["state"] == state
    assert report["status"] == "BUILT_UNSOLVED" and report["geometry_gate"] == "PASS"
    assert report["aedt_version"] == "2025.2"
    assert report["counts_before"] == {"metal_bodies": 1000, "diode_sheets": 240,
                                        "via_bodies": 480, "substrates": 2, "total": 1244}
    thickness = {"V": .8, "T8": .8, "T5": .508, "T2": .254, "T1": .127}[variant]
    assert report["target_bias_thickness_mm"] == thickness
    assert abs(report["bias_after"][2] - (-4.010 + .8 - thickness)) < 1e-5
    assert report["material_readback"]["SUBSTRATE"] == "RT5880_STAGE1"
    assert report["material_readback"]["SUBSTRATE_FR4"] == ("vacuum" if variant == "V" else "RT5880_STAGE1")
    assert len(report["impedance_readback"]["Objects"]) == 240
    rs, xs = ((3.885, -58.70) if state == "OFF" else (.28, 5.11))
    assert float(report["impedance_readback"]["Resistance"]) == rs
    assert float(report["impedance_readback"]["Reactance"]) == xs
    assert report["setup_readback"] == {"Frequency": "5.8GHz", "MaximumPasses": 4, "MaxDeltaS": .02}


def private_temp(path):
    path.mkdir(parents=True, exist_ok=True)
    probe = path / "write_probe.tmp"
    probe.write_text("ok")
    assert probe.read_text() == "ok"
    probe.unlink()  # exact file created above
    os.environ["TEMP"] = str(path)
    os.environ["TMP"] = str(path)
    os.environ["AEDT_PRIVATE_TEMP"] = str(path)


def rebuild_from_step(run, variant, state, folder, temp):
    assert STEP.is_file(), "STEP_MISSING"
    pre = WORK / "rebuilt_preflight"
    pre.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.update(TEMP=str(temp), TMP=str(temp), AEDT_PRIVATE_TEMP=str(temp))
    cmds = [
        [PY, ROOT / "scripts/research/sub6_planA_anchor_preflight.py", STEP, pre],
        [PY, ROOT / "scripts/research/sub6_planA_anchor_build.py", "--base", pre / "anchor_import_preflight.aedt",
         "--out", folder, "--variant", variant, "--state", state],
    ]
    for i, cmd in enumerate(cmds):
        with (folder / f"rebuild_{i}.log").open("w", encoding="utf-8") as f:
            p = subprocess.run([str(x) for x in cmd], cwd=ROOT, env=env, stdout=f, stderr=subprocess.STDOUT,
                               creationflags=subprocess.CREATE_NO_WINDOW)
        if p.returncode:
            raise RuntimeError(f"REBUILD_FROM_STEP_FAILED step={i} exit={p.returncode}")
    report = json.loads((folder / f"{run}_build.json").read_text())
    check_build(run, variant, state, report)
    return folder / f"{run}.aedt", report


def profile_metrics(path):
    data = path.read_text(errors="replace")
    blocks = re.split(r"Name='Adaptive Pass (\d+)'", data)
    passes = []
    for i in range(1, len(blocks), 2):
        n, part = int(blocks[i]), blocks[i + 1]
        tet = re.search(r"ProfileItem\('Matrix Assembly'.*?Tetrahedra\\', (\d+)", part)
        ds = re.findall(r"Max Mag\. Delta S\\', ([\d.]+)", part)
        passes.append({"pass": n, "tetrahedra": int(tet.group(1)) if tet else None,
                       "delta_s": float(ds[-1]) if ds else None})
    assert passes and passes[-1]["tetrahedra"] is not None, "PROFILE_METRICS_MISSING"
    return {"adaptive_passes": len(passes), "final_tetrahedra": passes[-1]["tetrahedra"],
            "final_delta_s": passes[-1]["delta_s"], "pass_metrics": passes}


def extract(h, run, out, report):
    from sub6_planA_restart_extract import ffd_data, point
    import numpy as np

    out.mkdir(parents=True, exist_ok=True)
    report["profile_file"] = str(h.export_profile("Setup1", output_file=str(out / "solve.prof")))
    report["convergence_file"] = str(h.export_convergence("Setup1", output_file=str(out / "solve.conv")))
    report.update(profile_metrics(Path(report["profile_file"])))
    far = out / "farfield"
    far.mkdir(exist_ok=True)
    assert h.export_antenna_metadata(frequencies=[5.8e9], setup="Setup1 : LastAdaptive",
                                     sphere="Global", output_dir=str(far), export_element_pattern=True,
                                     export_objects=False, export_touchstone=True, export_power=True), "FFD_EXPORT_FAILED"
    ffd = list(far.rglob("*.ffd"))
    assert len(ffd) == 1, f"FFD_COUNT_{len(ffd)}"
    report["ffd_file"] = str(ffd[0])
    report["s11_db"] = point(h, "dB(S(Port1,Port1))")
    for theta, phi in ((25, 0), (25, 180), (24, 180)):
        report[f"th{theta}_ph{phi}"] = {
            "gain_rhcp_db": point(h, "dB(GainRHCP)", theta, phi),
            "gain_lhcp_db": point(h, "dB(GainLHCP)", theta, phi)}
    th, ph, et, ep = ffd_data(ffd[0])
    assert th[0] == 0 and th[-1] == 180 and ph[0] == -180 and ph[-1] == 180
    intensity = np.abs(et)**2 + np.abs(ep)**2
    weights = np.sin(np.deg2rad(th))[:, None]
    total = float(np.sum(intensity[:, :-1] * weights))
    below = float(np.sum(intensity[th > 90, :-1] * weights[th > 90]))
    assert total > 0
    report["radiated_power_fraction_theta_gt_90"] = below / total
    plus = np.abs(et + 1j * ep)**2
    minus = np.abs(et - 1j * ep)**2
    i24 = int(np.argmin(abs(th - 24)))
    i0 = int(np.argmin(abs(ph)))
    rh_db = point(h, "dB(GainRHCP)", 24, 0)
    lh_db = point(h, "dB(GainLHCP)", 24, 0)
    ratio = 10 * math.log10(max(plus[i24, i0], 1e-300) / max(minus[i24, i0], 1e-300))
    handed = "plus" if abs(ratio - (rh_db-lh_db)) < abs(-ratio - (rh_db-lh_db)) else "minus"
    pattern = plus if handed == "plus" else minus
    idx = np.unravel_index(int(np.argmax(pattern[:, :-1])), pattern[:, :-1].shape)
    peak_th, peak_ph = float(th[idx[0]]), float(ph[idx[1]])
    report["peak_rhcp"] = {"theta_deg": peak_th, "phi_deg": peak_ph,
                           "gain_db": point(h, "dB(GainRHCP)", peak_th, peak_ph),
                           "ffd_handedness": handed}


def push_pair(index, run_names):
    def git(*args):
        p = subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True)
        if p.returncode:
            raise RuntimeError(f"GIT_{args[0]}_FAILED: {p.stderr.strip()}")
        return p.stdout.strip()
    if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=ROOT).returncode:
        raise RuntimeError("UNRELATED_STAGED_CHANGES")
    paths = [str((OUT / f"{run}.json").relative_to(ROOT)) for run in run_names]
    paths += [str((OUT / run / "farfield").relative_to(ROOT)) for run in run_names]
    paths += [str((OUT / run / "solve.prof").relative_to(ROOT)) for run in run_names]
    paths += [str((OUT / run / "solve.conv").relative_to(ROOT)) for run in run_names]
    if index == 1:
        paths += ["scripts/research/sub6_planA_run3.py", "scripts/research/sub6_planA_run3_watchdog.py",
                  "scripts/research/sub6_planA_anchor_preflight.py", "scripts/research/sub6_planA_anchor_build.py"]
    git("add", "--", *paths)
    git("commit", "-m", f"Stage 1 Run 3 pair {index}")
    git("push", "gh", "compute/sub6-cell-param-lab-b")
    return git("rev-parse", "HEAD")


def main():
    WORK.mkdir(exist_ok=True)
    TEMP.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    state = {"status": "STARTING", "controller_pid": os.getpid(), "started_utc": utc(),
             "run_order": [f"ANCHOR_{v}_{s}" for v, s in ORDER], "rows": {}}
    stop = threading.Event()

    def heartbeat():
        while not stop.is_set():
            state["heartbeat_utc"] = utc()
            write(STATUS, state)
            stop.wait(60)

    hb = threading.Thread(target=heartbeat, daemon=True)
    hb.start()
    parent = psutil.Process()
    with (WORK / "watchdog.log").open("w", encoding="utf-8") as wlog:
        subprocess.Popen([str(PY), str(Path(__file__).with_name("sub6_planA_run3_watchdog.py")),
                          str(parent.pid), str(parent.create_time()), str(STATUS), str(LOG)],
                         cwd=ROOT, stdout=wlog, stderr=subprocess.STDOUT,
                         creationflags=subprocess.DETACHED_PROCESS | subprocess.CREATE_NO_WINDOW)
    h = None
    try:
        from ansys.aedt.core import Hfss
        assert STEP.is_file()
        assert len(psutil.Process().cpu_affinity()) == 48, "EXPECTED_48_LOGICAL_CPUS"
        for k, (variant, state_name) in enumerate(ORDER, 1):
            run = f"ANCHOR_{variant}_{state_name}"
            folder = WORK / run
            folder.mkdir(exist_ok=True)
            run_temp = TEMP / run
            private_temp(run_temp)
            src = SOURCE / f"{run}.aedt"
            build = json.loads((SOURCE / f"{run}_build.json").read_text())
            check_build(run, variant, state_name, build)
            project = folder / f"{run}.aedt"
            assert src.is_file() and src.stat().st_size > 10_000_000
            if not project.exists():
                shutil.copy2(src, project)
            row = {"run_id": run, "status": "OPENING", "started_utc": utc(),
                   "project": str(project), "build_report": str(SOURCE / f"{run}_build.json"),
                   "requested_cores": 48}
            state["rows"][run] = row
            state["active_run"] = run
            state["status"] = "OPENING"
            write(STATUS, state)
            try:
                h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
                         version="2025.2", non_graphical=True, new_desktop=True,
                         close_on_exit=False)
                if len(h.modeler.object_names) != 1244:
                    raise RuntimeError("PROJECT_READBACK_COUNTS_OR_FILE_FORMAT")
            except Exception as exc:
                if h is not None:
                    h.release_desktop(close_projects=True, close_desktop=True)
                    h = None
                if "format" not in str(exc).lower() and "readback" not in str(exc).lower():
                    raise
                project, build = rebuild_from_step(run, variant, state_name, folder, run_temp)
                h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
                         version="2025.2", non_graphical=True, new_desktop=True,
                         close_on_exit=False)
                row["rebuilt_from_step"] = True
                assert len(h.modeler.object_names) == 1244, "REBUILT_READBACK_COUNTS"
            row["aedt_version"] = str(h.desktop_class.aedt_version_id)
            state["own_aedt_pid"] = int(h.desktop_class.aedt_process_id)
            state["own_aedt_created"] = psutil.Process(state["own_aedt_pid"]).create_time()
            assert row["aedt_version"] == "2025.2"
            h.desktop_class.temp_directory = run_temp
            row["aedt_temp_readback"] = str(h.desktop_class.temp_directory)
            assert Path(row["aedt_temp_readback"]).resolve() == run_temp.resolve(), "PRIVATE_TEMP_READBACK"
            names = list(h.modeler.object_names)
            assert len(names) == 1244
            b = next(x for x in h.boundaries if x.name == "BAR90_02LRH")
            assert len(b.props["Objects"]) == 240
            assert float(b.props["Resistance"]) == float(build["impedance_readback"]["Resistance"])
            assert float(b.props["Reactance"]) == float(build["impedance_readback"]["Reactance"])
            row["setup_names_before"] = [s.name for s in h.setups]
            assert "Setup1" in row["setup_names_before"]
            if "Auto1" in row["setup_names_before"]:
                assert h.delete_setup("Auto1"), "AUTO1_DELETE_FAILED"
            row["setup_names_after"] = [s.name for s in h.setups]
            assert row["setup_names_after"] == ["Setup1"], row["setup_names_after"]
            setup = h.setups[0]
            assert str(setup.props["Frequency"]) == "5.8GHz"
            assert int(setup.props["MaximumPasses"]) == 4
            validation, valid = h.validate_full_design(output_dir=str(folder), ports=1)
            row["validation"] = validation
            assert valid, "FULL_DESIGN_VALIDATION_FAILED"
            h.save_project()
            row["status"] = "SOLVING"
            row["solve_start_utc"] = utc()
            state["status"] = "SOLVING"
            write(STATUS, state)
            peak = [0]
            done = threading.Event()

            def ram_monitor():
                while not done.wait(5):
                    try:
                        proc = psutil.Process(state["own_aedt_pid"])
                        tree = [proc] + proc.children(recursive=True)
                        peak[0] = max(peak[0], sum(p.memory_info().rss for p in tree if p.is_running()))
                    except psutil.Error:
                        pass

            monitor = threading.Thread(target=ram_monitor, daemon=True)
            monitor.start()
            start = time.monotonic()
            ok = False
            for attempt in (1, 2):
                row["solve_attempt"] = attempt
                write(STATUS, state)
                ok = bool(h.analyze_setup("Setup1", cores=48, use_auto_settings=False,
                                          revert_to_initial_mesh=(attempt == 2), blocking=True))
                if ok:
                    break
            done.set()
            monitor.join(timeout=10)
            row["solve_seconds"] = round(time.monotonic() - start, 2)
            row["peak_own_process_tree_ram_gb"] = round(peak[0] / 2**30, 3)
            row["solve_end_utc"] = utc()
            assert ok, "SOLVE_FAILED_AFTER_RETRY"
            h.save_project()
            row["status"] = "SOLVED_EXTRACTING"
            state["status"] = "SOLVED_EXTRACTING"
            write(STATUS, state)
            extract(h, run, OUT / run, row)
            row["status"] = "COMPLETE"
            row["ended_utc"] = utc()
            write(OUT / f"{run}.json", row)
            h.release_desktop(close_projects=True, close_desktop=True)
            h = None
            state["own_aedt_pid"] = None
            state["own_aedt_created"] = None
            state["status"] = "ROW_COMPLETE"
            write(STATUS, state)
            if k in (2, 4):
                pair = k // 2
                state[f"pair_{pair}_push_head"] = push_pair(pair, [f"ANCHOR_{variant}_{s}" for s in ("OFF", "ON")])
                state[f"pair_{pair}_pushed_utc"] = utc()
                state["status"] = f"PAIR_{pair}_PUSHED"
                write(STATUS, state)
        state["final_push_head"] = push_pair(5, [f"ANCHOR_{v}_{s}" for v in ("T5", "T1", "T8") for s in ("OFF", "ON")])
        state["status"] = "COMPLETE"
    except BaseException:
        state["status"] = "ABORTED"
        state["traceback"] = traceback.format_exc()
        with LOG.open("a", encoding="utf-8") as f:
            f.write(state["traceback"])
        raise
    finally:
        if h is not None:
            try:
                h.release_desktop(close_projects=True, close_desktop=True)
            except BaseException:
                with LOG.open("a", encoding="utf-8") as f:
                    f.write(traceback.format_exc())
        pid, created = state.get("own_aedt_pid"), state.get("own_aedt_created")
        if pid and created:
            try:
                proc = psutil.Process(pid)
                if abs(proc.create_time() - created) < 1:
                    for child in proc.children(recursive=True):
                        child.terminate()
                    proc.terminate()
            except psutil.Error:
                pass
        state["ended_utc"] = utc()
        write(STATUS, state)
        stop.set()
        hb.join(timeout=2)


if __name__ == "__main__":
    main()
