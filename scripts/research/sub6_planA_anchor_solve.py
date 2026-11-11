"""Solve one validated Stage 1 anchor project in a new local AEDT 2025.2 instance."""
import argparse
import json
import os
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import psutil
from ansys.aedt.core import Hfss


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--cores", type=int, required=True)
    a = p.parse_args()
    a.project = a.project.resolve()
    a.out.mkdir(parents=True, exist_ok=True)
    run = a.project.stem
    build = json.loads((a.project.parent / f"{run}_build.json").read_text())
    if build.get("status") != "BUILT_UNSOLVED" or build.get("geometry_gate") != "PASS":
        raise RuntimeError("BUILD_GATE_NOT_PASSED")
    report = a.out / f"{run}.json"
    d = {"run_id": run, "status": "STARTING", "own_python_pid": os.getpid(),
         "started_utc": now(), "cores": a.cores, "build_report": str(a.project.parent / f"{run}_build.json")}
    report.write_text(json.dumps(d, indent=2))
    h = None
    stop = threading.Event()
    peak = {"bytes": 0}
    try:
        h = Hfss(project=str(a.project), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True,
                 close_on_exit=False)
        d["aedt_version"] = str(h.desktop_class.aedt_version_id)
        d["own_aedt_pid"] = h.desktop_class.aedt_process_id
        if d["aedt_version"] != "2025.2":
            raise RuntimeError("AEDT_VERSION")
        validation, valid = h.validate_full_design(output_dir=str(a.out), ports=1)
        d["validation"] = validation
        if not valid:
            raise RuntimeError("AEDT_DESIGN_VALIDATION")
        d["status"] = "SOLVING"
        report.write_text(json.dumps(d, indent=2))

        def monitor():
            while not stop.wait(5):
                try:
                    proc = psutil.Process(d["own_aedt_pid"])
                    own = [proc] + proc.children(recursive=True)
                    rss = sum(q.memory_info().rss for q in own if q.is_running())
                    peak["bytes"] = max(peak["bytes"], rss)
                except psutil.Error:
                    pass

        thread = threading.Thread(target=monitor, daemon=True)
        thread.start()
        for attempt in (1, 2):
            d["attempt"] = attempt
            d["solve_start_utc"] = now()
            report.write_text(json.dumps(d, indent=2))
            t0 = time.monotonic()
            try:
                ok = h.analyze_setup("Setup1", cores=a.cores, use_auto_settings=False,
                                     revert_to_initial_mesh=(attempt == 2), blocking=True)
            except Exception as exc:
                ok = False
                d[f"attempt_{attempt}_error"] = repr(exc)
            d[f"attempt_{attempt}_seconds"] = round(time.monotonic() - t0, 2)
            d[f"attempt_{attempt}_ok"] = bool(ok)
            d[f"attempt_{attempt}_end_utc"] = now()
            report.write_text(json.dumps(d, indent=2))
            if ok:
                break
        stop.set()
        thread.join(timeout=10)
        d["peak_own_process_tree_ram_gb"] = round(peak["bytes"] / 2**30, 3)
        if not ok:
            raise RuntimeError("SOLVE_FAILED_AFTER_RETRY")
        h.save_project()
        d["status"] = "SOLVED_EXPORTING"
        report.write_text(json.dumps(d, indent=2))
        d["profile_file"] = str(h.export_profile("Setup1", output_file=str(a.out / f"{run}.prof")))
        d["convergence_file"] = str(h.export_convergence("Setup1", output_file=str(a.out / f"{run}.conv")))
        d["ffd_export_ok"] = bool(h.export_antenna_metadata(
            frequencies=[5.8e9], setup="Setup1 : LastAdaptive", sphere="Global",
            output_dir=str(a.out / "farfield"), export_element_pattern=True,
            export_objects=False, export_touchstone=True, export_power=True))
        if not d["ffd_export_ok"]:
            raise RuntimeError("FFD_EXPORT_FAILED")
        d["exported_files"] = [str(x) for x in (a.out / "farfield").rglob("*") if x.is_file()]
        d["status"] = "SOLVED_EXPORTED"
    except Exception as exc:
        d["status"] = "FAILED"
        d["error"] = repr(exc)
    finally:
        stop.set()
        d["ended_utc"] = now()
        report.write_text(json.dumps(d, indent=2, default=str))
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
