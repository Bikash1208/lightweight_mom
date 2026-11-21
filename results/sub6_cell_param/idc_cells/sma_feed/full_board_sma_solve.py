"""Solve an existing SMA-feed board project (built natively on THIS machine) to convergence and record S11 at 5.68 GHz.
Saves the project and captures the AEDT messages and profile BEFORE any pass/fail check, so a failed or non-converged solve never loses its result.
Usage: python full_board_sma_solve.py --project X.aedt --design SMA_B_ALL_DIODE --out run_dir [--zr 24.0 --zx -234.0] [--aedt-version 2026.1] [--cores 16]
"""
import argparse
import json
import os
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import psutil


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", type=Path, required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--zr", default=None, help="set the Zr design variable (ohm/sq) before solving")
    ap.add_argument("--zx", default=None)
    ap.add_argument("--aedt-version", default="2025.1")
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--max-passes", type=int, default=12)
    ap.add_argument("--min-passes", type=int, default=6)
    a = ap.parse_args()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    temp = out / "aedt_temp"; temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp))
    root_env = "ANSYSEM_ROOT" + a.aedt_version.replace(".", "")
    if root_env not in os.environ:
        os.environ[root_env] = str(Path("C:/Program Files/ANSYS Inc") / ("v" + a.aedt_version.replace(".", "")) / "AnsysEM")
    from ansys.aedt.core import Hfss
    sp = out / "status.json"
    st = {"project": str(a.project), "design": a.design, "phase": "START", "start_utc": now(), "controller_pid": os.getpid(), "aedt_pid": None}
    sp.write_text(json.dumps(st, indent=1))
    done = threading.Event()

    def beat():
        while not done.wait(60):
            st["heartbeat_utc"] = now()
            try:
                st["free_memory_gib"] = round(psutil.virtual_memory().available / 2**30, 3)
            except psutil.Error:
                pass
            sp.write_text(json.dumps(st, indent=1, default=str))

    threading.Thread(target=beat, daemon=True).start()
    h = None
    try:
        h = Hfss(project=str(a.project.resolve()), design=a.design, solution_type="Modal", version=a.aedt_version,
                 non_graphical=True, new_desktop=True, close_on_exit=False)
        st["aedt_pid"] = h.desktop_class.aedt_process_id
        st["aedt_version"] = str(h.desktop_class.aedt_version_id)
        if st["aedt_version"] != a.aedt_version:
            raise RuntimeError("WRONG_AEDT_VERSION")
        h.desktop_class.temp_directory = temp
        if a.zr is not None:
            h["Zr"] = f"{float(a.zr)}ohm"; h["Zx"] = f"{float(a.zx)}ohm"
        st["Zr"] = str(h["Zr"]) if "Zr" in h.variable_manager.variables else None
        st["Zx"] = str(h["Zx"]) if "Zx" in h.variable_manager.variables else None
        setup = h.get_setup("Setup1")
        setup.props["MaximumPasses"] = a.max_passes; setup.props["MinimumPasses"] = a.min_passes
        setup.props["MaxDeltaS"] = 0.02; setup.props["BasisOrder"] = 1; setup.props["DrivenSolverType"] = "Direct"
        if not setup.update():
            raise RuntimeError("SETUP_UPDATE")
        st["setup_readback"] = {k: str(setup.props.get(k)) for k in ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS", "BasisOrder", "DrivenSolverType")}
        h.save_project()
        st["phase"] = "SOLVING"; st["solve_start_utc"] = now()
        sp.write_text(json.dumps(st, indent=1, default=str))
        t0 = time.monotonic()
        ok = h.analyze_setup("Setup1", cores=a.cores, use_auto_settings=False, blocking=True)
        st["wall_s"] = round(time.monotonic() - t0, 1)
        h.save_project()                                   # save BEFORE any check
        try:
            st["aedt_messages"] = [str(x)[:300] for x in list(h.odesktop.GetMessages("", "", 0))[-40:]]
        except BaseException as e:
            st["aedt_messages_error"] = repr(e)
        try:
            st["profile_path"] = str(Path(h.export_profile("Setup1", output_file=str(out / "solve.prof"))))
            text = Path(st["profile_path"]).read_text(errors="replace")
            st["converged"] = text.rfind("Adaptive Passes converged") > text.rfind("Adaptive Passes did not converge")
            st["normal_completion"] = "Status: Normal Completion" in text
        except BaseException as e:
            st["profile_error"] = repr(e)
        try:
            sd = h.post.get_solution_data(expressions=["S(Port1,Port1)"], setup_sweep_name="Setup1 : LastAdaptive")
            fr, re = sd.get_expression_data(formula="real"); _, im = sd.get_expression_data(formula="imag")
            st["S11"] = {"f_ghz": [float(x) for x in fr], "re": [float(x) for x in re], "im": [float(x) for x in im]}
        except BaseException as e:
            st["S11_error"] = repr(e)
        st["analyze_ok"] = bool(ok)
        st["phase"] = "COMPLETE" if ok else "SOLVE_FAILED"
        sp.write_text(json.dumps(st, indent=1, default=str))
        print(st["phase"], flush=True)
    except BaseException:
        st["phase"] = "FAILED"; st["traceback"] = traceback.format_exc()
        sp.write_text(json.dumps(st, indent=1, default=str))
        print(st["traceback"], file=sys.stderr, flush=True)
        raise
    finally:
        done.set()
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
