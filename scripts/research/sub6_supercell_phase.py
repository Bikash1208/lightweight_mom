"""Rotated 8-cell eigenmode solve; run only after the unit-cell PML gate passes."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import threading
import time
import traceback

from sub6_rotated_supercell_hfss import (BIAS, CELLS, CU, P, RESULTS, ROOT,
                                         RT, WORK, launch, modes, need, now, write)


def build(h, cell, phase, min_freq_ghz=5.3, num_modes=20):
    spec = CELLS[cell]
    m = h.modeler
    m.model_units = "mm"
    mat = h.materials.add_material("RT5880_SUPERCELL")
    mat.permittivity = 2.2
    mat.dielectric_loss_tangent = .0009
    h["$ph_s"] = f"{phase}deg"
    h["$ph_y"] = "0deg"
    for key, val in dict(period=P, cu_t=CU, rt_t=RT, bias_t=BIAS,
                         air_top=40.035, air_bottom=-43.499).items():
        h[key] = f"{val:.12g}mm"

    def box(name, pos, size, material="pec"):
        return need(m.create_box(pos, size, name=name, material=material), f"CREATE_{name}")

    cells = []
    for n in range(8):
        disc = need(m.create_cylinder("Z", ["0mm", "0mm", "0mm"],
                                      f"{spec['r']}mm", "cu_t",
                                      name=f"DISC_{n}", material="pec"), f"DISC_{n}")
        cuts = []

        def cut(name, x0, y0, dx, dy):
            cuts.append(box(f"{name}_{n}", [f"{x0:.9f}mm", f"{y0:.9f}mm", "-cu_t"],
                            [f"{dx:.9f}mm", f"{dy:.9f}mm", "3*cu_t"], "vacuum"))

        if spec["kind"] == "c4x":
            cut("SLOT_X", -4.36, -1.455, 8.72, 2.91)
            if spec["on"]:
                cut("SLOT_Y_N", -1.455, 1.955, 2.91, 2.405)
                cut("SLOT_Y_S", -1.455, -4.36, 2.91, 2.405)
            else:
                cut("SLOT_Y", -1.455, -4.36, 2.91, 8.72)
        elif spec["kind"] == "st":
            cut("SLOT", -spec["slot_w"] / 2, -spec["slot_l"] / 2,
                spec["slot_w"], spec["slot_l"])
        else:
            w, depth, bridge = spec["w"], spec["depth"], spec["bridge"]
            for i in range(5):
                cut(f"LEG_{i}", -2.4 + i * 1.2 - w / 2,
                    -depth / 2 - w / 2, w, depth + w)
            for i in range(4):
                yc = depth / 2 if i % 2 == 0 else -depth / 2
                cut(f"BAR_{i}", -2.4 + i * 1.2 - w / 2, yc - w / 2, 1.2 + w, w)
            rimx = math.sqrt(spec["r"] ** 2 - (depth / 2 + w / 2) ** 2)
            inner_rimx = math.sqrt(spec["r"] ** 2 - (depth / 2 - w / 2) ** 2)
            # The open cut must cross the rim over its full 0.4-mm width.
            # rimx leaves a 0.136-mm bridge at the inner y edge for D3p6.
            end_abs = spec["r"] + 0.1 if bridge is None else rimx - bridge
            need(end_abs > 2.2, "END_RUN_LENGTH_NONPOSITIVE")
            cut("END_LEFT", -end_abs, -depth / 2 - w / 2, end_abs - 2.2, w)
            cut("END_RIGHT", 2.2, depth / 2 - w / 2, end_abs - 2.2, w)
        need(m.subtract(disc.name, [c.name for c in cuts], keep_originals=False),
             f"SLOT_SUBTRACT_{n}")
        if spec["kind"] == "mn" and bridge is None:
            pre = [{"name": name, "material": m[name].material_name,
                    "volume_mm3": float(m[name].volume)}
                   for name in m.object_names if name.startswith(f"DISC_{n}")]
            separated = need(m.separate_bodies(disc.name), f"SEPARATE_TWIN_{n}")
            names = [part.name for part in separated]
            parts = [m[name] for name in names]
            need(len(parts) == 2, f"TWIN_PART_COUNT_{n}: {names}")
            need(all(part.material_name.lower() == "pec" for part in parts),
                 f"TWIN_MATERIAL_{n}")
            areas = [float(part.volume) / CU for part in parts]
            need(all(math.isclose(a, spec["area"] / 2, rel_tol=.02) for a in areas),
                 f"TWIN_HALF_AREA_{n}: {areas}")
            need(math.isclose(sum(areas), spec["area"], rel_tol=.02),
                 f"TWIN_TOTAL_AREA_{n}: {sum(areas)}")
            geometry_audit = {"objects_after_subtract": pre,
                              "end_run_original_mm": rimx,
                              "rim_at_inner_y_mm": inner_rimx,
                              "end_run_used_mm": end_abs,
                              "separated_names": names}
        else:
            parts = [disc]
            areas = [float(disc.volume) / CU]
            need(math.isclose(areas[0], spec["area"], rel_tol=.02),
                 f"CELL_COPPER_AREA_{n}: {areas[0]} expected {spec['area']}")
            if spec["kind"] == "c4x":
                need(abs(areas[0] - spec["area"]) < .03,
                     f"C4X_CELL_AREA_{n}: {areas[0]} expected {spec['area']}")
            geometry_audit = None
        for part in parts:
            need(part.rotate("Z", n * 22.5), f"ROTATE_{n}_{part.name}")
            need(part.move([(n - 3.5) * P, 0, 0]), f"MOVE_{n}_{part.name}")
        cells.append({"names": [part.name for part in parts],
                      "rotation_deg": n * 22.5,
                      "centre_x_mm": (n - 3.5) * P,
                      "areas_mm2": areas, "area_mm2": sum(areas),
                      "geometry_audit": geometry_audit})

    box("GROUND", ["-4*period", "-period/2", "-rt_t-cu_t"],
        ["8*period", "period", "cu_t"])
    substrate = box("SUBSTRATE", ["-4*period", "-period/2", "-rt_t"],
                    ["8*period", "period", "rt_t"], "RT5880_SUPERCELL")
    box("BIAS_RT5880", ["-4*period", "-period/2", "-rt_t-cu_t-bias_t"],
        ["8*period", "period", "bias_t"], "RT5880_SUPERCELL")
    box("AIR", ["-4*period", "-period/2", "air_bottom"],
        ["8*period", "period", "air_top-air_bottom"], "vacuum")
    face = lambda xyz: need(m.get_faceid_from_position(xyz, assignment="AIR"), "AIR_FACE")
    zm = (40.035 - 43.499) / 2
    need(h.assign_lattice_pair([face([-4 * P, 0, zm]), face([4 * P, 0, zm])],
                               phase_delay="InputPhaseDelay", phase_delay_param1="$ph_s",
                               name="Lattice_X"), "LATTICE_X")
    need(h.assign_lattice_pair([face([0, -P / 2, zm]), face([0, P / 2, zm])],
                               phase_delay="InputPhaseDelay", phase_delay_param1="$ph_y",
                               name="Lattice_Y"), "LATTICE_Y")
    h.oboundary.CreatePML(["NAME:PMLCreationSettings", "UserDrawnGroup:=", False,
                           "PMLFaces:=", [face([0, 0, 40.035])],
                           "CreateJoiningObjs:=", False, "Thickness:=", "10mm",
                           "RadDist:=", "20mm", "UseFreq:=", True,
                           "MinFreq:=", "5.8GHz"])
    need(h.assign_perfect_e([face([0, 0, -43.499])], name="PEC_BOTTOM"), "PEC_BOTTOM")
    setup = need(h.create_setup("EigenCell"), "SETUP")
    setup.props.update(MinimumFrequency=f"{min_freq_ghz}GHz", NumModes=num_modes, MaxDeltaFreq=.2,
                       MinimumPasses=4, MinimumConvergedPasses=2, MaximumPasses=20,
                       ConvergeOnRealFreq=True, SaveFields=True)
    need(setup.update(), "SETUP_UPDATE")
    validation_file = RESULTS / "validation" / f"SUPER_{cell}_{phase}.log"
    validation_file.parent.mkdir(parents=True, exist_ok=True)
    need(h.validate_simple(validation_file) == 1, "DESIGN_VALIDATION")
    boundaries = {b.name: dict(b.props) for b in h.boundaries}
    need(all(k in boundaries for k in ("Lattice_X", "Lattice_Y", "PEC_BOTTOM")),
         "BOUNDARY_READBACK")
    need(any("pml" in n.lower() for n in m.object_names), "PML_READBACK")
    need(substrate.material_name.lower() == "rt5880_supercell", "MATERIAL_READBACK")
    need(abs(float(substrate.bounding_box[5]) - float(substrate.bounding_box[2]) - RT) < 1e-6,
         "THICKNESS_READBACK")
    need(h.save_project(), "SAVE_PROJECT")
    return {"cell": cell, "phase_deg": phase, "cells": cells,
            "boundaries": boundaries, "objects": m.object_names,
            "setup": dict(setup.props), "aedt_version": h.desktop_class.aedt_version_id}


def one(cell, phase, cores, tag="", min_freq_ghz=5.3, num_modes=20,
        timeout_s=1500, build_only=False):
    label = f"SUPER_{cell}_ph{phase:g}_PML40{('_'+tag) if tag else ''}"
    project = WORK / "projects" / f"{label}.aedt"
    need(not project.exists(), f"REFUSE_EXISTING_PROJECT {project}")
    RESULTS.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    state = {"status": "STARTING", "cell": cell, "phase_deg": phase,
             "solve_timeout_s": timeout_s, "minimum_frequency_ghz": min_freq_ghz,
             "requested_modes": num_modes,
             "controller_pid": os.getpid(), "started_utc": now()}
    status_path = RESULTS / f"{label}_status.json"
    write(status_path, state)
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(60):
            state["heartbeat_utc"] = now()
            if state.get("aedt_pid"):
                try:
                    import psutil
                    state["own_child_pids"] = [p.pid for p in psutil.Process(int(state["aedt_pid"])).children(recursive=True)]
                except psutil.NoSuchProcess:
                    state["own_child_pids"] = []
            write(status_path, state)

    threading.Thread(target=heartbeat, daemon=True).start()
    h = None
    try:
        h = launch(project)
        state.update(status="BUILDING", aedt_pid=h.desktop_class.aedt_process_id)
        write(status_path, state)
        built = build(h, cell, phase, min_freq_ghz, num_modes)
        write(RESULTS / f"{label}_build.json", built)
        if build_only:
            state.update(status="BUILT", build_end_utc=now())
            write(status_path, state)
            return
        state.update(status="SOLVING", solve_start_utc=now())
        write(status_path, state)
        watchdog_log = open(RESULTS / f"{label}_external_watchdog.log", "w", encoding="utf-8")
        external = subprocess.Popen([sys.executable, str(Path(__file__).with_name("sub6_external_watchdog.py")),
                                     str(status_path)], cwd=str(ROOT), stdout=watchdog_log,
                                    stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW)
        state["external_watchdog_pid"] = external.pid
        write(status_path, state)

        def timeout_own_pid():
            import psutil
            state.update(status="TIMEOUT", error=f"{timeout_s}s solve timeout", timeout_utc=now())
            write(status_path, state)
            try:
                p = psutil.Process(int(h.desktop_class.aedt_process_id))
                targets = p.children(recursive=True) + [p]
                state["terminated_own_pids"] = [q.pid for q in targets]
                write(status_path, state)
                for q in targets:
                    if q.is_running():
                        q.terminate()
            except psutil.NoSuchProcess:
                pass

        watchdog = threading.Timer(timeout_s, timeout_own_pid)
        watchdog.start()
        try:
            need(h.analyze_setup("EigenCell", cores=cores, use_auto_settings=True,
                                 blocking=True), "SOLVER_FAILURE")
            need(state["status"] != "TIMEOUT", "SOLVE_TIMEOUT")
        finally:
            watchdog.cancel()
        result = {"status": "COMPLETE", "cell": cell, "phase_deg": phase,
                  "aedt_version": "2025.2", "aedt_pid": h.desktop_class.aedt_process_id,
                  "modes": modes(h), "seconds": time.monotonic() - start,
                  "solve_end_utc": now(), "project": str(project)}
        try:
            result["profile_path"] = h.export_profile("EigenCell", output_file=str(RESULTS / f"{label}_profile.prof"))
            result["convergence_path"] = h.export_convergence("EigenCell", output_file=str(RESULTS / f"{label}_convergence.conv"))
        except Exception as exc:
            result["profile_export_error"] = f"{type(exc).__name__}: {exc}"
        need(h.save_project(), "POSTSOLVE_SAVE_PROJECT")
        result["mode_span_reaches_6p5"] = result["modes"][-1]["f_ghz"] >= 6.5
        if num_modes >= 20:
            need(result["mode_span_reaches_6p5"], "MODE_SPAN_BELOW_6P5")
        write(RESULTS / f"{label}_modes.json", result)
        state.update(status="COMPLETE", seconds=result["seconds"], solve_end_utc=now())
        write(status_path, state)
        print(json.dumps({"cell": cell, "phase_deg": phase,
                          "seconds": result["seconds"], "n_modes": len(result["modes"]),
                          "highest_f_ghz": result["modes"][-1]["f_ghz"]}), flush=True)
    except BaseException as exc:
        state.update(status="FAILED", error=f"{type(exc).__name__}: {exc}",
                     traceback=traceback.format_exc(), failed_utc=now())
        write(status_path, state)
        traceback.print_exc()
        raise
    finally:
        stop.set()
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", choices=("ANCHOR", "B2p2_D3p0", "B2p2_D3p6",
                                       "D3p6_OPEN", "C4X_ON", "C4X_OFF"), required=True)
    ap.add_argument("--phase", type=int, choices=(15, 30, 45), required=True)
    ap.add_argument("--cores", type=int, default=48)
    ap.add_argument("--tag", default="")
    ap.add_argument("--min-freq-ghz", type=float, default=5.3)
    ap.add_argument("--num-modes", type=int, default=20)
    ap.add_argument("--timeout-s", type=int, default=1500)
    ap.add_argument("--build-only", action="store_true")
    args = ap.parse_args()
    one(args.cell, args.phase, args.cores, args.tag, args.min_freq_ghz,
        args.num_modes, args.timeout_s, args.build_only)
