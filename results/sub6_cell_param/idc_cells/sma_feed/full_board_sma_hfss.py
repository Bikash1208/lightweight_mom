"""Build (and optionally solve) a v3.1 R150 board with the SMA connector feed (sma.step) in HFSS.
Shared script for every machine: pass --aedt-version (2025.1 or 2026.1) and --root-env if the install path differs.
SMA frame z=0 is the flange face; it is placed against the board ground bottom (z -3.21). Pin tip is a parameter (board z, mm).
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

Z_FLANGE = -3.21           # board ground bottom
SMA_PORT_Z = -7.62         # SMA frame
PIN_R = 0.3048
CAP_T = 0.3


def now():
    return datetime.now(timezone.utc).isoformat()


def need(value, label):
    if not value:
        raise RuntimeError(label)
    return value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", type=Path, required=True)
    ap.add_argument("--sma", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--state", required=True)
    ap.add_argument("--pin-tip", type=float, default=8.925, help="pin tip z in the board frame (mm); board top is z 0")
    ap.add_argument("--aedt-version", default="2025.1")
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--max-passes", type=int, default=12)
    ap.add_argument("--min-passes", type=int, default=6)
    ap.add_argument("--sweep", action="store_true", help="add an interpolating 4.5-7 GHz sweep")
    ap.add_argument("--build-only", action="store_true")
    a = ap.parse_args()
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    temp = out / "aedt_temp"
    temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp))
    ver = a.aedt_version
    root_env = "ANSYSEM_ROOT" + ver.replace(".", "")
    if root_env not in os.environ:
        os.environ[root_env] = str(Path("C:/Program Files/ANSYS Inc") / ("v" + ver.replace(".", "")) / "AnsysEM")
    from ansys.aedt.core import Hfss
    status_path = out / "status.json"
    st = {"state": a.state, "phase": "START", "start_utc": now(), "pin_tip_mm": a.pin_tip, "controller_pid": os.getpid(),
          "aedt_pid": None, "aedt_version": None}
    status_path.write_text(json.dumps(st, indent=1))
    done = threading.Event()

    def beat():
        while not done.wait(60):
            st["heartbeat_utc"] = now()
            try:
                if st["aedt_pid"]:
                    p = psutil.Process(st["aedt_pid"])
                    st["own_process_tree_gib"] = round(sum(x.memory_info().rss for x in [p] + p.children(recursive=True) if x.is_running()) / 2**30, 3)
                st["free_memory_gib"] = round(psutil.virtual_memory().available / 2**30, 3)
            except psutil.Error:
                pass
            status_path.write_text(json.dumps(st, indent=1, default=str))

    threading.Thread(target=beat, daemon=True).start()
    h = None
    try:
        project = out / f"BOARD_R150_V31_{a.state}_SMA.aedt"
        need(not project.exists(), "PROJECT_EXISTS")
        h = Hfss(project=str(project), design=f"SMA_{a.state}", solution_type="Modal", version=ver,
                 non_graphical=True, new_desktop=True, close_on_exit=False)
        st["aedt_pid"] = h.desktop_class.aedt_process_id
        st["aedt_version"] = str(h.desktop_class.aedt_version_id)
        need(st["aedt_version"] == ver, "WRONG_AEDT_VERSION")
        h.desktop_class.temp_directory = temp
        m = h.modeler
        m.model_units = "mm"
        rt = h.materials.add_material("RT5880_GATE"); rt.permittivity = 2.2; rt.dielectric_loss_tangent = 0.0009
        ab = h.materials.add_material("ABSORBER_15_10"); ab.permittivity = 1.5; ab.dielectric_loss_tangent = 1 / 1.5
        ab.permeability = 1.5; ab.magnetic_loss_tangent = 1 / 1.5
        pt = h.materials.add_material("PTFE_SMA"); pt.permittivity = 2.1; pt.dielectric_loss_tangent = 0.0002
        st["phase"] = "IMPORT_BOARD"
        need(m.import_3d_cad(str(a.step.resolve())), "BOARD_IMPORT_FAILED")
        board_names = set(m.object_names)
        metals = [n for n in board_names if n.startswith("METAL_PEC")]
        need(len(metals) > 0 and {"SUBSTRATE", "BIAS_RT5880", "AIR"} <= board_names and "PORT" not in board_names, "BOARD_READBACK")
        for n in metals:
            m[n].material_name = "pec"
        m["SUBSTRATE"].material_name = "RT5880_GATE"; m["BIAS_RT5880"].material_name = "RT5880_GATE"; m["AIR"].material_name = "vacuum"
        st["phase"] = "IMPORT_SMA"
        need(m.import_3d_cad(str(a.sma.resolve())), "SMA_IMPORT_FAILED")
        new = [n for n in m.object_names if n not in board_names]
        st["sma_objects"] = new
        find = lambda key: next(n for n in new if key.lower() in n.lower())
        pin, dielec, gnd, circ = find("pin"), find("dielectric"), find("ground"), find("circle")
        need(m.move([pin, dielec, gnd, circ], [0, 0, Z_FLANGE]), "SMA_MOVE_FAILED")
        m[pin].material_name = "pec"; m[gnd].material_name = "pec"; m[dielec].material_name = "PTFE_SMA"
        # pin tip: STEP pin protrudes 7.794 above the flange face -> board z Z_FLANGE + 7.794
        step_tip = Z_FLANGE + 7.794
        st["step_pin_tip_z"] = step_tip
        if a.pin_tip > step_tip + 1e-6:
            cyl = need(m.create_cylinder("Z", [0, 0, step_tip - 0.5], PIN_R, a.pin_tip - step_tip + 0.5, name="PIN_EXT", material="pec"), "PIN_EXT")
            need(m.unite([pin, cyl.name]), "PIN_UNITE")
        elif a.pin_tip < step_tip - 1e-6:
            box = need(m.create_box([-1, -1, a.pin_tip], [2, 2, step_tip - a.pin_tip + 1], name="PIN_TRIM", material="vacuum"), "PIN_TRIM")
            need(m.subtract(pin, box.name, keep_originals=False), "PIN_SUBTRACT")
        # PEC cap behind the wave-port sheet
        zp = Z_FLANGE + SMA_PORT_Z
        cap = need(m.create_cylinder("Z", [0, 0, zp - CAP_T], 2.6926, CAP_T, name="SMA_CAP", material="pec"), "CAP")
        wp = h.wave_port(assignment=circ, integration_line=[[0.635, 0, zp], [2.0574, 0, zp]], modes=1, impedance=50,
                         name="Port1", renormalize=True, deembed=0)
        need(wp, "WAVE_PORT_FAILED")
        need(h.assign_radiation_boundary_to_objects("AIR"), "RADIATION_FAILED")
        outer = need(m.create_cylinder("Z", [0, 0, 0], 186, 30.035, name="ABS_OUT", material="ABSORBER_15_10"), "ABS_OUT")
        inner = need(m.create_cylinder("Z", [0, 0, -1], 156, 32, name="ABS_IN", material="vacuum"), "ABS_IN")
        need(m.subtract(outer, inner, keep_originals=False), "ABS_SUB")
        m[outer.name].name = "ABSORBER"; m["ABSORBER"].material_name = "ABSORBER_15_10"
        mp = 1 if a.build_only else a.max_passes
        setup = need(h.create_setup("Setup1", Frequency="5.68GHz", MaximumPasses=a.max_passes, MaxDeltaS=0.02, MinimumPasses=a.min_passes), "SETUP")
        setup.props["BasisOrder"] = 1
        setup.props["DrivenSolverType"] = "Direct"
        need(setup.update(), "SETUP_UPDATE")
        if a.sweep:
            h.create_linear_count_sweep(setup="Setup1", units="GHz", start_frequency=4.5, stop_frequency=7.0, num_of_freq_points=26,
                                        name="Sweep1", sweep_type="Interpolating", save_fields=False)
        need(h.insert_infinite_sphere(theta_start=0, theta_stop=180, theta_step=2, phi_start=0, phi_stop=360, phi_step=2, name="Global"), "SPHERE")
        st["objects"] = len(m.object_names)
        st["boundaries"] = [b.name for b in h.boundaries]
        v = h.validate_simple(str(out / "validation.log"))
        st["validation"] = v
        need(v == 1 or v is True, "VALIDATION_FAILED")
        need(h.save_project(), "SAVE")
        if a.build_only:
            st["phase"] = "BUILD_COMPLETE"
        else:
            st["phase"] = "SOLVING"; st["solve_start_utc"] = now()
            status_path.write_text(json.dumps(st, indent=1, default=str))
            t0 = time.monotonic()
            need(h.analyze_setup("Setup1", cores=a.cores, use_auto_settings=False, blocking=True), "SOLVE_FAILED")
            st["wall_s"] = round(time.monotonic() - t0, 1)
            need(h.save_project(), "POST_SAVE")
            try:
                sw = "Setup1 : Sweep1" if a.sweep else "Setup1 : LastAdaptive"
                sd = h.post.get_solution_data(expressions=["S(Port1,Port1)"], setup_sweep_name=sw)
                fr, re = sd.get_expression_data(formula="real"); _, im = sd.get_expression_data(formula="imag")
                st["S11"] = {"f_ghz": [float(x) for x in fr], "re": [float(x) for x in re], "im": [float(x) for x in im]}
            except BaseException as e:
                st["S11_error"] = repr(e)
            try:
                prof = Path(h.export_profile("Setup1", output_file=str(out / "solve.prof"))); st["profile_path"] = str(prof)
            except BaseException as e:
                st["profile_error"] = repr(e)
            st["phase"] = "COMPLETE"
        status_path.write_text(json.dumps(st, indent=1, default=str))
        print(st["phase"], flush=True)
    except BaseException:
        st["phase"] = "FAILED"; st["traceback"] = traceback.format_exc()
        status_path.write_text(json.dumps(st, indent=1, default=str))
        print(st["traceback"], file=sys.stderr, flush=True)
        raise
    finally:
        done.set()
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
