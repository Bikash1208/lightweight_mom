"""SUB6 gate test: driven-modal HFSS strip model builder + solver + extractor (one row per worker).

Launches its OWN new AEDT session(s) on an explicit gRPC port; never attaches to, or stops, a process it did not
start. Not imported by the tests.  Usage:
    python gate_hfss.py worker --run-id REF_H20 --port 50301 --cores 8
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time
import traceback
from pathlib import Path

import gate_spec as S

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "gate_work"
ROWS = REPO / "results" / "sub6_gate" / "rows"
FAILED = REPO / "results" / "sub6_gate" / "failed"
AEDT_VERSION = os.environ.get("SUB6_AEDT_VERSION", "2025.1")
DESIGN, SETUP, SWEEP = "Gate", "Gate", "Sweep"
ROW_TIMEOUT_S = 6 * 3600
TAG = time.strftime("%m%d%H%M%S")
UNIT_SI = {"mm": 1e-3, "ohm": 1.0, "": 1.0}


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def require(value, context):
    if value is False or value is None:
        raise RuntimeError(context)
    return value


def connect(project: Path, port: int):
    from ansys.aedt.core import Hfss
    return Hfss(project=str(project.resolve()), design=DESIGN, solution_type="Modal", version=AEDT_VERSION,
                non_graphical=True, new_desktop=True, close_on_exit=False, port=port)


def values_match(want, got, unit):
    if not isinstance(got, (int, float)) or isinstance(got, bool) or not math.isfinite(got):
        return False
    return any(math.isclose(got, c, rel_tol=1e-9, abs_tol=1e-12) for c in (want * UNIT_SI.get(unit, 1.0), want))


def set_and_read_back(h, row):
    spec = S.variables_for(row)
    for k, (val, unit) in spec.items():
        h[k] = f"{val:.12g}{unit}"
    rb, bad = {}, {}
    for k, (val, unit) in spec.items():
        got = h.get_evaluated_value(k)
        rb[k] = got if isinstance(got, (int, float)) and not isinstance(got, bool) else None
        if not values_match(val, got, unit):
            bad[k] = {"intended": val, "read_back": got}
    for k in ("ground_bottom", "lid_top"):
        rb[k] = h.get_evaluated_value(k)
    if bad:
        raise RuntimeError(f"READBACK_MISMATCH {row['run_id']}: {json.dumps(bad)}")
    return rb


def build(h, row, out: Path):
    """Build the strip. Copper/via/ground follow opt4_cell.py (G-ISL-noline) exactly; see FLAGS in the report."""
    cells, total = S.layout(row["n_gate"])
    gate_cells = [c for c in cells if c["kind"] == "gate"]
    spec = S.variables_for(row)
    for k, (val, unit) in spec.items():
        h[k] = f"{val:.12g}{unit}"
    h["ground_bottom"] = "-rt_t-cu_t"
    h["lid_top"] = "cu_t+lid_h"
    h.modeler.model_units = "mm"
    mat = h.materials.add_material("RT5880_gate")
    mat.permittivity = 2.2
    mat.dielectric_loss_tangent = 0.0009
    m = h.modeler
    x0 = -total / 2

    def box(name, origin, size, material="pec"):
        return require(m.create_box(origin, size, name=name, material=material), "create " + name)

    def cyl(name, x, y, z, r, hgt, material="pec"):
        return require(m.create_cylinder("Z", [x, y, z], r, hgt, name=name, material=material), "create " + name)

    for c in cells:
        i, xc = c["index"], f"{c['xc']:.9f}mm"
        disc = cyl(f"DISC_{i}", xc, "0mm", "0mm", "disc_r", "cu_t")
        slot = box(f"SLOT_{i}", [f"{xc}-slot_w/2", "-slot_l/2", "-cu_t"], ["slot_w", "slot_l", "3*cu_t"], "vacuum")
        require(m.subtract(disc.name, slot.name, keep_originals=False), f"slot subtraction {i}")
    substrate = box("RT5880", [f"{x0:.9f}mm", "-period/2", "-rt_t"], [f"{total:.9f}mm", "period", "rt_t"], "RT5880_gate")
    ground = box("GROUND", [f"{x0:.9f}mm", "-period/2", "ground_bottom"], [f"{total:.9f}mm", "period", "cu_t"])
    air = box("AIR", [f"{x0:.9f}mm", "-period/2", "0mm"], [f"{total:.9f}mm", "period", "cu_t+lid_h"], "vacuum")
    glue = []
    for c in gate_cells:
        i, xc = c["index"], f"{c['xc']:.9f}mm"
        vx = f"{xc}-slot_w/2+diode_gap+island_w/2"           # island centroid x = via x
        isl = box(f"ISLAND_{i}", [f"{xc}-slot_w/2+diode_gap", "-island_l/2", "0mm"], ["island_w", "island_l", "cu_t"])
        via = cyl(f"VIA_{i}", vx, "0mm", "-rt_t", "via_r", "rt_t+cu_t")
        pad = cyl(f"PAD_{i}", vx, "0mm", "0mm", "pad_r", "cu_t")
        require(m.subtract(substrate.name, via.name, keep_originals=True), f"drill via {i}")
        # validation rejects the united ground+island body overlapping AIR: cut the island/pad/via-top volume out of AIR
        require(m.subtract(air.name, isl.name, keep_originals=True), f"cut island volume from AIR {i}")
        glue += [isl.name, via.name, pad.name]
    if glue:
        require(m.unite([ground.name] + glue), "island + via + pad + solid ground united")
    sheets = []
    for c in gate_cells:
        sh = require(m.create_rectangle("XY", [f"{c['xc']:.9f}mm-slot_w/2", "-diode_w/2", "cu_t/2"],
                                       ["diode_gap", "diode_w"], name=f"DIODE_{c['index']}", material="vacuum"),
                     "diode sheet")
        sheets.append(sh.name)
    if sheets:
        require(h.assign_impedance_to_sheet(sheets, name="BAR90_02LRH", resistance="Rs", reactance="Xs"),
                "diode impedance")
    lid_top = S.CU_T + row["lid_mm"]
    face = lambda pos, obj: require(m.get_faceid_from_position(pos, assignment=obj), f"{obj} face {pos}")
    P2 = S.PERIOD / 2
    require(h.assign_perfect_e([face([0, 0, lid_top], "AIR")], name="PEC_LID"), "PEC lid")
    require(h.assign_perfect_e([face([0, 0, -S.RT_T - S.CU_T], "GROUND")], name="PEC_GROUND_BOTTOM"), "PEC ground bottom")
    # Mirror symmetry about y = 0 (decision 1): y faces are even-symmetry planes -> Perfect H, full height.
    # STOP if any copper (disc, island, pad, via) reaches a y face.
    for n in (c["index"] for c in cells):
        bb = m[f"DISC_{n}"].bounding_box          # xmin, ymin, zmin, xmax, ymax, zmax (mm)
        if bb[1] <= -P2 + 1e-6 or bb[4] >= P2 - 1e-6:
            raise RuntimeError(f"STOP: DISC_{n} touches a y face, bbox {bb}")
    if gate_cells and S.ISLAND_LENGTHS and row["island_l"] / 2 >= P2 - 1e-6:
        raise RuntimeError("STOP: island touches a y face")
    yfaces = [face([0, -P2, lid_top / 2], "AIR"), face([0, P2, lid_top / 2], "AIR")]
    yfaces += [face([0, -P2, -S.RT_T / 2], "RT5880"), face([0, P2, -S.RT_T / 2], "RT5880")]
    require(h.assign_perfect_h(yfaces, name="PMC_Y"), "perfect H y faces")
    ports = {}
    for n, x in ((1, -total / 2), (2, total / 2)):
        sheet = require(m.create_rectangle("YZ", [f"{x:.9f}mm", "-period/2", "-rt_t"], ["period", "rt_t+cu_t+lid_h"],
                                           name=f"PORT{n}_SHEET", material="vacuum"), "port sheet")
        port = h.wave_port(sheet.name, integration_line=[[f"{x:.9f}mm", "0mm", f"{-S.RT_T:.9f}mm"],
                                          [f"{x:.9f}mm", "0mm", f"{S.CU_T + row['lid_mm']:.9f}mm"]],
                           modes=1, impedance=50, name=f"P{n}", renormalize=True, deembed=0)
        ports[f"P{n}"] = bool(port)
        require(port, f"wave port {n} rejected")
    setup = require(h.create_setup(SETUP, setup_type="HFSSDriven", Frequency=f"{S.SOLVE_FREQ_GHZ}GHz",
                                   MaximumPasses=S.MAX_PASSES, MinimumPasses=S.MIN_PASSES,
                                   MaxDeltaS=S.MAX_DELTA_S), "driven setup")
    require(h.create_linear_count_sweep(setup=SETUP, unit="GHz", start_frequency=S.FREQ_START,
                                        stop_frequency=S.FREQ_STOP, num_of_freq_points=S.N_FREQ, name=SWEEP,
                                        save_fields=False, sweep_type="Discrete"), "discrete sweep")
    info, ok = h.validate_full_design(output_dir=str(out))
    objects = list(m.object_names)
    disc_names = sorted(n for n in objects if n.startswith("DISC_"))
    summary = {"status": "BUILT" if ok else "VALIDATION_FAIL", "aedt": h.desktop_class.aedt_version_id,
               "objects": objects, "boundaries": [b.name for b in h.boundaries], "ports": ports,
               "n_discs": len(disc_names), "n_gate": len(gate_cells), "total_length_mm": total,
               "mesh_operations": [], "setup_props": {k: setup.props.get(k) for k in
                                                       ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS",
                                                        "BasisOrder", "SolveType", "SaveFields")},
               "validation": info}
    write_json(out / "build_summary.json", summary)
    require(ok, "native design validation failed")
    # HFSS-side DC check: one united ground body, every disc a separate body, one diode sheet per gate cell.
    require(len(disc_names) == len(cells), f"disc bodies {len(disc_names)} != {len(cells)}")
    require("GROUND" in objects and all(f"DIODE_{c['index']}" in objects for c in gate_cells), "ground/diode objects")
    require(all(f"ISLAND_{c['index']}" not in objects for c in gate_cells), "island not united with ground")
    require(h.save_project(), "save project")
    return summary


# ---------------------------------------------------------------------------------------
def parse_profile(text):
    """Best-effort parse of the HFSS profile: adaptive passes (tets, delta S) kept alongside the raw lines."""
    passes = []
    for line in text.splitlines():
        m = re.search(r"Adaptive Pass\s+(\d+)", line, re.I)
        if m:
            passes.append({"pass": int(m.group(1)), "line": line.strip()})
    return passes


def _re(sd, e):
    return [float(v) for v in sd.get_expression_data(e, formula="real")[1]]


def _im(sd, e):
    return [float(v) for v in sd.get_expression_data(e, formula="imag")[1]]


def extract(h, row, out: Path, rb, cores, seconds):
    post = h.post
    sweep = f"{SETUP} : {SWEEP}"
    cats = {}
    for cat in ("Modal Solution Data",):
        try:
            cats[cat] = post.available_quantities_categories(report_category=cat)
        except Exception as exc:
            cats[cat] = f"ERR {exc}"
    expr = ["S(1,1)", "S(2,1)", "S(1,2)", "S(2,2)"]
    sd = require(post.get_solution_data(expressions=expr, setup_sweep_name=sweep, report_category="Modal Solution Data"),
                 "S-parameter solution data")
    freqs = [float(f) for f in sd.primary_sweep_values]
    if len(freqs) != S.N_FREQ or any(abs(a - b) > 1e-6 for a, b in zip(freqs, S.FREQS_GHZ)):
        raise RuntimeError(f"unexpected frequency list {freqs}")
    s = {e: [complex(a, b) for a, b in zip(_re(sd, e), _im(sd, e))] for e in expr}
    # port gamma
    gamma, gerr = {}, None
    found = []
    for c in (cats.get("Modal Solution Data") or []):
        try:
            found += [q for q in post.available_report_quantities(report_category="Modal Solution Data",
                                                                   quantities_category=c) if "amma" in q]
        except Exception:
            pass
    for gname in found + ["Gamma(P1:1)", "Gamma(1:1)", "Gamma(P1,1)", "Gamma(1,1)"]:
        try:
            gd = post.get_solution_data(expressions=[gname], setup_sweep_name=sweep, report_category="Modal Solution Data")
            if gd and len(_re(gd, gname)):
                gamma = {"expression": gname, "re": list(map(float, _re(gd, gname))),
                         "im": list(map(float, _im(gd, gname)))}
                break
        except Exception as exc:
            gerr = f"{gname}: {exc}"
    quantities = {}
    try:
        quantities = {c: post.available_report_quantities(report_category="Modal Solution Data", quantities_category=c)
                      for c in ("Gamma", "Port Gamma") if c in (cats.get("Modal Solution Data") or [])}
    except Exception as exc:
        quantities = {"ERR": str(exc)}
    prof = out / f"{row['run_id']}.prof"
    h.export_profile(SETUP, output_file=str(prof))
    conv = out / f"{row['run_id']}.conv"
    try:
        h.export_convergence(SETUP, output_file=str(conv))
    except Exception as exc:
        conv = None
    mesh_raw = {"profile_text": prof.read_text(errors="replace") if prof.exists() else None,
                "convergence_text": conv.read_text(errors="replace") if conv and conv.exists() else None}
    per_f = []
    for k, f in enumerate(freqs):
        s11, s21 = s["S(1,1)"][k], s["S(2,1)"][k]
        d = {"f_ghz": f, "S11": [s11.real, s11.imag], "S21": [s21.real, s21.imag],
             "S12": [s["S(1,2)"][k].real, s["S(1,2)"][k].imag], "S22": [s["S(2,2)"][k].real, s["S(2,2)"][k].imag],
             "absorbed_fraction": 1 - abs(s11) ** 2 - abs(s21) ** 2}
        if gamma:
            d["port_gamma_1_per_m"] = [gamma["re"][k], gamma["im"][k]]
        per_f.append(d)
    return {"per_frequency": per_f, "gamma_expression": gamma.get("expression"), "gamma_error": gerr,
            "report_categories": cats, "gamma_quantities": quantities, "gamma_candidates": found, "mesh_raw": mesh_raw}


def solve_once(row, port, cores, attempt, out_dir):
    run_id = row["run_id"]
    proj_dir = WORK / "projects" / f"{run_id}_{TAG}_a{attempt}"
    project = proj_dir / f"{run_id}.aedt"
    if project.exists():
        raise RuntimeError(f"refusing to reuse existing project {project}")
    proj_dir.mkdir(parents=True, exist_ok=True)
    h = None
    try:
        h = connect(project, port)
        (WORK / "pids").mkdir(parents=True, exist_ok=True)
        (WORK / "pids" / f"{run_id}_{TAG}_a{attempt}.pid").write_text(
            str(h.desktop_class.aedt_process_id) + "\n", encoding="utf-8")
        bsum = build(h, row, out_dir)
        rb = set_and_read_back(h, row)
        h.save_project()
        t0 = time.time()
        require(h.analyze_setup(SETUP, cores=cores, use_auto_settings=True, blocking=False), "solve failed to start")
        time.sleep(20)
        while h.desktop_class.are_there_simulations_running:
            if time.time() - t0 > ROW_TIMEOUT_S:
                try:
                    h.stop_simulations(clean_stop=False)
                except Exception as exc:
                    print("stop_simulations:", exc, flush=True)
                raise RuntimeError(f"ROW_TIMEOUT {run_id}")
            time.sleep(15)
        seconds = time.time() - t0
        data = extract(h, row, out_dir, rb, cores, seconds)
        result = {"run_id": run_id, "row": row, "variables": {k: v[0] for k, v in S.variables_for(row).items()},
                  "variable_units": {k: v[1] for k, v in S.variables_for(row).items()}, "readback_si": rb,
                  "readback_ok": True, "build": bsum, "aedt": h.desktop_class.aedt_version_id, "cores": cores,
                  "solve_seconds": seconds, "attempt": attempt, "dc_gate": S.dc_gate(row["n_gate"]) if row["n_gate"] else None,
                  "label": "HFSS_SIMULATED / NOT_MEASURED", **data}
        return result
    finally:
        if h is not None:
            try:
                h.save_project()
            except Exception:
                pass
            try:
                h.release_desktop(close_projects=True, close_desktop=True)
            except Exception:
                pass


def worker(run_id, port, cores):
    row = next(r for r in S.cases() if r["run_id"] == run_id)
    final = ROWS / f"{run_id}.json"
    if final.exists():
        print("already done", run_id)
        return 0
    errors = []
    for attempt in (1, 2):
        out_dir = WORK / "out" / f"{run_id}_{TAG}_a{attempt}"
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            res = solve_once(row, port + 1000 * (attempt - 1), cores, attempt, out_dir)
            write_json(final, res)
            print("DONE", run_id, flush=True)
            return 0
        except Exception as exc:
            errors.append({"attempt": attempt, "error": f"{type(exc).__name__}: {str(exc)[:600]}",
                           "trace_tail": traceback.format_exc()[-1200:]})
            print(f"ROW_FAILED {run_id} attempt {attempt}: {exc}", flush=True)
    write_json(FAILED / f"{run_id}.json", {"run_id": run_id, "attempts": 2, "errors": errors, "row": row})
    print("ROW_SKIPPED_AS_FAILED", run_id, flush=True)
    return 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("worker")
    w.add_argument("--run-id", required=True)
    w.add_argument("--port", type=int, required=True)
    w.add_argument("--cores", type=int, required=True)
    a = ap.parse_args()
    if a.cmd == "worker":
        sys.exit(worker(a.run_id, a.port, a.cores))


if __name__ == "__main__":
    main()
