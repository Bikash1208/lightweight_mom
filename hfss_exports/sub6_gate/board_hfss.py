"""Open-board gate test: HFSS builder / solver / extractor (one run per worker, own AEDT session).

    python board_hfss.py worker --run-id BOARD_REF --port 51101 --cores 21
Not imported by tests. Never touches AEDT processes it did not start.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback
from pathlib import Path

import board_spec as BS
import gate_spec as G
import gate_hfss as GH

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = REPO / "results" / "sub6_gate" / "board"
ROWS = OUT / "rows"
FAILED = OUT / "failed"
DEV = HERE / "gate_work" / "board_dev"
SETUP, SWEEP = "Board", "Sweep"
TAG = time.strftime("%m%d%H%M%S")
ROW_TIMEOUT_S = 24 * 3600
ETA0 = 376.730313668


def dev_row():
    return {"kind": "BOARD", "run_id": "BOARD_DEV", "state": "OFF", "gate": True, "Rs": G.STATES["OFF"][0],
            "Xs": G.STATES["OFF"][1], "lid_mm": 0.0, "island_l": 0.0, "n_gate": 0, "modes": 1, "renorm": True,
            "r_board": 45.0, "gate_cols": [1, 2], "probe_x": 38.0, "probe_half_y": 20.0, "max_passes": 2}


def all_rows():
    return {r["run_id"]: r for r in BS.cases()} | {"BOARD_DEV": dev_row()}


def nm(prefix, c):
    return f"{prefix}_{c['i']}_{c['j']}".replace("-", "m")


def set_and_read_back(h, row):
    spec = BS.variables_for(row)
    for k, (val, unit) in spec.items():
        h[k] = f"{val:.12g}{unit}"
    rb, bad = {}, {}
    for k, (val, unit) in spec.items():
        got = h.get_evaluated_value(k)
        rb[k] = got if isinstance(got, (int, float)) and not isinstance(got, bool) else None
        if not GH.values_match(val, got, unit):
            bad[k] = {"intended": val, "read_back": got}
    if bad:
        raise RuntimeError(f"READBACK_MISMATCH {row['run_id']}: {json.dumps(bad)}")
    return rb


def build(h, row, out: Path):
    require = GH.require
    r_board = row.get("r_board", BS.R_BOARD)
    gate_cols = tuple(row.get("gate_cols", BS.GATE_COLS))
    cs = BS.cells(r_board, gate_cols, gate=row["gate"])
    cols = BS.columns(cs)
    exp = BS.expected(r_board, gate_cols)["gate" if row["gate"] else "ref"]
    dcg = BS.dc_gate(cs) if row["gate"] else None
    if row["gate"]:
        require(dcg["ok"], f"STOP: analytic DC gate failed {dcg}")
    for k, (val, unit) in BS.variables_for(row).items():
        h[k] = f"{val:.12g}{unit}"
    h["ground_bottom"] = "-rt_t-cu_t"
    h.modeler.model_units = "mm"
    mat = h.materials.add_material("RT5880_board")
    mat.permittivity = 2.2
    mat.dielectric_loss_tangent = 0.0009
    m = h.modeler

    def box(name, origin, size, material="pec"):
        return require(m.create_box(origin, size, name=name, material=material), "create " + name)

    def cyl(name, x, y, z, r, hgt, material="pec"):
        return require(m.create_cylinder("Z", [x, y, z], r, hgt, name=name, material=material), "create " + name)

    # --- board stack: substrate, solid ground with the probe opening, probe pin
    sub = cyl("SUBSTRATE", "0mm", "0mm", "-rt_t", "r_board", "rt_t", "RT5880_board")
    gnd = cyl("GROUND", "0mm", "0mm", "ground_bottom", "r_board", "cu_t")
    hole = cyl("GHOLE", "0mm", "0mm", "ground_bottom-cu_t", "hole_r", "3*cu_t", "vacuum")
    require(m.subtract(gnd.name, hole.name, keep_originals=False), "ground opening")
    pin_h = BS.PIN_ABOVE_GROUND_TOP + G.RT_T + G.CU_T + 0.0     # from ground bottom (z=-3.21) to 12.1 above ground top
    pin_top = -G.RT_T + BS.PIN_ABOVE_GROUND_TOP
    pin = cyl("PIN", "0mm", "0mm", "ground_bottom", "pin_r", f"{pin_top - BS.GROUND_BOTTOM:.9f}mm")
    require(m.subtract(sub.name, pin.name, keep_originals=True), "pin hole in substrate")

    # --- copper: static discs, and gate cells
    for c in cs:
        if c["gate"]:
            continue
        xc, yc = f"{c['x']:.9f}mm", f"{c['y']:.9f}mm"
        d = cyl(nm("DISC", c), xc, yc, "0mm", "disc_r", "cu_t")
        s = box(nm("SLOT", c), [f"{xc}-slot_w/2", f"{yc}-slot_l/2", "-cu_t"], ["slot_w", "slot_l", "3*cu_t"], "vacuum")
        require(m.subtract(d.name, s.name, keep_originals=False), "slot " + d.name)
    sheets = []
    column_names = []
    for i, col in cols.items():
        parts = {"L": [], "R": []}
        for c in col:
            xc, yc = f"{c['x']:.9f}mm", f"{c['y']:.9f}mm"
            for tag, sign in (("L", -1), ("R", 1)):
                d = cyl(nm("SD" + tag, c), xc, yc, "0mm", "disc_r", "cu_t")
                if sign < 0:
                    origin = [f"{xc}-split_w/2", f"{yc}-period/2-0.1mm", "-cu_t"]
                else:
                    origin = [f"{xc}-period/2-0.1mm", f"{yc}-period/2-0.1mm", "-cu_t"]
                cut = box(nm("SDCUT" + tag, c), origin, ["period/2+split_w/2+0.1mm", "period+0.2mm", "3*cu_t"], "vacuum")
                require(m.subtract(d.name, cut.name, keep_originals=False), "split " + d.name)
                tx = f"{xc}-split_w/2" if sign < 0 else f"{xc}+diode_gap/2"
                tab = box(nm("SDTAB" + tag, c), [tx, f"{yc}-diode_w/2", "0mm"],
                          ["split_w/2-diode_gap/2", "diode_w", "cu_t"])
                parts[tag] += [d.name, tab.name]
            sh = require(m.create_rectangle("XY", [f"{xc}-diode_gap/2", f"{yc}-diode_w/2", "cu_t/2"],
                                            ["diode_gap", "diode_w"], name=nm("DIODE", c), material="vacuum"), "diode sheet")
            sheets.append(sh.name)
        for tag, sign in (("L", -1), ("R", 1)):
            xl, lo, hi = BS.line_span(col, sign, r_board)
            line = box(f"LINE_{tag}_{i}", [f"{xl - G.LINE_W / 2:.9f}mm", f"{lo:.9f}mm", "0mm"],
                       ["line_w", f"{hi - lo:.9f}mm", "cu_t"])
            require(m.unite(parts[tag] + [line.name]), f"column {i} {tag} net union")
            m[parts[tag][0]].name = f"COL{i}_{tag}"
            column_names.append(f"COL{i}_{tag}")
    if sheets:
        require(h.assign_impedance_to_sheet(sheets, name="BAR90_02LRH", resistance="Rs", reactance="Xs"),
                "diode impedance (one boundary on all sheets)")

    # --- probe port: annulus sheet r 0.635..2.05 at the ground bottom plane, 50 ohm lumped
    ps = require(m.create_circle("XY", ["0mm", "0mm", "ground_bottom"], "hole_r", name="PORT_SHEET", material="vacuum"),
                 "port sheet")
    pc = require(m.create_circle("XY", ["0mm", "0mm", "ground_bottom"], "pin_r", name="PORT_CUT", material="vacuum"),
                 "port cut")
    require(m.subtract(ps.name, pc.name, keep_originals=False), "annulus")
    zg = BS.GROUND_BOTTOM
    require(h.lumped_port("PORT_SHEET", integration_line=[[BS.PIN_R, 0.0, zg], [BS.HOLE_R, 0.0, zg]],
                          impedance=50, name="Probe", renormalize=True), "lumped probe port")

    # --- air box (board solids cut out), radiation on all six outer faces
    half = r_board + BS.PAD_LAT
    zmin, zmax = zg - BS.PAD_BOT, G.CU_T + BS.PAD_TOP
    air = box("AIR", [f"{-half:.9f}mm", f"{-half:.9f}mm", f"{zmin:.9f}mm"],
              [f"{2 * half:.9f}mm", f"{2 * half:.9f}mm", f"{zmax - zmin:.9f}mm"], "vacuum")
    require(m.subtract(air.name, [sub.name, gnd.name], keep_originals=True), "cut board solids from AIR")
    face = lambda pos: require(m.get_faceid_from_position(pos, assignment="AIR"), f"AIR face {pos}")
    faces = [face([half, 0, 0]), face([-half, 0, 0]), face([0, half, 0]), face([0, -half, 0]),
             face([30, 30, zmin]), face([30, 30, zmax])]
    require(h.assign_radiation_boundary_to_faces(faces, name="RAD"), "radiation boundary")

    # --- field probe lines (non-model) at z = 3 mm above copper
    px, py = row.get("probe_x", BS.PROBE_X), row.get("probe_half_y", BS.PROBE_HALF_Y)
    pz = G.CU_T + BS.PROBE_Z_ABOVE_CU
    for nm_, x in (("PROBE_P", px), ("PROBE_M", -px)):
        require(m.create_polyline([[x, -py, pz], [x, py, pz]], name=nm_, non_model=True), "probe line " + nm_)

    sweep_ok = True
    setup = require(h.create_setup(SETUP, setup_type="HFSSDriven", Frequency=f"{BS.SOLVE_GHZ}GHz",
                                   MaximumPasses=row.get("max_passes", BS.MAX_PASSES), MinimumPasses=1,
                                   MaxDeltaS=BS.MAX_DELTA_S), "driven setup")
    require(h.create_linear_count_sweep(setup=SETUP, unit="GHz", start_frequency=BS.SWEEP[0], stop_frequency=BS.SWEEP[1],
                                        num_of_freq_points=len(BS.FREQS), name=SWEEP, save_fields=False,
                                        sweep_type="Discrete"), "discrete sweep")
    require(h.insert_infinite_sphere(theta_start=0, theta_stop=180, theta_step=2, phi_start=0, phi_stop=360,
                                     phi_step=2, name="FF"), "infinite sphere")
    info, ok = h.validate_full_design(output_dir=str(out))
    objects = list(m.object_names)
    counts = {"discs": sum(1 for o in objects if o.startswith("DISC_")),
              "gate_cells": len({o[len("SDL_"):] for o in objects if False}) if False else sum(len(c) for c in cols.values()),
              "diode_sheets": sum(1 for o in objects if o.startswith("DIODE_")),
              "column_nets": sum(1 for o in objects if o.startswith("COL"))}
    summary = {"status": "BUILT" if ok else "VALIDATION_FAIL", "aedt": h.desktop_class.aedt_version_id,
               "expected_counts": exp, "object_counts": counts, "n_objects": len(objects),
               "boundaries": [b.name for b in h.boundaries], "dc_gate": dcg, "validation": info,
               "setup_props": {k: setup.props.get(k) for k in ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS",
                                                               "BasisOrder")},
               "columns": {str(i): {"x": col[0]["x"], "rows": [col[0]["j"], col[-1]["j"]], "cells": len(col),
                                    "line_L": BS.line_span(col, -1, r_board), "line_R": BS.line_span(col, 1, r_board)}
                           for i, col in cols.items()}}
    GH.write_json(out / "build_summary.json", summary)
    require(ok, "native design validation failed")
    for k in exp:
        require(counts[k] == exp[k], f"STOP: object count {k}: built {counts[k]} != spec {exp[k]}")
    if row["gate"]:
        require(all(f"COL{i}_{t}" in objects for i in cols for t in "LR"), "STOP: column net object missing/merged")
    require(h.save_project(), "save project")
    return summary


def _col(sd, e, formula="real"):
    return [float(v) for v in sd.get_expression_data(e, formula=formula)[1]]


def extract(h, row, out: Path):
    import numpy as np
    post = h.post
    res = {}
    # S11 sweep
    sd = GH.require(post.get_solution_data(expressions=["S(1,1)"], setup_sweep_name=f"{SETUP} : {SWEEP}",
                                           report_category="Modal Solution Data"), "S11")
    f = [float(v) for v in sd.primary_sweep_values]
    re_, im_ = _col(sd, "S(1,1)", "real"), _col(sd, "S(1,1)", "imag")
    res["S11"] = [{"f_ghz": a, "S11": [b, c], "abs_S11": math.hypot(b, c),
                   "accepted_fraction": 1 - (b * b + c * c)} for a, b, c in zip(f, re_, im_)]
    # field lines at 5.8 GHz (adaptive solution)
    fields = {}
    for nm_ in ("PROBE_P", "PROBE_M"):
        gd = GH.require(post.get_solution_data(expressions=["Mag_E"], setup_sweep_name=f"{SETUP} : LastAdaptive",
                                               report_category="Fields", context=nm_,
                                               primary_sweep_variable="Distance",
                                               variations={"Freq": [f"{BS.SOLVE_GHZ}GHz"], "Phase": ["0deg"]}),
                        "Mag_E on " + nm_)
        x = [float(v) for v in gd.primary_sweep_values]
        y = _col(gd, "Mag_E", "real")
        fields[nm_] = {"n": len(y), "mean_mag_E_V_per_m": float(np.mean(y)), "min": float(np.min(y)), "max": float(np.max(y)),
                       "dist_first_last_mm": [x[0], x[-1]]}
    res["field_lines"] = fields
    res["behind_over_mirror_dB"] = 20 * math.log10(fields["PROBE_P"]["mean_mag_E_V_per_m"] / fields["PROBE_M"]["mean_mag_E_V_per_m"])
    trap = getattr(np, "trapezoid", None) or np.trapz
    # far field: U = |rETotal|^2 / (2 eta0) for 1 W incident (rETotal reported in mV); one cut per phi
    th = None
    Uphi = []
    phis = [2.0 * k for k in range(180)]            # 0..358 deg; 360 duplicates 0
    unit_scale = None
    for p_ in phis:
        sd2 = GH.require(post.get_solution_data(expressions=["rETotal"], setup_sweep_name=f"{SETUP} : LastAdaptive",
                                                report_category="Far Fields", context="FF", primary_sweep_variable="Theta",
                                                variations={"Freq": [f"{BS.SOLVE_GHZ}GHz"], "Phi": [f"{p_:g}deg"]}),
                         f"far field phi {p_}")
        if unit_scale is None:
            unit = sd2.units_data.get("rETotal", "")
            unit_scale = {"mV": 1e-3, "V": 1.0, "uV": 1e-6, "kV": 1e3}[unit]
            res["rETotal_unit"] = unit
        x, y = sd2.get_expression_data("rETotal", formula="real")
        th = np.array([float(v) for v in x])
        Uphi.append((np.array([float(v) for v in y]) * unit_scale) ** 2 / (2 * ETA0))
    U = np.array(Uphi)                                # [phi, theta]
    dph = np.radians(2.0)
    wth = np.sin(np.radians(th))
    P_phi = trap(U * wth[None, :], np.radians(th), axis=1) * dph      # per phi sample, full theta 0..180
    up = th <= 90 + 1e-9
    P_phi_up = trap((U * wth[None, :])[:, up], np.radians(th[up]), axis=1) * dph
    ph = np.array(phis)
    wf = np.where((ph == 90) | (ph == 270), 0.5, np.where((ph < 90) | (ph > 270), 1.0, 0.0))
    wb = np.where((ph == 90) | (ph == 270), 0.5, np.where((ph > 90) & (ph < 270), 1.0, 0.0))
    p_tot, p_up = float(P_phi.sum()), float(P_phi_up.sum())
    pf, pb = float((P_phi_up * wf).sum()), float((P_phi_up * wb).sum())
    acc = res["S11"][[abs(r["f_ghz"] - BS.SOLVE_GHZ) < 1e-6 for r in res["S11"]].index(True)]["accepted_fraction"]
    res["far_field"] = {"n_theta": int(len(th)), "n_phi": len(phis), "radiated_power_W_per_1W_incident": p_tot,
                        "accepted_power_W_per_1W_incident": acc, "radiated_over_accepted": p_tot / acc,
                        "upper_hemisphere_power_W": p_up, "upper_fraction_of_radiated": p_up / p_tot,
                        "upper_phi_m90_90_fraction": pf / (pf + pb), "upper_phi_90_270_fraction": pb / (pf + pb),
                        "front_to_back_upper_dB": 10 * math.log10(pf / pb), "peak_U_W_per_sr": float(U.max())}
    return res


def solve_once(row, port, cores, attempt, out_dir, root):
    run_id = row["run_id"]
    proj_dir = GH.WORK / "projects" / f"{run_id}_{TAG}_a{attempt}"
    project = proj_dir / f"{run_id}.aedt"
    if project.exists():
        raise RuntimeError(f"refusing to reuse {project}")
    proj_dir.mkdir(parents=True, exist_ok=True)
    h = None
    try:
        h = GH.connect(project, port)
        (GH.WORK / "pids").mkdir(parents=True, exist_ok=True)
        (GH.WORK / "pids" / f"{run_id}_{TAG}_a{attempt}.pid").write_text(str(h.desktop_class.aedt_process_id) + "\n")
        bsum = build(h, row, out_dir)
        rb = set_and_read_back(h, row)
        h.save_project()
        t0 = time.time()
        GH.require(h.analyze_setup(SETUP, cores=cores, use_auto_settings=True, blocking=False), "solve failed to start")
        time.sleep(30)
        while h.desktop_class.are_there_simulations_running:
            if time.time() - t0 > ROW_TIMEOUT_S:
                h.stop_simulations(clean_stop=False)
                raise RuntimeError(f"ROW_TIMEOUT {run_id}")
            time.sleep(20)
        seconds = time.time() - t0
        prof, conv = out_dir / f"{run_id}.prof", out_dir / f"{run_id}.conv"
        h.export_profile(SETUP, output_file=str(prof))
        try:
            h.export_convergence(SETUP, output_file=str(conv))
        except Exception:
            pass
        data = extract(h, row, out_dir)
        return {"run_id": run_id, "row": row, "variables": {k: v[0] for k, v in BS.variables_for(row).items()},
                "readback_si": rb, "readback_ok": True, "build": bsum, "aedt": h.desktop_class.aedt_version_id,
                "cores": cores, "solve_seconds": seconds, "attempt": attempt, **data,
                "mesh_raw": {"profile_text": prof.read_text(errors="replace") if prof.exists() else None,
                             "convergence_text": conv.read_text(errors="replace") if conv.exists() else None},
                "label": "HFSS_SIMULATED / NOT_MEASURED"}
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
    row = all_rows()[run_id]
    root = DEV if run_id == "BOARD_DEV" else OUT
    final = (root / "rows" / f"{run_id}.json")
    if final.exists():
        print("already done", run_id)
        return 0
    errors = []
    for attempt in (1, 2):
        out_dir = GH.WORK / "out" / f"{run_id}_{TAG}_a{attempt}"
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            res = solve_once(row, port + 1000 * (attempt - 1), cores, attempt, out_dir, root)
            GH.write_json(final, res)
            print("DONE", run_id, flush=True)
            return 0
        except Exception as exc:
            errors.append({"attempt": attempt, "error": f"{type(exc).__name__}: {str(exc)[:600]}",
                           "trace_tail": traceback.format_exc()[-1500:]})
            print(f"ROW_FAILED {run_id} attempt {attempt}: {exc}", flush=True)
    GH.write_json(root / "failed" / f"{run_id}.json", {"run_id": run_id, "attempts": 2, "errors": errors, "row": row})
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
    sys.exit(worker(a.run_id, a.port, a.cores))


if __name__ == "__main__":
    main()
