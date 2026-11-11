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


def rows_60():
    return [dict(r, run_id=r["run_id"] + "_6p0", solve_ghz=6.0) for r in BS.cases()]


def rows_ctrl():
    base = {"kind": "BOARD", "lid_mm": 0.0, "island_l": 0.0, "n_gate": 0, "modes": 1, "renorm": True, "gate": False, "state": "CTRL"}
    rs0, xs0 = G.STATES["OFF"]
    return [dict(base, run_id="BOARD_C1", control="C1", Rs=rs0, Xs=xs0), dict(base, run_id="BOARD_C2", control="C2", Rs=rs0, Xs=xs0),
            dict(base, run_id="BOARD_C3", control="C3", Rs=377.0, Xs=0.0)]


def dev_ctrl_rows():
    out = []
    for r in rows_ctrl():
        out.append(dict(r, run_id=r["run_id"].replace("BOARD_", "BOARD_DEV") , r_board=45.0, gate_cols=[1, 2], wall_x=17.4,
                        probe_x=38.0, probe_half_y=20.0, max_passes=2, flux_xs=[38.0, 25.0]))
    return out


def all_rows():
    return ({r["run_id"]: r for r in BS.cases()} | {r["run_id"]: r for r in rows_60()} | {r["run_id"]: r for r in rows_ctrl()} | {r["run_id"]: r for r in dev_ctrl_rows()} | {"BOARD_DEV": dev_row()})


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


WALL_T = 0.5
LINE_DEFS = {   # name: ((x0, y0), (x1, y1), z_above_copper)
    "PROBE_P3": ((BS.PROBE_X, -BS.PROBE_HALF_Y), (BS.PROBE_X, BS.PROBE_HALF_Y), 3.0),
    "PROBE_M3": ((-BS.PROBE_X, -BS.PROBE_HALF_Y), (-BS.PROBE_X, BS.PROBE_HALF_Y), 3.0),
    "PROBE_P1": ((BS.PROBE_X, -BS.PROBE_HALF_Y), (BS.PROBE_X, BS.PROBE_HALF_Y), 1.0),
    "PROBE_M1": ((-BS.PROBE_X, -BS.PROBE_HALF_Y), (-BS.PROBE_X, BS.PROBE_HALF_Y), 1.0),
    "RAD_P": ((20.0, 0.0), (95.0, 0.0), 3.0),
    "RAD_M": ((-20.0, 0.0), (-95.0, 0.0), 3.0),
}


def make_lines(h, row):
    px, py = row.get("probe_x", BS.PROBE_X), row.get("probe_half_y", BS.PROBE_HALF_Y)
    existing = set(h.modeler.line_names) if hasattr(h.modeler, "line_names") else set()
    for name, ((x0, y0), (x1, y1), zc) in LINE_DEFS.items():
        if name.startswith("PROBE_"):
            sgn = 1 if name.endswith("P1") or name.endswith("P3") else -1
            x0, x1, y0, y1 = sgn * px, sgn * px, -py, py
        if h.modeler.does_object_exist(name) if hasattr(h.modeler, "does_object_exist") else name in h.modeler.object_names:
            continue
        GH.require(h.modeler.create_polyline([[x0, y0, G.CU_T + zc], [x1, y1, G.CU_T + zc]], name=name, non_model=True),
                   "line " + name)


def ensure_exprs(h, anchor):
    fc = h.post.fields_calculator
    base = {"description": "phasor field magnitude", "design_type": ["HFSS"], "fields_type": ["Fields"], "solution_type": "",
            "primary_sweep": "Freq", "assignment": "", "assignment_type": ["Line"], "report": ["Data Table", "Rectangular Plot"]}
    for name, ops in (("CMag_E", ["Fundamental_Quantity('E')", "Operation('Smooth')", "Operation('CmplxMag')", "Operation('Mag')"]),
                      ("Mag_Ez", ["Fundamental_Quantity('E')", "Operation('Smooth')", "Operation('CmplxMag')", "Operation('ScalarZ')"])):
        if not fc.is_expression_defined(name):
            GH.require(fc.add_expression(dict(base, name=name, operations=ops), anchor), "named expression " + name)


def line_series(h, line, qty, fs, setup_sol):
    sd = GH.require(h.post.get_solution_data(expressions=[qty], setup_sweep_name=setup_sol, report_category="Fields",
                                             context=line, primary_sweep_variable="Distance",
                                             variations={"Freq": [f"{fs:g}GHz"], "Phase": ["0deg"]}), f"{qty} on {line}")
    x, y = sd.get_expression_data(sd.expressions[0], formula="real")
    return [float(v) for v in x], [float(v) for v in y]


def field_metrics(h, fs, setup_sol):
    """Phasor |E| (CMag_E), |Ez| (Mag_Ez) and the phase-0 Mag_E on all lines; radial profiles on a 1 mm grid."""
    import numpy as np
    ensure_exprs(h, "PROBE_P3")
    out = {"lines": {}, "radial": {}}
    for name in ("PROBE_P3", "PROBE_M3", "PROBE_P1", "PROBE_M1"):
        rec = {}
        for q in ("CMag_E", "Mag_E") + (("Mag_Ez",) if name.endswith("1") else ()):
            x, y = line_series(h, name, q, fs, setup_sol)
            rec[q] = {"mean": float(np.nanmean(y)), "min": float(np.nanmin(y)), "max": float(np.nanmax(y)), "n": len(y), "n_nan": int(np.isnan(np.array(y)).sum())}
        out["lines"][name] = rec
    for name in ("RAD_P", "RAD_M"):
        d = np.arange(0, 76, 1.0)
        rec = {"x_mm": [(20.0 + k) * (1 if name == "RAD_P" else -1) for k in d]}
        for q in ("CMag_E", "Mag_E"):
            x, y = line_series(h, name, q, fs, setup_sol)
            rec[q] = [None if v != v else float(v) for v in np.interp(d, np.array(x), np.array(y))]
        out["radial"][name] = rec
    L = out["lines"]
    for z in ("3", "1"):
        for q in ("CMag_E", "Mag_E"):
            out[f"behind_over_mirror_dB_{q}_z{z}"] = 20 * math.log10(L[f"PROBE_P{z}"][q]["mean"] / L[f"PROBE_M{z}"][q]["mean"])
    return out


FLUX_TMP = GH.WORK / "flux_tmp"
FLUX_BANDS = {"lo": (0.0, 15.0), "hi": (15.0, 39.8)}      # mm above the copper top (hi stops 0.2 mm short of the radiation boundary)


def create_nonmodel_yz(h, name, x, y0, z0, w, hgt):
    """Non-model rectangle in a YZ plane (flag NonModel# at creation so the solution is not invalidated)."""
    args = ["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", f"{x:.9f}mm", "YStart:=", f"{y0:.9f}mm",
            "ZStart:=", f"{z0:.9f}mm", "Width:=", f"{w:.9f}mm", "Height:=", f"{hgt:.9f}mm", "WhichAxis:=", "X"]
    att = ["NAME:Attributes", "Name:=", name, "Flags:=", "NonModel#", "Color:=", "(143 175 143)", "Transparency:=", 0,
           "PartCoordinateSystem:=", "Global", "UDMId:=", "", "MaterialValue:=", '""', "SurfaceMaterialValue:=", '""',
           "SolveInside:=", False, "ShellElement:=", False, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", True,
           "UseMaterialAppearance:=", False, "IsLightweight:=", False]
    h.modeler.oeditor.CreateRectangle(args, att)


def surface_flux(h, sheet, sgn, fs, setup_sol):
    """Integral over `sheet` of Real(Poynting) . (sgn*x_hat) in W (field calculator; sheet is a non-model rectangle)."""
    f = h.post.ofieldsreporter
    f.CalcStack("clear")
    f.EnterQty("Poynting")
    f.CalcOp("Real")
    f.EnterVector([str(sgn), "0", "0"])
    f.CalcOp("Dot")
    f.EnterSurf(sheet)
    f.CalcOp("Integrate")
    FLUX_TMP.mkdir(parents=True, exist_ok=True)
    tmp = FLUX_TMP / f"flux_{sheet}_{os.getpid()}.fld"
    f.CalculatorWrite(str(tmp), ["Solution:=", setup_sol], ["Freq:=", f"{fs:g}GHz", "Phase:=", "0deg"])
    val = float(tmp.read_text().strip().splitlines()[-1])
    tmp.unlink()
    f.CalcStack("clear")
    return val


def flux_metrics(h, fs, setup_sol, xs=(75.0, 45.0, 25.0), half_y=40.0):
    """Real Poynting flux (W, 1 W incident) through non-model YZ rectangles: +x flux at x=+X, -x flux at x=-X."""
    fc = h.post.fields_calculator
    existing = set(h.modeler.oeditor.GetObjectsInGroup("Sheets")) | set(h.modeler.oeditor.GetObjectsInGroup("Non Model"))
    out = {"planes": {}, "note": "Real(Poynting) . (+-x) integrated over y +-%g mm; lo = 0..15 mm, hi = 15..39.8 mm above copper" % half_y}
    for x in xs:
        for side, sgn in (("P", 1), ("M", -1)):
            for band, (za, zb) in FLUX_BANDS.items():
                obj_name, ex_name = f"FLXS_{side}{x:g}_{band}", f"Flux_{side}{x:g}_{band}"
                if obj_name not in existing:
                    create_nonmodel_yz(h, obj_name, sgn * x, -half_y, G.CU_T + za, 2 * half_y, zb - za)
                    existing.add(obj_name)
                y_ = [surface_flux(h, obj_name, sgn, fs, setup_sol)]
                out["planes"][f"{side}{x:g}_{band}"] = float(y_[0])
    for x in xs:
        for band in FLUX_BANDS:
            p_, m_ = out["planes"][f"P{x:g}_{band}"], out["planes"][f"M{x:g}_{band}"]
            out[f"ratio_dB_P_over_M_{x:g}_{band}"] = (10 * math.log10(p_ / m_) if p_ > 0 and m_ > 0 else None)
    return out


def build(h, row, out: Path):
    require = GH.require
    r_board = row.get("r_board", BS.R_BOARD)
    gate_cols = tuple(row.get("gate_cols", BS.GATE_COLS))
    ctrl = row.get("control")
    cs = BS.cells(r_board, gate_cols, gate=row["gate"])
    if ctrl:                                   # positive controls: discs of the gate band removed
        cs = [c for c in cs if c["i"] not in gate_cols]
    cols = BS.columns(cs)
    exp = (BS.expected(r_board, gate_cols)["gate" if row["gate"] else "ref"] if not ctrl else
           {"discs": len(cs), "gate_cells": 0, "diode_sheets": 0, "column_nets": 0})
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
    wall = None
    if ctrl:
        wx, wt = row.get("wall_x", 46.4), WALL_T
        half_w = math.sqrt(r_board ** 2 - (wx + wt / 2) ** 2) - 0.1         # chord, 0.1 mm short of the rim
        top = (0.0 if ctrl == "C1" else G.CU_T + 15.0)                        # z of the wall top (substrate top, or 15 mm above copper)
        if ctrl in ("C1", "C2"):
            wall = box("WALL", [f"{wx - wt / 2:.9f}mm", f"{-half_w:.9f}mm", f"{-G.RT_T:.9f}mm"],
                       [f"{wt:.9f}mm", f"{2 * half_w:.9f}mm", f"{top + G.RT_T:.9f}mm"])
            require(m.subtract(sub.name, wall.name, keep_originals=True), "wall slot in substrate")
        else:
            wall = require(m.create_rectangle("YZ", [f"{wx:.9f}mm", f"{-half_w:.9f}mm", f"{-G.RT_T:.9f}mm"],
                                              [f"{2 * half_w:.9f}mm", f"{top + G.RT_T:.9f}mm"], name="WALL_SHEET",
                                              material="vacuum"), "wall sheet")
            require(h.assign_impedance_to_sheet(wall.name, name="WALL_R377", resistance="Rs", reactance="Xs"),
                    "wall resistive sheet")

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
    cut = [sub.name, gnd.name] + ([wall.name] if (ctrl == "C2") else [])
    require(m.subtract(air.name, cut, keep_originals=True), "cut board solids from AIR")
    face = lambda pos: require(m.get_faceid_from_position(pos, assignment="AIR"), f"AIR face {pos}")
    faces = [face([half, 0, 0]), face([-half, 0, 0]), face([0, half, 0]), face([0, -half, 0]),
             face([30, 30, zmin]), face([30, 30, zmax])]
    require(h.assign_radiation_boundary_to_faces(faces, name="RAD"), "radiation boundary")

    # --- field probe lines (non-model): y-lines at x = +-75 (z = 1 and 3 mm above copper), radial lines y = 0
    make_lines(h, row)
    sweep_ok = True
    setup = require(h.create_setup(SETUP, setup_type="HFSSDriven", Frequency=f"{row.get('solve_ghz', BS.SOLVE_GHZ):g}GHz",
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
    if ctrl:
        require(("WALL" in objects) if ctrl in ("C1", "C2") else ("WALL_SHEET" in objects), "STOP: control wall object missing")
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
    fs = row.get("solve_ghz", BS.SOLVE_GHZ)
    res["solve_ghz"] = fs
    res["field_metrics"] = field_metrics(h, fs, f"{SETUP} : LastAdaptive")
    res["behind_over_mirror_dB"] = res["field_metrics"]["behind_over_mirror_dB_CMag_E_z3"]
    res["flux"] = flux_metrics(h, fs, f"{SETUP} : LastAdaptive", xs=tuple(row.get("flux_xs", (75.0, 45.0, 25.0))))
    trap = getattr(np, "trapezoid", None) or np.trapz
    # far field: U = |rETotal|^2 / (2 eta0) for 1 W incident (rETotal reported in mV); one cut per phi
    th = None
    Uphi = []
    phis = [2.0 * k for k in range(180)]            # 0..358 deg; 360 duplicates 0
    unit_scale = None
    for p_ in phis:
        sd2 = GH.require(post.get_solution_data(expressions=["rETotal"], setup_sweep_name=f"{SETUP} : LastAdaptive",
                                                report_category="Far Fields", context="FF", primary_sweep_variable="Theta",
                                                variations={"Freq": [f"{fs:g}GHz"], "Phi": [f"{p_:g}deg"]}),
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
    acc = res["S11"][[abs(r["f_ghz"] - fs) < 1e-6 for r in res["S11"]].index(True)]["accepted_fraction"]
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
    root = DEV if run_id.startswith("BOARD_DEV") else OUT
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
