"""Phase 28 (HFSS-first path): writes the two HFSS packages for the cell G and the G state library report.

    python scripts/research/phase28_export_g_packages.py            # needs results/phase28/f10/g_library/g_library.json

Outputs (all MODEL_ONLY / NOT_HFSS_VALIDATED, diode package inductance 0 nH assumed):
    hfss_exports/phase28_cell_G/        single-cell Floquet check: build_cell.py, set_state.py, mom_reference.csv, bare_grounded_slab_reference.csv, README.md
    hfss_exports/phase28_coded_small/   70 mm (across flats) hexagon: build_board.py, set_coding.py, layout.csv, codings/*.csv, mom_reference/*, README.md
    results/phase28/f10/g_library/G_LIBRARY.md
"""

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_g_geometry as geo  # noqa: E402
import phase28_small_board as sb  # noqa: E402

from src.dispersion.sheet_model import _y_slab  # noqa: E402
from src.frequency import angular_frequency, free_space_wavenumber  # noqa: E402
from src.hologram import coded_feasibility as cf  # noqa: E402
from src.pattern.vector_aperture_pattern import copol_component  # noqa: E402

LIB = ROOT / "results/phase28/f10/g_library/g_library.json"
CELL_DIR = ROOT / "hfss_exports/phase28_cell_G"
BOARD_DIR = ROOT / "hfss_exports/phase28_coded_small"
H_MM = 1.575
P_MM = 6.4
BEAM_NAMES = {(0.0, 0.0): "beam_000", (30.0, 225.0): "beam_030_225", (45.0, 225.0): "beam_045_225"}
VERIFIED = "VERIFIED in AEDT (3x3 smoke build, 2026-10-24)"
HEADER_CALLS = (
    "# AEDT calls: " + VERIFIED + ": CreatePolyline sheets, Unite, AssignPerfectE, AssignLumpedRLC (ON: R only, OFF: C only, one element type per sheet), far-field sphere.\n"
    "# Copied from the Phase 26 hex150 script (built the 150 mm board): substrate/ground/monopole/lumped port/radiation box/setup.\n"
    "# NOT yet run in AEDT: AssignMaster/AssignSlave with a scan angle, AssignFloquetPort, EditLumpedRLC / DeleteBoundaries. Each is in a try/except that prints what to do by hand.\n"
)

COMMON = '''import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()

def attrs(name, material, solve_inside, color="(255 128 0)", transp=0):
    return ["NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=", color, "Transparency:=", transp,
            "PartCoordinateSystem:=", "Global", "UDMId:=", "", "MaterialValue:=", '"' + material + '"',
            "SurfaceMaterialValue:=", '""', "SolveInside:=", solve_inside, "IsMaterialEditable:=", True,
            "UseMaterialAppearance:=", False, "IsLightweight:=", False]

def mm(v):
    return repr(float(v)) + "mm"

def poly_sheet(oEditor, name, pts, z, color="(255 128 0)"):
    pl = ["NAME:PolylineParameters", "IsPolylineCovered:=", True, "IsPolylineClosed:=", True,
          ["NAME:PolylinePoints"] + [["NAME:PLPoint", "X:=", mm(x), "Y:=", mm(y), "Z:=", z] for (x, y) in list(pts) + [pts[0]]],
          ["NAME:PolylineSegments"] + [["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", k, "NoOfPoints:=", 2] for k in range(len(pts))],
          ["NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm",
           "XSectionTopWidth:=", "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner"]]
    oEditor.CreatePolyline(pl, attrs(name, "vacuum", True, color))

def rlc_args(name, on, line, z):
    # ON = R only, OFF = C only: one single-type Lumped RLC boundary per sheet (the AEDT-verified form)
    return ["NAME:" + name + "_RLC", "Objects:=", [name],
            ["NAME:CurrentLine", "Start:=", [mm(line[0][0]), mm(line[0][1]), z], "End:=", [mm(line[1][0]), mm(line[1][1]), z]],
            "RLC Type:=", "Parallel",
            "UseResist:=", bool(on), "Resistance:=", R_ON if on else "0ohm",
            "UseInduct:=", False, "Inductance:=", "0nH",
            "UseCap:=", (not on), "Capacitance:=", "0pF" if on else C_OFF]
'''


def _lit(obj):
    """Python 2/3 compatible literal of nested lists/tuples of numbers/strings (repr; floats rounded to 6 decimals)."""
    if isinstance(obj, dict):
        return "{" + ", ".join(f"{_lit(k)}: {_lit(v)}" for k, v in obj.items()) + "}"
    if isinstance(obj, (list, tuple)):
        return "[" + ", ".join(_lit(v) for v in obj) + "]"
    if isinstance(obj, float):
        return repr(round(obj, 6))
    return repr(obj)


def _piece_table():
    return {k: [(round(x, 6), round(y, 6)) for x, y in v] for k, v in geo.cell_polygons_mm().items()}


def _neck_table():
    return [
        {
            "site": n["site"],
            "corners": [
                (n["rect_mm"][0], n["rect_mm"][2]),
                (n["rect_mm"][1], n["rect_mm"][2]),
                (n["rect_mm"][1], n["rect_mm"][3]),
                (n["rect_mm"][0], n["rect_mm"][3]),
            ],
            "line": n["current_line_mm"],
        }
        for n in geo.necks_mm()
    ]


# ----------------------------------------------------------------------------------------------------------- reference numbers
def _c(pair):
    return complex(pair[0], pair[1])


def _db(z):
    return 20 * math.log10(max(abs(z), 1e-300))


def _deg(z):
    return math.degrees(math.atan2(z.imag, z.real))


def circular(s):
    """`s` = 2 x 2 [received, incident] in the (TM, TE) order. Returns r_RR, r_RL, r_LR, r_LL (reflected, incident; propagation-relative; a metal plate is r_LR = r_RL = -1)."""
    a_r, a_l = np.array([1, 1j]) / math.sqrt(2), np.array([1, -1j]) / math.sqrt(2)
    b_r, b_l = np.array([1, -1j]) / math.sqrt(2), np.array([1, 1j]) / math.sqrt(2)
    return {
        "r_RR": complex(np.conj(b_r) @ s @ a_r),
        "r_RL": complex(np.conj(b_r) @ s @ a_l),
        "r_LR": complex(np.conj(b_l) @ s @ a_r),
        "r_LL": complex(np.conj(b_l) @ s @ a_l),
    }


def reference_rows(lib):
    rows = []
    best = 96
    for r in lib["rows"]:
        if r["kind"] == "normal" and r["mesh"] in (64, best):
            g = np.array([[_c(r["gamma_xy"][i][j]) for j in range(2)] for i in range(2)])  # [x, y] x [x, y]
            s = np.array([[g[0, 0], g[0, 1]], [g[1, 0], g[1, 1]]])  # TM = x, TE = y at phi = 0: order (TM, TE)
            ent = {"S_TM_TM": s[0, 0], "S_TE_TE": s[1, 1], "S_TM_TE": s[0, 1], "S_TE_TM": s[1, 0]}
            ent.update(circular(s))
            for q, z in ent.items():
                rows.append((q, r["state"], 0.0, 0.0, r["f_ghz"], r["mesh"], _db(z), _deg(z)))
        elif r["kind"] == "oblique":
            s11 = np.array([[_c(r["s11_te_tm"][i][j]) for j in range(2)] for i in range(2)])  # [TE, TM] x [TE, TM]
            s = np.array([[s11[1, 1], s11[1, 0]], [s11[0, 1], s11[0, 0]]])  # reorder to (TM, TE)
            ent = {"S_TM_TM": s[0, 0], "S_TE_TE": s[1, 1], "S_TM_TE": s[0, 1], "S_TE_TM": s[1, 0]}
            ent.update(circular(s))
            for q, z in ent.items():
                rows.append((q, r["state"], r["theta_deg"], r["phi_deg"], r["f_ghz"], r["mesh"], _db(z), _deg(z)))
    return sorted(rows, key=lambda t: (t[2], t[1], t[0], t[4], t[5]))


def bare_slab_rows():
    from src.hologram import tensor_cp_design as t26

    rows = []
    for f in (9.5, 9.75, 10.0, 10.25, 10.5):
        omega, k0 = angular_frequency(f * 1e9), free_space_wavenumber(f * 1e9)
        y = _y_slab(0.0, omega, t26.SUB.eps_complex, k0, t26.SUB.thickness)
        y0 = 1.0 / 376.730313668
        g = (y0 - y) / (y0 + y)
        rows.append((f, _db(g), _deg(g)))
    return rows


# ----------------------------------------------------------------------------------------------------------- cell package
def write_cell_package(lib):
    CELL_DIR.mkdir(parents=True, exist_ok=True)
    dv = geo.diode_values()
    head = (
        "# -*- coding: utf-8 -*-\n"
        "# HFSS (AEDT) IronPython script -- single G cell (6.4 mm rational-sector hub, 8 MADP diode sites), Floquet check. Generated by lwmom (Phase 28 HFSS-first path).\n"
        "# MODEL_ONLY / NOT_HFSS_VALIDATED. Diode package inductance 0 nH (NOT_AVAILABLE_FROM_CURRENT_DATA, assumed zero).\n" + HEADER_CALLS
    )
    consts = (
        f"STATE = 0                 # 0..3 = slot axis 0, 45, 90, 135 deg; change here, or run set_state.py on the built model\n"
        f"P = {P_MM}; H = {H_MM}; LID = 20.0\n"
        f"R_ON = \"{dv['r_on_ohm']}ohm\"; C_OFF = \"{dv['c_off_pf']}pF\"\n"
        f"PIECES = {_lit(_piece_table())}\n"
        f"NECKS = {_lit(_neck_table())}\n"
        f"STATE_ON = {_lit(geo.state_on_sites())}\n"
    )
    build = (
        head
        + COMMON
        + consts
        + '''
oProject = oDesktop.NewProject()
oProject.InsertDesign("HFSS", "phase28_cell_G", "DrivenModal", "")
oDesign = oProject.SetActiveDesign("phase28_cell_G")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])
oProject.GetDefinitionManager().AddMaterial(["NAME:RT5880_board", "CoordinateSystemType:=", "Cartesian",
    "BulkOrSurfaceType:=", 1, ["NAME:PhysicsTypes", "set:=", ["Electromagnetic"]],
    "permittivity:=", "2.2", "dielectric_loss_tangent:=", "0.0009"])
for n, v in (("scan_theta", "0deg"), ("scan_phi", "225deg")):
    oDesign.ChangeProperty(["NAME:AllTabs", ["NAME:LocalVariableTab", ["NAME:PropServers", "LocalVariables"],
        ["NAME:NewProps", ["NAME:" + n, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", v]]]])
Z = mm(H)
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", mm(-P / 2), "YPosition:=", mm(-P / 2), "ZPosition:=", "0mm",
    "XSize:=", mm(P), "YSize:=", mm(P), "ZSize:=", Z], attrs("Substrate", "RT5880_board", True, "(0 128 0)", 0.6))
oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", mm(-P / 2), "YStart:=", mm(-P / 2),
    "ZStart:=", "0mm", "Width:=", mm(P), "Height:=", mm(P), "WhichAxis:=", "Z"], attrs("Ground", "vacuum", True, "(128 128 128)"))
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", mm(-P / 2), "YPosition:=", mm(-P / 2), "ZPosition:=", Z,
    "XSize:=", mm(P), "YSize:=", mm(P), "ZSize:=", mm(LID)], attrs("Air", "air", True, "(143 175 143)", 0.9))
names = []
for nm in sorted(PIECES):
    poly_sheet(oEditor, "M_" + nm, PIECES[nm], Z)
    names.append("M_" + nm)
oEditor.Unite(["NAME:Selections", "Selections:=", ",".join(names)], ["NAME:UniteParameters", "KeepOriginals:=", False])
oEditor.ChangeProperty(["NAME:AllTabs", ["NAME:Geometry3DAttributeTab", ["NAME:PropServers", names[0]],
    ["NAME:ChangedProps", ["NAME:Name", "Value:=", "Metal"]]]])
for k, nk in enumerate(NECKS):
    poly_sheet(oEditor, "D%d" % k, nk["corners"], Z, "(255 0 0)")
oBnd = oDesign.GetModule("BoundarySetup")
oBnd.AssignPerfectE(["NAME:PerfE_Ground", "Objects:=", ["Ground"], "InfGroundPlane:=", False])
oBnd.AssignPerfectE(["NAME:PerfE_Metal", "Objects:=", ["Metal"], "InfGroundPlane:=", False])
for k, nk in enumerate(NECKS):
    oBnd.AssignLumpedRLC(rlc_args("D%d" % k, k in STATE_ON[STATE], nk["line"], Z))

def faces_at(x, y):
    ids = []
    for body, z in (("Substrate", H / 2.0), ("Air", H + LID / 2.0)):
        ids.append(int(oEditor.GetFaceByPosition(["NAME:FaceParameters", "BodyName:=", body,
                                                  "XPosition:=", mm(x), "YPosition:=", mm(y), "ZPosition:=", mm(z)])))
    return ids

# periodic walls (scan angle linked to the Floquet port) -- NOT yet run in AEDT; if it fails assign by hand (README step 2)
try:
    pairs = (("X", (-P / 2, 0.0), (P / 2, 0.0), ((-P / 2, -P / 2), (-P / 2, P / 2)), ((P / 2, -P / 2), (P / 2, P / 2))),
             ("Y", (0.0, -P / 2), (0.0, P / 2), ((-P / 2, -P / 2), (P / 2, -P / 2)), ((-P / 2, P / 2), (P / 2, P / 2))))
    for tag, mc, sc, mu, su in pairs:
        oBnd.AssignMaster(["NAME:Master" + tag, "Faces:=", faces_at(mc[0], mc[1]),
            "CoordSysVector:=", ["NAME:Vectors", "Origin:=", [mm(mu[0][0]), mm(mu[0][1]), "0mm"], "UPos:=", [mm(mu[1][0]), mm(mu[1][1]), "0mm"]],
            "ReverseV:=", False])
        oBnd.AssignSlave(["NAME:Slave" + tag, "Faces:=", faces_at(sc[0], sc[1]),
            "CoordSysVector:=", ["NAME:Vectors", "Origin:=", [mm(su[0][0]), mm(su[0][1]), "0mm"], "UPos:=", [mm(su[1][0]), mm(su[1][1]), "0mm"]],
            "ReverseU:=", False, "Master:=", "Master" + tag, "UseScanAngle:=", True, "Phi:=", "scan_phi", "Theta:=", "scan_theta"])
except Exception as exc:
    oDesktop.AddMessage(oProject.GetName(), "phase28_cell_G", 2, "Master/slave assignment failed (" + str(exc) + "): assign the X and Y pairs by hand with the scan angle, README step 2.")
try:
    top = int(oEditor.GetFaceByPosition(["NAME:FaceParameters", "BodyName:=", "Air", "XPosition:=", "0mm", "YPosition:=", "0mm", "ZPosition:=", mm(H + LID)]))
    oBnd.AssignFloquetPort(["NAME:FloquetPort1", "Faces:=", [top], "NumModes:=", 2, "RenormalizeAllTerminals:=", True,
        "UseScanAngle:=", True, "Phi:=", "scan_phi", "Theta:=", "scan_theta", "DoDeembed:=", True, "DeembedDist:=", mm(LID),
        ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", False, "CharImp:=", "Zpi"],
                       ["NAME:Mode2", "ModeNum:=", 2, "UseIntLine:=", False, "CharImp:=", "Zpi"]],
        "ShowReporterFilter:=", False, "ReporterFilter:=", [True, True]])
except Exception as exc:
    oDesktop.AddMessage(oProject.GetName(), "phase28_cell_G", 2, "Floquet port assignment failed (" + str(exc) + "): assign it by hand on the top face of Air, 2 modes, de-embed %.1f mm, README step 2." % LID)
oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven", ["NAME:Setup1", "Frequency:=", "10GHz", "MaxDeltaS:=", 0.01, "MaximumPasses:=", 15, "MinimumPasses:=", 2,
    "MinimumConvergedPasses:=", 2, "PercentRefinement:=", 30, "IsEnabled:=", True, "BasisOrder:=", 1, "DoLambdaRefine:=", True,
    "DoMaterialLambda:=", True, "SetLambdaTarget:=", False, "Target:=", 0.3333, "UseMaxTetIncrease:=", False])
oAna.InsertFrequencySweep("Setup1", ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearStep", "RangeStart:=", "9.5GHz", "RangeEnd:=", "10.5GHz",
    "RangeStep:=", "0.25GHz", "Type:=", "Discrete", "SaveFields:=", False, "SaveRadFields:=", False])
oEditor.FitAll()
oProject.SaveAs(oDesktop.GetProjectDirectory() + "/phase28_cell_G.aedt", True)
oDesktop.AddMessage(oProject.GetName(), "phase28_cell_G", 0, "Built state %d. Check the 8 RLC boundaries, the periodic pairs and the Floquet port, set scan_theta/scan_phi, Analyze." % STATE)
'''
    )
    (CELL_DIR / "build_cell.py").write_text(build, encoding="utf-8")
    setstate = (
        head.replace("single G cell", "set the diode state of the built single G cell")
        + COMMON
        + consts.replace("STATE = 0                 # 0..3 = slot axis 0, 45, 90, 135 deg; change here, or run set_state.py on the built model", "STATE = 1                 # 0..3: the state to put on the open model")
        + '''
oProject = oDesktop.GetActiveProject()
oDesign = oProject.GetActiveDesign()
oBnd = oDesign.GetModule("BoundarySetup")
Z = mm(H)
for k, nk in enumerate(NECKS):
    args = rlc_args("D%d" % k, k in STATE_ON[STATE], nk["line"], Z)
    try:
        oBnd.EditLumpedRLC("D%d_RLC" % k, args)
    except Exception:
        oBnd.DeleteBoundaries(["D%d_RLC" % k])
        oBnd.AssignLumpedRLC(args)
oDesktop.AddMessage(oProject.GetName(), oDesign.GetName(), 0, "State %d set (ON sites %s)." % (STATE, STATE_ON[STATE]))
'''
    )
    (CELL_DIR / "set_state.py").write_text(setstate, encoding="utf-8")
    rows = reference_rows(lib)
    with (CELL_DIR / "mom_reference.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "state", "theta_deg", "phi_deg", "f_ghz", "mesh", "mag_db", "phase_deg"])
        for r in rows:
            w.writerow([r[0], r[1], f"{r[2]:g}", f"{r[3]:g}", f"{r[4]:g}", r[5], f"{r[6]:.4f}", f"{r[7]:.3f}"])
    with (CELL_DIR / "bare_grounded_slab_reference.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["f_ghz", "r_mag_db", "r_phase_deg"])
        for f, m, p in bare_slab_rows():
            w.writerow([f"{f:g}", f"{m:.4f}", f"{p:.3f}"])
    return rows


# ----------------------------------------------------------------------------------------------------------- board package
def write_board_package(alpha_g):
    BOARD_DIR.mkdir(parents=True, exist_ok=True)
    (BOARD_DIR / "codings").mkdir(exist_ok=True)
    (BOARD_DIR / "mom_reference").mkdir(exist_ok=True)
    dv = geo.diode_values()
    preds = {}
    green = cf.green_function()
    pitch = sb.PITCH_MM * 1e-3
    for beam in sb.BEAMS:
        p = sb.predict(beam, alpha_g)
        preds[beam] = p
        name = BEAM_NAMES[beam]
        e = p["table_capped_at_G"]
        with (BOARD_DIR / "codings" / f"coding_{name}.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["cell_id", "x_mm", "y_mm", "state_0to3", "axis_deg"])
            for i, (x, y, s) in enumerate(zip(e["x_mm"], e["y_mm"], e["state_index"], strict=True)):
                w.writerow([i, f"{x:.4f}", f"{y:.4f}", s, 45 * s])
        # pattern cuts of the capped, quantised design (same currents as the prediction)
        x, y, d, march = cf.design_cells(
            cf.capped_table(sb.SupercellTable.from_files(sb.TABLE_GLOB), alpha_g), cf.hexagon_m(sb.circumradius_mm() * 1e-3), pitch, beam[0], beam[1]
        )
        idx, axis = cf.cp_states(d, 2)
        jx, jy = cf.cp_currents(d, axis)
        launch = cf.LAUNCH_BARE * float(march["eta_conv"])
        rows = []
        for phi in (beam[1], (beam[1] + 180.0) % 360.0):
            th = np.arange(0.0, 90.0 + 1e-9, 1.0)
            et, ep = cf.far_field(jx, jy, x, y, green, th, np.full_like(th, phi), pitch)
            rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
            lh = np.abs(copol_component(et, ep, "lhcp")) ** 2
            rows.append((phi, th, rh, lh))
        pk = max(float(np.max(r[2])) for r in rows)
        e["quantised_4_axes"]["_cut_reference_peak"] = pk
        d_ap = e["quantised_4_axes"]["d_ap_dbi"]
        with (BOARD_DIR / "mom_reference" / f"cut_{name}.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["phi_deg", "theta_deg", "rhcp_gain_dbi", "lhcp_gain_dbi"])
            for phi, th, rh, lh in rows:
                for t, a, b in zip(th, rh, lh, strict=True):
                    g0 = d_ap + 10 * math.log10(launch)
                    w.writerow([f"{phi:g}", f"{t:g}", f"{g0 + 10 * math.log10(max(a, 1e-30) / pk):.2f}", f"{g0 + 10 * math.log10(max(b, 1e-30) / pk):.2f}"])
        e["quantised_4_axes"].pop("_cut_reference_peak")
        (BOARD_DIR / "mom_reference" / f"predictions_{name}.json").write_text(json.dumps({**p, "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"]}, indent=1), encoding="utf-8")
    # layout from the 30/225 coding (cells are identical for all codings)
    e0 = preds[(30.0, 225.0)]["table_capped_at_G"]
    cells = [(i, x, y) for i, (x, y) in enumerate(zip(e0["x_mm"], e0["y_mm"], strict=True))]
    with (BOARD_DIR / "layout.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["cell_id", "x_mm", "y_mm"])
        for i, x, y in cells:
            w.writerow([i, f"{x:.4f}", f"{y:.4f}"])
    for beam in sb.BEAMS:
        a = preds[beam]["table_capped_at_G"]
        assert [round(v, 4) for v in a["x_mm"]] == [round(c[1], 4) for c in cells]
    init = preds[(30.0, 225.0)]["table_capped_at_G"]["state_index"]
    rad = sb.circumradius_mm()
    wx = sb.ACROSS_FLATS_MM + 10.0
    wy = 2 * rad + 10.0
    consts = (
        f"H = {H_MM}; P = {P_MM}\n"
        f"WX = {wx:.3f}; WY = {wy:.3f}          # substrate/ground: hexagon (70 mm across flats, circumradius {rad:.3f}) + 5 mm margin\n"
        f"R_ON = \"{dv['r_on_ohm']}ohm\"; C_OFF = \"{dv['c_off_pf']}pF\"\n"
        f"PIECES = {_lit(_piece_table())}\n"
        f"NECKS = {_lit(_neck_table())}\n"
        f"STATE_ON = {_lit(geo.state_on_sites())}\n"
    )
    cells_lit = "[" + ", ".join(f"({i}, {x:.4f}, {y:.4f})" for i, x, y in cells) + "]"
    head = (
        "# -*- coding: utf-8 -*-\n"
        f"# HFSS (AEDT) IronPython script -- small coded board: 70 mm (across flats) hexagon, {len(cells)} G cells, {8 * len(cells)} diode boundaries, centre monopole. Generated by lwmom (Phase 28 HFSS-first path).\n"
        "# MODEL_ONLY / NOT_HFSS_VALIDATED. Diode package inductance 0 nH (NOT_AVAILABLE_FROM_CURRENT_DATA, assumed zero). Builds the geometry ONCE; the diode states of any coding are set by set_coding.py.\n"
        + HEADER_CALLS
    )
    build = (
        head
        + "import math\n"
        + COMMON
        + consts
        + f"CELLS = {cells_lit}   # (cell_id, x, y) mm\n"
        + f"INIT_STATE = {_lit(init)}   # coding_beam_030_225.csv\n"
        + '''
oProject = oDesktop.NewProject()
oProject.InsertDesign("HFSS", "phase28_coded_small", "DrivenModal", "")
oDesign = oProject.SetActiveDesign("phase28_coded_small")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])
oProject.GetDefinitionManager().AddMaterial(["NAME:HMSA_substrate", "CoordinateSystemType:=", "Cartesian",
    "BulkOrSurfaceType:=", 1, ["NAME:PhysicsTypes", "set:=", ["Electromagnetic"]],
    "permittivity:=", "2.2", "dielectric_loss_tangent:=", "0.0009"])
def var(name, value):
    oDesign.ChangeProperty(["NAME:AllTabs", ["NAME:LocalVariableTab", ["NAME:PropServers", "LocalVariables"],
        ["NAME:NewProps", ["NAME:" + name, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", value]]]])
# the Phase 26 launch geometry (hex150_tensor_cp_single): pin 6.0 mm above the substrate = mono_h 7.575 mm from the ground, r_pin 0.635 mm, coax annulus r_out 2.05 mm
for n, v in (("mono_h", "7.575mm"), ("r_pin", "0.635mm"), ("r_out", "2.05mm"), ("pad_xy", "15mm"), ("pad_top", "30mm"), ("pad_bot", "15mm")):
    var(n, v)
FEED = (0.0, 0.0)
Z = mm(H)
x_lo = "(-" + mm(WX / 2) + ")"; y_lo = "(-" + mm(WY / 2) + ")"
def rect(name, x, y, z, w, h, color="(255 128 0)"):
    oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", x, "YStart:=", y, "ZStart:=", z,
        "Width:=", w, "Height:=", h, "WhichAxis:=", "Z"], attrs(name, "vacuum", True, color))
def circle(name, x, y, z, r, color="(255 128 0)"):
    oEditor.CreateCircle(["NAME:CircleParameters", "IsCovered:=", True, "XCenter:=", x, "YCenter:=", y, "ZCenter:=", z,
        "Radius:=", r, "WhichAxis:=", "Z", "NumSegments:=", "0"], attrs(name, "vacuum", True, color))
def subtract(blank, tools, keep):
    oEditor.Subtract(["NAME:Selections", "Blank Parts:=", blank, "Tool Parts:=", ",".join(tools)],
                     ["NAME:SubtractParameters", "KeepOriginals:=", keep])
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", x_lo, "YPosition:=", y_lo, "ZPosition:=", "0mm",
    "XSize:=", mm(WX), "YSize:=", mm(WY), "ZSize:=", Z], attrs("Substrate", "HMSA_substrate", True, "(0 128 0)", 0.6))
rect("Ground", x_lo, y_lo, "0mm", mm(WX), mm(WY), "(128 128 128)")
oBnd = oDesign.GetModule("BoundarySetup")
oEditor.CreateCylinder(["NAME:CylinderParameters", "XCenter:=", mm(FEED[0]), "YCenter:=", mm(FEED[1]), "ZCenter:=", "0mm",
    "Radius:=", "r_pin", "Height:=", "mono_h", "WhichAxis:=", "Z", "NumSides:=", "0"], attrs("Monopole1", "pec", False, "(255 255 0)"))
circle("Hole1", mm(FEED[0]), mm(FEED[1]), "0mm", "r_out")
circle("PortSheet1", mm(FEED[0]), mm(FEED[1]), "0mm", "r_out")
circle("PinCut1", mm(FEED[0]), mm(FEED[1]), "0mm", "r_pin")
subtract("PortSheet1", ["PinCut1"], False)
subtract("Ground", ["Hole1"], False)
subtract("Substrate", ["Monopole1"], True)
oBnd.AssignLumpedPort(["NAME:Port1", "Objects:=", ["PortSheet1"], "DoDeembed:=", False, "RenormalizeAllTerminals:=", True,
    ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,
        ["NAME:IntLine", "Start:=", [mm(FEED[0]) + "+r_pin", mm(FEED[1]), "0mm"], "End:=", [mm(FEED[0]) + "+r_out", mm(FEED[1]), "0mm"]],
        "AlignmentGroup:=", 0, "CharImp:=", "Zpi", "RenormImp:=", "50ohm"]],
    "ShowReporterFilter:=", False, "ReporterFilter:=", [True], "Impedance:=", "50ohm"])
# metal: nine pieces per cell, zero-thickness sheets at the pattern plane, united in blocks of 100
names = []
for cid, cx, cy in CELLS:
    for nm in sorted(PIECES):
        poly_sheet(oEditor, "M%d_%s" % (cid, nm), [(cx + x, cy + y) for (x, y) in PIECES[nm]], Z)
        names.append("M%d_%s" % (cid, nm))
blocks = []
for b in range(0, len(names), 100):
    grp = names[b:b + 100]
    if len(grp) > 1:
        oEditor.Unite(["NAME:Selections", "Selections:=", ",".join(grp)], ["NAME:UniteParameters", "KeepOriginals:=", False])
    blocks.append(grp[0])
if len(blocks) > 1:
    oEditor.Unite(["NAME:Selections", "Selections:=", ",".join(blocks)], ["NAME:UniteParameters", "KeepOriginals:=", False])
oEditor.ChangeProperty(["NAME:AllTabs", ["NAME:Geometry3DAttributeTab", ["NAME:PropServers", names[0]],
    ["NAME:ChangedProps", ["NAME:Name", "Value:=", "Metal"]]]])
# diode sheets D<cell>_<site>
for cid, cx, cy in CELLS:
    for k, nk in enumerate(NECKS):
        poly_sheet(oEditor, "D%d_%d" % (cid, k), [(cx + x, cy + y) for (x, y) in nk["corners"]], Z, "(255 0 0)")
oBnd.AssignPerfectE(["NAME:PerfE_Ground", "Objects:=", ["Ground"], "InfGroundPlane:=", False])
oBnd.AssignPerfectE(["NAME:PerfE_Metal", "Objects:=", ["Metal"], "InfGroundPlane:=", False])
failed = []
for cid, cx, cy in CELLS:
    for k, nk in enumerate(NECKS):
        try:
            line = [(cx + nk["line"][0][0], cy + nk["line"][0][1]), (cx + nk["line"][1][0], cy + nk["line"][1][1])]
            oBnd.AssignLumpedRLC(rlc_args("D%d_%d" % (cid, k), k in STATE_ON[INIT_STATE[cid]], line, Z))
        except Exception as e:
            failed.append("D%d_%d" % (cid, k))
if failed:
    oDesktop.AddMessage(oProject.GetName(), "phase28_coded_small", 2, "Lumped RLC assignment failed for %d sheets (first: %s)." % (len(failed), failed[0]))
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", x_lo + "-pad_xy", "YPosition:=", y_lo + "-pad_xy", "ZPosition:=", "-pad_bot",
    "XSize:=", mm(WX) + "+2*pad_xy", "YSize:=", mm(WY) + "+2*pad_xy", "ZSize:=", "pad_bot+" + Z + "+pad_top"],
    attrs("Air", "air", True, "(143 175 143)", 0.9))
oBnd.AssignRadiation(["NAME:Rad1", "Objects:=", ["Air"], "IsFssReference:=", False, "IsForPML:=", False])
oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven", ["NAME:Setup1", "Frequency:=", "10GHz", "MaxDeltaS:=", 0.02, "MaximumPasses:=", 15,
    "MinimumPasses:=", 2, "MinimumConvergedPasses:=", 2, "PercentRefinement:=", 30, "IsEnabled:=", True, "BasisOrder:=", 1,
    "DoLambdaRefine:=", True, "DoMaterialLambda:=", True, "SetLambdaTarget:=", False, "Target:=", 0.3333, "UseMaxTetIncrease:=", False])
oAna.InsertFrequencySweep("Setup1", ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearStep", "RangeStart:=", "9.5GHz",
    "RangeEnd:=", "10.5GHz", "RangeStep:=", "0.05GHz", "Type:=", "Interpolating", "SaveFields:=", False, "SaveRadFields:=", False,
    "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 250, "InterpMinSolns:=", 0, "InterpMinSubranges:=", 1])
oRad = oDesign.GetModule("RadField")
oRad.InsertFarFieldSphereSetup(["NAME:3D", "UseCustomRadiationSurface:=", False, "ThetaStart:=", "0deg", "ThetaStop:=", "180deg",
    "ThetaStep:=", "1deg", "PhiStart:=", "-180deg", "PhiStop:=", "180deg", "PhiStep:=", "1deg", "UseLocalCS:=", False])
oEditor.FitAll()
oProject.SaveAs(oDesktop.GetProjectDirectory() + "/phase28_coded_small.aedt", True)
oDesktop.AddMessage(oProject.GetName(), "phase28_coded_small", 0, "Built: %d cells, %d diode sheets (coding beam_030_225), 1 monopole, 1 lumped port. Validate, Analyze; switch beams with set_coding.py." % (len(CELLS), 8 * len(CELLS)))
'''
    )
    (BOARD_DIR / "build_board.py").write_text(build, encoding="utf-8")
    setcoding = (
        head.replace("Builds the geometry ONCE; the diode states of any coding are set by set_coding.py.", "Sets the diode states of the built model from a coding CSV.")
        + "import csv\n"
        + COMMON
        + f"R_ON = \"{dv['r_on_ohm']}ohm\"; C_OFF = \"{dv['c_off_pf']}pF\"; Z = \"{H_MM}mm\"\n"
        + f"NECKS = {_lit(_neck_table())}\n"
        + f"STATE_ON = {_lit(geo.state_on_sites())}\n"
        + f"CELLS = {cells_lit}\n"
        + '''
# EDIT: the full path of the coding CSV (columns cell_id, x_mm, y_mm, state_0to3, axis_deg), e.g. codings/coding_beam_045_225.csv
CODING_CSV = r"codings/coding_beam_030_225.csv"
state = {}
with open(CODING_CSV) as fh:
    for row in csv.DictReader(fh):
        state[int(row["cell_id"])] = int(row["state_0to3"])
assert sorted(state) == [c[0] for c in CELLS], "coding does not match the 102-cell layout"
oProject = oDesktop.GetActiveProject()
oDesign = oProject.GetActiveDesign()
oBnd = oDesign.GetModule("BoundarySetup")
used_delete = 0
for cid, cx, cy in CELLS:
    for k, nk in enumerate(NECKS):
        name = "D%d_%d" % (cid, k)
        line = [(cx + nk["line"][0][0], cy + nk["line"][0][1]), (cx + nk["line"][1][0], cy + nk["line"][1][1])]
        args = rlc_args(name, k in STATE_ON[state[cid]], line, Z)
        try:
            oBnd.EditLumpedRLC(name + "_RLC", args)
        except Exception:
            used_delete += 1
            oBnd.DeleteBoundaries([name + "_RLC"])
            oBnd.AssignLumpedRLC(args)
oDesktop.AddMessage(oProject.GetName(), oDesign.GetName(), 0, "Coding %s applied to %d cells (%d via delete+assign)." % (CODING_CSV, len(CELLS), used_delete))
'''
    )
    (BOARD_DIR / "set_coding.py").write_text(setcoding, encoding="utf-8")
    return preds, len(cells)




def main():
    lib = json.loads(LIB.read_text(encoding="utf-8"))
    ok = lib["convergence_mesh_64_to_96"]
    alphas = [r["alpha_beam_np_per_m"] for r in lib["rows"] if r["kind"] == "normal" and r["mesh"] == 96 and r["f_ghz"] == 10.0]
    alpha_g = float(np.mean(alphas))
    print("alpha_G (mean of 4 states, mesh 96, 10.0 GHz) =", alpha_g, "converged:", {k: v["converged"] for k, v in ok.items()})
    rows = write_cell_package(lib)
    preds, ncell = write_board_package(alpha_g)
    summary = {
        "alpha_g_np_per_m": alpha_g,
        "n_cells": ncell,
        "n_diode_boundaries": 8 * ncell,
        "across_flats_mm": sb.ACROSS_FLATS_MM,
        "circumradius_mm": sb.circumradius_mm(),
        "predictions": {BEAM_NAMES[b]: {k: v for k, v in p.items() if k.startswith("table_")} for b, p in preds.items()},
        "n_reference_rows": len(rows),
    }
    for b, p in summary["predictions"].items():
        for t in p.values():
            t.pop("x_mm", None), t.pop("y_mm", None), t.pop("state_index", None)
    (BOARD_DIR / "mom_reference" / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    import phase28_export_readmes as rd

    (CELL_DIR / "README.md").write_text(rd.cell_readme(lib), encoding="utf-8")
    (BOARD_DIR / "README.md").write_text(rd.board_readme(preds, alpha_g, ncell, lib), encoding="utf-8")
    (LIB.parent / "G_LIBRARY.md").write_text(rd.library_report(lib), encoding="utf-8")
    print("written; cells", ncell, "diodes", 8 * ncell)


if __name__ == "__main__":
    main()
