# -*- coding: utf-8 -*-
# SUB6 final board, RF-ONLY, built NATIVELY in AEDT (no STEP import).  MODEL_ONLY / NOT_HFSS_VALIDATED.
# AEDT 2026 R1, IronPython 2.7:  Tools > Run Script...  Run it in an EMPTY new HFSS Driven Modal design that is ACTIVE (the script uses the active project/design
# if one is open, otherwise it creates project + design).  Data: rf_board_cells.csv next to this script (758 cells; vias of the 240 switched cells).  Units mm.
#
# Geometry (all dimensions from the validated layout, thin stack):
#   top copper sheets at z 0: discs (circle r minus the rotated slot rectangle) METAL_PEC_DISC_<site>, islands METAL_PEC_ISL_<site>, diode sheets D<string>_<pos>
#   (rectangle EXACTLY spanning the 1.00 mm gap, cell-local y 0.45 ... 1.45, 0.35 mm along the slot; edges coincident with the island edge and the slot edge, NO overlap);
#   ground sheet z -3.175 (hexagon minus 480 anti-pad circles r 0.5 and the opening r 2.05) METAL_PEC_GROUND; PORT_ANNULUS r 0.635 ... 2.05 at z -3.175 (separate object);
#   RT5880 SUBSTRATE z -3.175 ... 0, FR4 SUBSTRATE_FR4 z -3.975 ... -3.175 (hexagonal prisms, circumradius 211.547), AIR box minus the two dielectrics;
#   480 vias: PEC solid cylinders r 0.15 SPLIT per layer (VIA_<s>_<p><I|D>_RT z -3.175 ... 0 and VIA_<s>_<p><I|D>_FR4 z -3.975 ... -3.175) + via pad circle r 0.30 at z -3.975 (PAD_...);
#   pin: METAL_PEC_PIN_SUB r 0.635 z -3.175 ... 0 (inside the RT5880) and METAL_PEC_PIN_AIR z 0 ... 8.925 (inside the AIR box).
# No Unite anywhere.  Subtract only: slot from disc (758), 481 circles from the ground, pin circle from the port annulus, dielectrics from AIR.
#
# API calls I am NOT sure about (recorded-script syntax written from memory, never run in AEDT): CreatePolyline / SweepAlongVector (hexagon prisms),
# Rotate / Move argument names, Subtract with a long comma separated tool list, oDesktop.AddMessage, GetFaceIDs / GetFaceCenter, AssignLumpedRLC with the
# boundary named like the object, AssignLumpedPort argument list, InsertFarFieldSphereSetup.  If one is rejected, record the same action in the GUI and compare.
import ScriptEnv
import os
import csv
import math

# ----------------------------------------------------------------- configuration
DIR = r"C:\Users\Bikash\Desktop\sub6_board_hfss"      # folder with rf_board_cells.csv
MAX_CELLS = 0                                         # 0 = all 758 cells; N > 0: only the first N cells (quick test of the script)
STATE = "OFF"                                         # "OFF" or "ON" initial state of the 240 diode boundaries (switch later with set_state.py)
USE_WAVE_PORT = False
DO_SWEEP = False
HEX_R = 211.547                                       # board hexagon circumradius (vertices up)
Z_TOP, Z_GND, Z_BOT = 0.0, -3.175, -3.975
R_BAR, R_PAD, R_ANTI, R_OPEN, R_PIN = 0.15, 0.30, 0.50, 2.05, 0.635
PIN_TOP = 8.925
AIR_PAD_XY, AIR_PAD_BOT, AIR_PAD_TOP = 15.0, 15.0, 40.0
ISL_L, ISL_W, ISL_OY = 4.348, 1.30, -0.20             # island 4.348 x 1.30, centre at cell-local (0, -0.20)  (y -0.85 ... 0.45)
GAP_Y0, GAP_Y1, SH_ALONG = 0.45, 1.45, 0.35           # diode gap between the island edge (y 0.45) and the slot edge (y 1.45)
OFF = {"R": "3876ohm", "C": "0.1703pF", "L": "100000nH"}   # L = 100 uH ~ open
ON = {"R": "162.4ohm", "L": "0.6113nH", "C": "0.00001pF"}  # C = 10 aF ~ open

ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.GetActiveProject()
if oProject is None:
    oProject = oDesktop.NewProject()
oDesign = oProject.GetActiveDesign()
if oDesign is None:
    oProject.InsertDesign("HFSS", "BOARD_N240_RF", "DrivenModal", "")
    oDesign = oProject.SetActiveDesign("BOARD_N240_RF")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])
pname, dname = oProject.GetName(), oDesign.GetName()
COUNT = {"n": 0}
GROUPS = {}


def log(msg):
    print(msg)
    try:
        oDesktop.AddMessage(pname, dname, 0, msg)
    except Exception, e:
        pass


def mm(v):
    return "%.9fmm" % v


def tick(group):
    GROUPS[group] = GROUPS.get(group, 0) + 1
    COUNT["n"] += 1
    if COUNT["n"] % 100 == 0:
        log("created %d objects (%s ...)" % (COUNT["n"], group))


def attrs(name, material="vacuum", solve_inside=True):
    return ["NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=", "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=", "Global", "UDMId:=", "",
            "MaterialValue:=", '"%s"' % material, "SurfaceMaterialValue:=", '""', "SolveInside:=", solve_inside, "ShellElement:=", False, "ShellElementThickness:=", "0mm",
            "IsMaterialEditable:=", True, "UseMaterialAppearance:=", False, "IsLightweight:=", False]


def circle(name, x, y, z, r, group):
    oEditor.CreateCircle(["NAME:CircleParameters", "IsCovered:=", True, "XCenter:=", mm(x), "YCenter:=", mm(y), "ZCenter:=", mm(z), "Radius:=", mm(r), "WhichAxis:=", "Z", "NumSegments:=", "0"], attrs(name))
    tick(group)


def rect_centered(name, w, h, group):
    """rectangle w x h centred on the global origin in the z = 0 plane (to be rotated and moved)"""
    oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", mm(-w / 2.0), "YStart:=", mm(-h / 2.0), "ZStart:=", mm(0.0), "Width:=", mm(w), "Height:=", mm(h), "WhichAxis:=", "Z"], attrs(name))
    tick(group)


def rotate_move(name, ang_deg, x, y, z):
    oEditor.Rotate(["NAME:Selections", "Selections:=", name, "NewPartsModelFlag:=", "Model"], ["NAME:RotateParameters", "RotateAxis:=", "Z", "RotateAngle:=", "%.6fdeg" % ang_deg])
    oEditor.Move(["NAME:Selections", "Selections:=", name, "NewPartsModelFlag:=", "Model"], ["NAME:TranslateParameters", "TranslateVectorX:=", mm(x), "TranslateVectorY:=", mm(y), "TranslateVectorZ:=", mm(z)])


def cylinder(name, x, y, z0, r, h, material, group):
    oEditor.CreateCylinder(["NAME:CylinderParameters", "XCenter:=", mm(x), "YCenter:=", mm(y), "ZCenter:=", mm(z0), "Radius:=", mm(r), "Height:=", mm(h), "WhichAxis:=", "Z", "NumSides:=", "0"],
                           attrs(name, material, False if material == "pec" else True))
    tick(group)


def subtract(blank, tools, keep=False):
    oEditor.Subtract(["NAME:Selections", "Blank Parts:=", blank, "Tool Parts:=", ",".join(tools)], ["NAME:SubtractParameters", "KeepOriginals:=", keep])


def hex_prism(name, z0, h, material, group):
    pts = []
    for i in range(7):
        a = math.radians(90 + 60 * (i % 6))
        pts.append(["NAME:PLPoint", "X:=", mm(HEX_R * math.cos(a)), "Y:=", mm(HEX_R * math.sin(a)), "Z:=", mm(z0)])
    segs = [["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", i, "NoOfPoints:=", 2] for i in range(6)]
    oEditor.CreatePolyline(["NAME:PolylineParameters", "IsPolylineCovered:=", True, "IsPolylineClosed:=", True, ["NAME:PolylinePoints"] + pts, ["NAME:PolylineSegments"] + segs,
                            ["NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=", "0mm",
                             "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner"]], attrs(name, material))
    if h != 0.0:
        oEditor.SweepAlongVector(["NAME:Selections", "Selections:=", name, "NewPartsModelFlag:=", "Model"],
                                 ["NAME:VectorSweepParameters", "DraftAngle:=", "0deg", "DraftType:=", "Round", "CheckFaceFaceIntersection:=", False, "SweepVectorX:=", "0mm", "SweepVectorY:=", "0mm", "SweepVectorZ:=", mm(h)])
    tick(group)


# ----------------------------------------------------------------- materials
oDef = oProject.GetDefinitionManager()
for nm, er, td in (("RT5880", 2.2, 0.0009), ("FR4_assumed", 4.4, 0.02)):
    if not oDef.DoesExist(nm):
        oDef.AddMaterial(["NAME:" + nm, "CoordinateSystemType:=", "Cartesian", "BulkOrSurfaceType:=", 1, ["NAME:PhysicsTypes", "set:=", ["Electromagnetic"]],
                          "permittivity:=", str(er), "dielectric_loss_tangent:=", str(td)])

# ----------------------------------------------------------------- read the cell data
cells = list(csv.DictReader(open(os.path.join(DIR, "rf_board_cells.csv"))))
if MAX_CELLS > 0:
    cells = cells[:MAX_CELLS]
log("cells read: %d (switched %d)" % (len(cells), len([c for c in cells if c["class"] == "SWITCHED"])))

# ----------------------------------------------------------------- dielectrics and air (boxes)
hex_prism("SUBSTRATE", Z_GND, Z_TOP - Z_GND, "RT5880", "dielectric")
hex_prism("SUBSTRATE_FR4", Z_BOT, Z_GND - Z_BOT, "FR4_assumed", "dielectric")
ax = HEX_R * math.cos(math.radians(30)) + AIR_PAD_XY
ay = HEX_R + AIR_PAD_XY
az0, az1 = Z_BOT - AIR_PAD_BOT, 0.035 + AIR_PAD_TOP
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", mm(-ax), "YPosition:=", mm(-ay), "ZPosition:=", mm(az0), "XSize:=", mm(2 * ax), "YSize:=", mm(2 * ay), "ZSize:=", mm(az1 - az0)], attrs("AIR", "vacuum"))
tick("air")
subtract("AIR", ["SUBSTRATE", "SUBSTRATE_FR4"], True)      # AIR minus both dielectrics (clones: the dielectrics stay)

# ----------------------------------------------------------------- top copper: discs, islands, diode sheets
sheet_rows = []      # (name, string, pos, x, y, axis) for the RLC boundaries
via_rows = []        # (tag, string, pos, x, y)
for c in cells:
    site = c["site"]
    x, y, R = float(c["x"]), float(c["y"]), float(c["r"])
    L_, W_, ang = float(c["slot_L"]), float(c["slot_W"]), float(c["axis_deg"])
    dn = "METAL_PEC_DISC_" + site
    circle(dn, x, y, Z_TOP, R, "disc")
    sn = "SLOT_" + site
    rect_centered(sn, L_, W_, "tmp")
    rotate_move(sn, ang, x, y, Z_TOP)
    subtract(dn, [sn], False)
    if c["class"] != "SWITCHED":
        continue
    s, p = int(c["string"]), int(c["pos"])
    th = math.radians(ang)
    cs, sn_ = math.cos(th), math.sin(th)
    # island: centre at cell-local (0, ISL_OY)
    ix, iy = x - ISL_OY * 0 + (-sn_) * ISL_OY, y + cs * ISL_OY
    iname = "METAL_PEC_ISL_" + site
    rect_centered(iname, ISL_L, ISL_W, "island")
    rotate_move(iname, ang, ix, iy, Z_TOP)
    # diode sheet: centre at cell-local (0, 0.95), 0.35 along the slot, 1.00 across (exactly the gap, no overlap)
    gy = 0.5 * (GAP_Y0 + GAP_Y1)
    gx_, gy_ = x + (-sn_) * gy, y + cs * gy
    shn = "D%d_%d" % (s, p)
    rect_centered(shn, SH_ALONG, GAP_Y1 - GAP_Y0, "diode_sheet")
    rotate_move(shn, ang, gx_, gy_, Z_TOP)
    a_pt = (x + (-sn_) * GAP_Y0, y + cs * GAP_Y0)     # anode edge (island side)
    k_pt = (x + (-sn_) * GAP_Y1, y + cs * GAP_Y1)     # cathode edge (disc side)
    sheet_rows.append((shn, s, p, a_pt, k_pt, ang))
    via_rows.append(("D", s, p, float(c["dsc_via_x"]), float(c["dsc_via_y"])))
    via_rows.append(("I", s, p, float(c["isl_via_x"]), float(c["isl_via_y"])))
log("top copper done: %s" % str(GROUPS))

# ----------------------------------------------------------------- vias (split per dielectric layer) and via pads
for tag, s, p, vx, vy in via_rows:
    base = "VIA_%d_%d%s" % (s, p, tag)
    cylinder(base + "_RT", vx, vy, Z_GND, R_BAR, Z_TOP - Z_GND, "pec", "via_solid")
    cylinder(base + "_FR4", vx, vy, Z_BOT, R_BAR, Z_GND - Z_BOT, "pec", "via_solid")
    circle("PAD_%d_%d%s" % (s, p, tag), vx, vy, Z_BOT, R_PAD, "via_pad")
log("vias done: %s" % str(GROUPS))

# ----------------------------------------------------------------- ground sheet (hexagon minus anti-pads and the opening), port annulus, pin
hex_prism("METAL_PEC_GROUND", Z_GND, 0.0, "vacuum", "ground")
tools = []
for k, (tag, s, p, vx, vy) in enumerate(via_rows):
    nm = "AP_%d" % k
    circle(nm, vx, vy, Z_GND, R_ANTI, "tmp")
    tools.append(nm)
circle("AP_OPEN", 0.0, 0.0, Z_GND, R_OPEN, "tmp")
tools.append("AP_OPEN")
subtract("METAL_PEC_GROUND", tools, False)                 # the opening circle is consumed: the port annulus is a NEW object
circle("PORT_ANNULUS", 0.0, 0.0, Z_GND, R_OPEN, "port")
circle("PIN_HOLE_TMP", 0.0, 0.0, Z_GND, R_PIN, "tmp")
subtract("PORT_ANNULUS", ["PIN_HOLE_TMP"], False)
cylinder("METAL_PEC_PIN_SUB", 0.0, 0.0, Z_GND, R_PIN, Z_TOP - Z_GND, "pec", "pin")
cylinder("METAL_PEC_PIN_AIR", 0.0, 0.0, Z_TOP, R_PIN, PIN_TOP - Z_TOP, "pec", "pin")
log("ground / port / pin done")

# ----------------------------------------------------------------- object counts per group
solids = list(oEditor.GetObjectsInGroup("Solids"))
sheets = list(oEditor.GetObjectsInGroup("Sheets"))


def cnt(lst, pre):
    return len([n for n in lst if n.startswith(pre)])


log("COUNTS: solids %d, sheets %d | DISC %d, ISL %d, D sheets %d, PAD %d, VIA_RT %d, VIA_FR4 %d, GROUND %d, PORT %d, PIN %d, dielectrics+air %d" % (
    len(solids), len(sheets), cnt(sheets, "METAL_PEC_DISC_"), cnt(sheets, "METAL_PEC_ISL_"), len([n for n in sheets if n.startswith("D") and "_" in n and n[1:2].isdigit()]),
    cnt(sheets, "PAD_"), len([n for n in solids if n.endswith("_RT") and n.startswith("VIA_")]), len([n for n in solids if n.endswith("_FR4") and n.startswith("VIA_")]),
    cnt(sheets, "METAL_PEC_GROUND"), cnt(sheets, "PORT_ANNULUS"), cnt(solids, "METAL_PEC_PIN"), len([n for n in solids if n in ("SUBSTRATE", "SUBSTRATE_FR4", "AIR")])))

# ----------------------------------------------------------------- boundaries
oBnd = oDesign.GetModule("BoundarySetup")
pec_sheets = [n for n in sheets if n.startswith("METAL_PEC_") or n.startswith("PAD_")]
oBnd.AssignPerfectE(["NAME:PEC_sheets", "Objects:=", pec_sheets, "InfGroundPlane:=", False])


def set_project_variable(name, value):
    oProject.ChangeProperty(["NAME:AllTabs", ["NAME:ProjectVariableTab", ["NAME:PropServers", "ProjectVariables"],
                            ["NAME:NewProps", ["NAME:" + name, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", value]]]])


st = ON if STATE == "ON" else OFF
set_project_variable("$Rd", st["R"])
set_project_variable("$Ld", st["L"])
set_project_variable("$Cd", st["C"])
set_project_variable("$STATE_ON", "1" if STATE == "ON" else "0")
# 240 lumped RLC boundaries, one per diode sheet, by NAME; the current line runs across the gap from the anode edge (island) to the cathode edge (disc)
for k, (shn, s, p, a_pt, k_pt, ang) in enumerate(sheet_rows):
    if k == 0:
        dx, dy = k_pt[0] - a_pt[0], k_pt[1] - a_pt[1]
        dot = dx * math.cos(math.radians(ang)) + dy * math.sin(math.radians(ang))
        log("CHECK sheet %s: line length %.6f mm (1.000), dot(line, slot axis) = %.2e (0)" % (shn, math.hypot(dx, dy), dot))
        if abs(math.hypot(dx, dy) - 1.0) > 1e-6 or abs(dot) > 1e-6:
            raise Exception("diode current line is not across the slot")
    params = ["NAME:" + shn, "Objects:=", [shn], ["NAME:CurrentLine", "Start:=", [mm(a_pt[0]), mm(a_pt[1]), mm(Z_TOP)], "End:=", [mm(k_pt[0]), mm(k_pt[1]), mm(Z_TOP)]],
              "UseResist:=", True, "Resistance:=", "$Rd", "UseInduct:=", True, "Inductance:=", "$Ld", "UseCap:=", True, "Capacitance:=", "$Cd"]
    try:
        oBnd.AssignLumpedRLC(params)
    except Exception, e:
        params[0] = "NAME:RLC_" + shn
        oBnd.AssignLumpedRLC(params)
log("assigned %d lumped RLC boundaries, state %s" % (len(sheet_rows), STATE))
# port on PORT_ANNULUS (in the ground plane): integration line from the pin (r 0.635) to the opening edge (r 2.05) along +x at z -3.175
line = ["NAME:IntLine", "Start:=", [mm(R_PIN), "0mm", mm(Z_GND)], "End:=", [mm(R_OPEN), "0mm", mm(Z_GND)]]
if USE_WAVE_PORT:
    oBnd.AssignWavePort(["NAME:1", "Objects:=", ["PORT_ANNULUS"], "NumModes:=", 1, "RenormalizeAllTerminals:=", True, "UseLineModeAlignment:=", False, "DoDeembed:=", False,
                         ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True, line, "AlignmentGroup:=", 0, "CharImp:=", "Zpi", "RenormImp:=", "50ohm"]], "ShowReporterFilter:=", False, "ReporterFilterRealOnly:=", False])
else:
    oBnd.AssignLumpedPort(["NAME:1", "Objects:=", ["PORT_ANNULUS"], "DoDeembed:=", False,
                           ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True, line, "CharImp:=", "Zpi", "AlignmentGroup:=", 0, "RenormImp:=", "50ohm"]], "ShowReporterFilter:=", False, "ReporterFilterRealOnly:=", False])
# radiation boundary on the 6 outer faces of AIR (face centres on the box limits)
ids = []
for fid in oEditor.GetFaceIDs("AIR"):
    ctr = oEditor.GetFaceCenter(fid)
    cx, cy, cz = float(ctr[0]), float(ctr[1]), float(ctr[2])
    if abs(abs(cx) - ax) < 1e-3 or abs(abs(cy) - ay) < 1e-3 or abs(cz - az1) < 1e-3 or abs(cz - az0) < 1e-3:
        ids.append(int(fid))
log("radiation faces: %d (expect 6)" % len(ids))
oBnd.AssignRadiation(["NAME:Rad1", "Faces:=", ids, "IsFssReference:=", False, "IsForPML:=", False])
# analysis setup and far field
oAn = oDesign.GetModule("AnalysisSetup")
oAn.InsertSetup("HfssDriven", ["NAME:Setup1", "SolveType:=", "Single", "Frequency:=", "5.8GHz", "MaxDeltaS:=", 0.02, "MaximumPasses:=", 10, "MinimumPasses:=", 2,
                               "MinimumConvergedPasses:=", 1, "PercentRefinement:=", 30, "IsEnabled:=", True, "BasisOrder:=", 1, "UseIterativeSolver:=", False])
if DO_SWEEP:
    oAn.InsertFrequencySweep("Setup1", ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearCount", "RangeStart:=", "5.6GHz", "RangeEnd:=", "6GHz", "RangeCount:=", 3,
                                        "Type:=", "Discrete", "SaveFields:=", True, "SaveRadFields:=", True, "ExtrapToDC:=", False])
oRad = oDesign.GetModule("RadField")
oRad.InsertFarFieldSphereSetup(["NAME:Infinite Sphere1", "UseCustomRadiationSurface:=", False, "CSDefinition:=", "Theta-Phi", "Polarization:=", "Linear",
                                "ThetaStart:=", "0deg", "ThetaStop:=", "90deg", "ThetaStep:=", "1deg", "PhiStart:=", "0deg", "PhiStop:=", "360deg", "PhiStep:=", "1deg", "UseLocalCS:=", False])
try:
    oProject.SaveAs(os.path.join(DIR, "BOARD_N240_RF_native.aedt"), True)
except Exception, e:
    log("SaveAs skipped: %s" % str(e))
log("DONE: native RF board built, state %s. Run Setup1, switch states with set_state.py." % STATE)
