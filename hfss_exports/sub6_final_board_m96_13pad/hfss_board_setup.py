# -*- coding: utf-8 -*-
# SUB6 final board (13-pad bias layout, SMA flat on the FR4 side), HFSS setup script.  MODEL_ONLY / NOT_HFSS_VALIDATED.
# AEDT 2026 R1, IronPython 2.7:  Tools > Run Script...   (no PyAEDT).
# NOT TESTED inside AEDT by the author (no AEDT licence in the build session): the API calls follow the recorded-script syntax of AEDT;
# if one call is rejected, record the same action in the GUI and compare.
#
# What it does: new Driven Modal design, imports the STEP (mm), assigns materials, PEC, 240 lumped-RLC diode boundaries D<string>_<pos>,
# lumped port (wave port optional), radiation boundary, setup 5.8 GHz, optional sweep, infinite sphere, variable-driven state OFF/ON.
import ScriptEnv
import os
import csv
import math

# ----------------------------------------------------------------- configuration (edit here)
DIR = r"C:\Users\Bikash\Desktop\sub6_board_hfss"          # folder with the STEP and diode_sheets.csv
THIN = True                                               # True: BOARD_N240_thin_clean.step (zero-thickness copper, recommended); False: BOARD_N240_clean.step (35 um copper solids)
STATE = "OFF"                                             # "OFF" or "ON" (initial state; switch later with set_state.py)
USE_WAVE_PORT = False                                     # False: lumped port 50 ohm on PORT_ANNULUS; True: wave port on the same sheet
DO_SWEEP = False                                          # discrete sweep 5.6 / 5.8 / 6.0 GHz
UNITE_METAL = True                                        # unite all METAL_PEC_* solids into one object METAL_PEC (falls back to separate objects if HFSS refuses)
DEFENSIVE_BOOLEANS = True                                 # only used if the validation check reports intersections: subtract (clones) metal from dielectrics, everything from AIR
AIR_BOX = {"ax": 198.205, "ay": 226.547}                  # AIR half sizes x / y (mm), from geometry_build_*.json -> meta.air_box
# thick STEP values (THIN overrides them below)
Z_AIR_LO, Z_AIR_HI = -19.045, 40.035
Z_SHEET = 0.0175                                          # diode sheet z
Z_PORT = -4.010                                           # PORT_ANNULUS z (bottom of the FR4)
R_PIN, R_OUT = 0.635, 2.05
# diode models at 5.8 GHz (parallel equivalents of the BAR50 model)
OFF = {"R": "3876ohm", "C": "0.1703pF", "L": "100000nH"}  # L = 100 uH ~ open
ON = {"R": "162.4ohm", "L": "0.6113nH", "C": "0.00001pF"}  # C = 10 aF ~ open
# ----------------------------------------------------------------- start AEDT
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
if THIN:
    Z_SHEET, Z_PORT, Z_AIR_LO = 0.0, -3.975, -18.975
step_path = os.path.join(DIR, "BOARD_N240_thin_clean.step" if THIN else "BOARD_N240_clean.step")
csv_path = os.path.join(DIR, "diode_sheets.csv")


def q(v):
    return "%.6fmm" % v


oProject = oDesktop.NewProject()
oProject.InsertDesign("HFSS", "BOARD_N240", "DrivenModal", "")
oDesign = oProject.SetActiveDesign("BOARD_N240")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])
pname = oProject.GetName()
dname = oDesign.GetName()

# ----------------------------------------------------------------- materials
oDef = oProject.GetDefinitionManager()


def add_material(name, er, tand):
    if oDef.DoesExist(name):
        return
    oDef.AddMaterial(["NAME:" + name, "CoordinateSystemType:=", "Cartesian", "BulkOrSurfaceType:=", 1,
                      ["NAME:PhysicsTypes", "set:=", ["Electromagnetic"]],
                      "permittivity:=", str(er), "dielectric_loss_tangent:=", str(tand)])


add_material("RT5880", 2.2, 0.0009)
add_material("FR4_assumed", 4.4, 0.02)

# ----------------------------------------------------------------- import STEP
oEditor.Import(["NAME:NativeBodyParameters", "HealOption:=", 0, "Options:=", "-1", "FileType:=", "UnclassifiedSystem", "MaxStitchTol:=", -1,
                "ImportFreeSurfaces:=", True, "GroupByAssembly:=", False, "CreateGroup:=", False, "STLFileUnit:=", "Auto", "MergeFacesAngle:=", -1,
                "HealingOption:=", 0, "ReduceSurfaces:=", False, "SourceFile:=", step_path])
solids = list(oEditor.GetObjectsInGroup("Solids"))
sheets = list(oEditor.GetObjectsInGroup("Sheets"))
print("imported solids %d, sheets %d" % (len(solids), len(sheets)))
metal_solids = [n for n in solids if n.startswith("METAL_PEC")]
metal_sheets = [n for n in sheets if n.startswith("METAL_PEC")]


def set_material(objs, mat, solve_inside):
    oEditor.ChangeProperty(["NAME:AllTabs", ["NAME:Geometry3DAttributeTab", ["NAME:PropServers"] + list(objs),
                           ["NAME:ChangedProps", ["NAME:Material", "Value:=", '"%s"' % mat], ["NAME:Solve Inside", "Value:=", solve_inside]]]])


set_material(["SUBSTRATE"], "RT5880", True)
set_material(["SUBSTRATE_FR4"], "FR4_assumed", True)
set_material(["AIR"], "vacuum", True)


# ----------------------------------------------------------------- validation (number of intersect errors)
def count_intersections():
    n_int = -1
    try:
        oDesktop.ClearMessages(pname, dname, 2)
        oDesign.ValidateDesign()          # "Validation Check"; if this call is rejected run HFSS > Validation Check manually
    except Exception, e:
        print("ValidateDesign call failed (%s): run Validation Check from the menu and read the message manager" % str(e))
    try:
        msgs = list(oDesktop.GetMessages(pname, dname, 2)) + list(oDesktop.GetMessages(pname, dname, 1))
        n_int = len([m for m in msgs if "intersect" in m.lower()])
    except Exception, e:
        print("GetMessages failed: %s" % str(e))
    return n_int


n_int = count_intersections()
print("INTERSECT ERRORS after import: %d" % n_int)
if n_int != 0 and DEFENSIVE_BOOLEANS:
    # clones of the tool objects (KeepOriginals True) so that the metal and the dielectrics survive
    try:
        oEditor.Subtract(["NAME:Selections", "Blank Parts:=", "SUBSTRATE", "Tool Parts:=", ",".join(metal_solids)], ["NAME:SubtractParameters", "KeepOriginals:=", True])
        oEditor.Subtract(["NAME:Selections", "Blank Parts:=", "SUBSTRATE_FR4", "Tool Parts:=", ",".join(metal_solids)], ["NAME:SubtractParameters", "KeepOriginals:=", True])
        oEditor.Subtract(["NAME:Selections", "Blank Parts:=", "AIR", "Tool Parts:=", ",".join(["SUBSTRATE", "SUBSTRATE_FR4"] + metal_solids)], ["NAME:SubtractParameters", "KeepOriginals:=", True])
    except Exception, e:
        print("defensive subtract failed: %s" % str(e))
    n_int = count_intersections()
    print("INTERSECT ERRORS after defensive booleans: %d" % n_int)

# ----------------------------------------------------------------- metal: unite (solids) and Perfect E
oBnd = oDesign.GetModule("BoundarySetup")
metal_names = list(metal_solids)
if metal_solids:
    if UNITE_METAL and len(metal_solids) > 1:
        try:
            oEditor.Unite(["NAME:Selections", "Selections:=", ",".join(metal_solids)], ["NAME:UniteParameters", "KeepOriginals:=", False])
            metal_names = [n for n in oEditor.GetObjectsInGroup("Solids") if n.startswith("METAL_PEC")]
            print("united METAL_PEC solids -> %s" % str(metal_names[:3]))
        except Exception, e:
            print("unite refused (%s): keeping separate METAL_PEC_* objects" % str(e))
    set_material(metal_names, "pec", False)
if metal_sheets:
    oBnd.AssignPerfectE(["NAME:PEC_sheets", "Objects:=", metal_sheets, "InfGroundPlane:=", False])


# ----------------------------------------------------------------- state variables (project variables drive all 240 RLC boundaries)
def set_project_variable(name, value):
    oProject.ChangeProperty(["NAME:AllTabs", ["NAME:ProjectVariableTab", ["NAME:PropServers", "ProjectVariables"],
                            ["NAME:NewProps", ["NAME:" + name, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", value]]]])


st = ON if STATE == "ON" else OFF
set_project_variable("$Rd", st["R"])
set_project_variable("$Ld", st["L"])
set_project_variable("$Cd", st["C"])
set_project_variable("$STATE_ON", "1" if STATE == "ON" else "0")

# ----------------------------------------------------------------- 240 diode sheets: lumped RLC boundaries D<string>_<pos>
# Current line = across the 1.00 mm diode gap, i.e. across the NARROW dimension of the 8.696 x 2.9 mm slot (cell-local +y, anode edge on the island -> cathode edge on the disc).
# This is the board-frame direction (current_dir_x, current_dir_y) in diode_sheets.csv (it is NOT the unit-cell variable $ph_x; no unit-cell variables are used).
rows = list(csv.DictReader(open(csv_path)))
for k, r in enumerate(rows):
    nm = r["name"]
    ax_, ay_, kx_, ky_ = float(r["anode_edge_x"]), float(r["anode_edge_y"]), float(r["cathode_edge_x"]), float(r["cathode_edge_y"])
    dx, dy = float(r["current_dir_x"]), float(r["current_dir_y"])
    th = math.radians(float(r["slot_axis_deg"]))
    if k == 0:   # check one sheet: the current line must be perpendicular to the slot axis and 1.00 mm long
        dot = dx * math.cos(th) + dy * math.sin(th)
        ln = math.hypot(kx_ - ax_, ky_ - ay_)
        print("CHECK sheet %s: current dir (%.4f, %.4f), slot axis %.3f deg, dot(dir, slot axis) = %.2e (must be 0), line length %.4f mm (must be 1.000)" % (nm, dx, dy, float(r["slot_axis_deg"]), dot, ln))
        if abs(dot) > 1e-3 or abs(ln - 1.0) > 1e-3:
            raise Exception("diode current line is not across the slot: stop")
    params = ["NAME:" + nm, "Objects:=", [nm],
              ["NAME:CurrentLine", "Start:=", [q(ax_), q(ay_), q(Z_SHEET)], "End:=", [q(kx_), q(ky_), q(Z_SHEET)]],
              "UseResist:=", True, "Resistance:=", "$Rd", "UseInduct:=", True, "Inductance:=", "$Ld", "UseCap:=", True, "Capacitance:=", "$Cd"]
    try:
        oBnd.AssignLumpedRLC(params)
    except Exception, e:
        params[0] = "NAME:RLC_" + nm      # boundary name clash with the object name: prefix
        oBnd.AssignLumpedRLC(params)
print("assigned %d lumped RLC boundaries, state %s" % (len(rows), STATE))

# ----------------------------------------------------------------- port on PORT_ANNULUS (integration line pin -> opening edge, radial along +x)
line = ["NAME:IntLine", "Start:=", [q(R_PIN), "0mm", q(Z_PORT)], "End:=", [q(R_OUT), "0mm", q(Z_PORT)]]
if USE_WAVE_PORT:
    oBnd.AssignWavePort(["NAME:1", "Objects:=", ["PORT_ANNULUS"], "NumModes:=", 1, "RenormalizeAllTerminals:=", True, "UseLineModeAlignment:=", False, "DoDeembed:=", False,
                         ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True, line, "AlignmentGroup:=", 0, "CharImp:=", "Zpi", "RenormImp:=", "50ohm"]],
                         "ShowReporterFilter:=", False, "ReporterFilterRealOnly:=", False])
else:
    oBnd.AssignLumpedPort(["NAME:1", "Objects:=", ["PORT_ANNULUS"], "DoDeembed:=", False,
                           ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True, line, "CharImp:=", "Zpi", "AlignmentGroup:=", 0, "RenormImp:=", "50ohm"]],
                           "ShowReporterFilter:=", False, "ReporterFilterRealOnly:=", False])

# ----------------------------------------------------------------- radiation boundary on all outer faces of AIR
ids = []
for fid in oEditor.GetFaceIDs("AIR"):
    c = oEditor.GetFaceCenter(fid)
    cx, cy, cz = float(c[0]), float(c[1]), float(c[2])
    if abs(abs(cx) - AIR_BOX["ax"]) < 1e-3 or abs(abs(cy) - AIR_BOX["ay"]) < 1e-3 or abs(cz - Z_AIR_HI) < 1e-3 or abs(cz - Z_AIR_LO) < 1e-3:
        ids.append(int(fid))
print("radiation faces: %d (expect 6)" % len(ids))
oBnd.AssignRadiation(["NAME:Rad1", "Faces:=", ids, "IsFssReference:=", False, "IsForPML:=", False])

# ----------------------------------------------------------------- analysis setup
oAn = oDesign.GetModule("AnalysisSetup")
oAn.InsertSetup("HfssDriven", ["NAME:Setup1", "SolveType:=", "Single", "Frequency:=", "5.8GHz", "MaxDeltaS:=", 0.02, "MaximumPasses:=", 10, "MinimumPasses:=", 2,
                               "MinimumConvergedPasses:=", 1, "PercentRefinement:=", 30, "IsEnabled:=", True, "BasisOrder:=", 1, "UseIterativeSolver:=", False])
if DO_SWEEP:
    oAn.InsertFrequencySweep("Setup1", ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearCount", "RangeStart:=", "5.6GHz", "RangeEnd:=", "6GHz", "RangeCount:=", 3,
                                        "Type:=", "Discrete", "SaveFields:=", True, "SaveRadFields:=", True, "ExtrapToDC:=", False])
oRad = oDesign.GetModule("RadField")
oRad.InsertFarFieldSphereSetup(["NAME:Infinite Sphere1", "UseCustomRadiationSurface:=", False, "CSDefinition:=", "Theta-Phi", "Polarization:=", "Linear",
                                "ThetaStart:=", "0deg", "ThetaStop:=", "90deg", "ThetaStep:=", "1deg", "PhiStart:=", "0deg", "PhiStop:=", "360deg", "PhiStep:=", "1deg", "UseLocalCS:=", False])
oProject.SaveAs(os.path.join(DIR, "BOARD_N240_%s.aedt" % ("thin" if THIN else "thick")), True)
print("DONE: setup built, state %s, intersect errors %d. Run Setup1, then set_state.py to switch OFF/ON." % (STATE, n_int))
