# -*- coding: utf-8 -*-
# HFSS (AEDT) IronPython script -- 4-port 1-bit HMSA, square cells, gap @@GAP@@ mm.
# Run inside AEDT: Tools > Run Script... (works in AEDT 2019 R1 and later; IronPython 2.7 or CPython).
# Builds everything natively with exact millimetre values, so no STEP import and no 1.99999 dimensions.
#
#   Substrate  Rogers RT/duroid 5880 (tm), @@W@@ x @@W@@ x @@H@@ mm, z = 0 .. @@H@@ mm
#   Ground     PEC sheet at z = 0 with a coax opening (radius r_out) at each feed
#   Patches    @@N@@ squares @@L@@ x @@L@@ mm at z = @@H@@ mm, united into one sheet, Perfect E
#   Monopoles  PEC pins (radius r_pin, height mono_h from the ground), one per port
#   Ports      lumped ports on the coax annulus at z = 0 (r_pin .. r_out), 50 ohm
#   Air        box padded by pad_xy / pad_top / pad_bot, radiation boundary
#   Setup      10.54 GHz driven modal + optional sweep; far-field sphere 1 deg
#
# Edit the design variables (mono_h, r_pin, ...) in HFSS after the build to match your monopole.
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
oProject.InsertDesign("HFSS", "HMSA_4port_gap@@GAP@@", "DrivenModal", "")
oDesign = oProject.SetActiveDesign("HMSA_4port_gap@@GAP@@")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])

BUILD_MONOPOLES = True      # False: geometry only (patches, substrate, ground), no feeds or ports
ADD_SWEEP = True            # 10.0 - 11.0 GHz interpolating sweep in addition to the 10.54 GHz solve

W = "@@W@@mm"; H = "@@H@@mm"; L = "@@L@@mm"
FEEDS = [@@FEEDS@@]   # (x, y) mm, ports 1..4
PATCHES = [  # lower-left corner (x, y) mm of every @@L@@ mm square
@@PATCHES@@
]

def var(name, value):
    oDesign.ChangeProperty(["NAME:AllTabs", ["NAME:LocalVariableTab",
        ["NAME:PropServers", "LocalVariables"],
        ["NAME:NewProps", ["NAME:" + name, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", value]]]])

for n, v in (("mono_h", "7.1mm"), ("r_pin", "@@RPIN@@mm"), ("r_out", "@@ROUT@@mm"),
             ("pad_xy", "15mm"), ("pad_top", "30mm"), ("pad_bot", "15mm")):
    var(n, v)

def attrs(name, material, solve_inside, color="(255 128 0)", transp=0):
    return ["NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=", color, "Transparency:=", transp,
            "PartCoordinateSystem:=", "Global", "UDMId:=", "", "MaterialValue:=", '"' + material + '"',
            "SurfaceMaterialValue:=", '""', "SolveInside:=", solve_inside, "IsMaterialEditable:=", True,
            "UseMaterialAppearance:=", False, "IsLightweight:=", False]

def mm(v):
    return repr(v) + "mm" if isinstance(v, float) else str(v) + "mm"

def rect(name, x, y, z, w, h, color="(255 128 0)"):
    oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True,
        "XStart:=", x, "YStart:=", y, "ZStart:=", z, "Width:=", w, "Height:=", h, "WhichAxis:=", "Z"],
        attrs(name, "vacuum", True, color))

def circle(name, x, y, z, r):
    oEditor.CreateCircle(["NAME:CircleParameters", "IsCovered:=", True, "XCenter:=", x, "YCenter:=", y,
        "ZCenter:=", z, "Radius:=", r, "WhichAxis:=", "Z", "NumSegments:=", "0"], attrs(name, "vacuum", True))

def unite(names):
    oEditor.Unite(["NAME:Selections", "Selections:=", ",".join(names)], ["NAME:UniteParameters", "KeepOriginals:=", False])

def subtract(blank, tools, keep):
    oEditor.Subtract(["NAME:Selections", "Blank Parts:=", blank, "Tool Parts:=", ",".join(tools)],
                     ["NAME:SubtractParameters", "KeepOriginals:=", keep])

half = "(-" + W + "/2)"
# substrate
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", half, "YPosition:=", half, "ZPosition:=", "0mm",
    "XSize:=", W, "YSize:=", W, "ZSize:=", H], attrs("Substrate", "Rogers RT/duroid 5880 (tm)", True, "(0 128 0)", 0.6))
# ground
rect("Ground", half, half, "0mm", W, W, "(128 128 128)")
# patches, united in blocks of 100 then into one body
names = []
for k, (x, y) in enumerate(PATCHES):
    nm = "P%d" % k
    rect(nm, mm(x), mm(y), H, L, L)
    names.append(nm)
blocks = []
for b in range(0, len(names), 100):
    grp = names[b:b + 100]
    if len(grp) > 1:
        unite(grp)
    blocks.append(grp[0])
if len(blocks) > 1:
    unite(blocks)
oEditor.ChangeProperty(["NAME:AllTabs", ["NAME:Geometry3DAttributeTab", ["NAME:PropServers", names[0]],
    ["NAME:ChangedProps", ["NAME:Name", "Value:=", "Patches"]]]])

oBnd = oDesign.GetModule("BoundarySetup")
if BUILD_MONOPOLES:
    pins = []; holes = []
    for k, (x, y) in enumerate(FEEDS):
        p = k + 1
        oEditor.CreateCylinder(["NAME:CylinderParameters", "XCenter:=", mm(x), "YCenter:=", mm(y), "ZCenter:=", "0mm",
            "Radius:=", "r_pin", "Height:=", "mono_h", "WhichAxis:=", "Z", "NumSides:=", "0"],
            attrs("Monopole%d" % p, "pec", False, "(255 255 0)"))
        pins.append("Monopole%d" % p)
        circle("Hole%d" % p, mm(x), mm(y), "0mm", "r_out"); holes.append("Hole%d" % p)
        circle("PortSheet%d" % p, mm(x), mm(y), "0mm", "r_out")
        circle("PinCut%d" % p, mm(x), mm(y), "0mm", "r_pin")
        subtract("PortSheet%d" % p, ["PinCut%d" % p], False)
    subtract("Ground", holes, False)
    subtract("Substrate", pins, True)
    for k, (x, y) in enumerate(FEEDS):
        p = k + 1
        oBnd.AssignLumpedPort(["NAME:Port%d" % p, "Objects:=", ["PortSheet%d" % p], "DoDeembed:=", False,
            "RenormalizeAllTerminals:=", True,
            ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,
                ["NAME:IntLine", "Start:=", [mm(x) + "+r_pin", mm(y), "0mm"], "End:=", [mm(x) + "+r_out", mm(y), "0mm"]],
                "AlignmentGroup:=", 0, "CharImp:=", "Zpi", "RenormImp:=", "50ohm"]],
            "ShowReporterFilter:=", False, "ReporterFilter:=", [True], "Impedance:=", "50ohm"])

oBnd.AssignPerfectE(["NAME:PerfE_Ground", "Objects:=", ["Ground"], "InfGroundPlane:=", False])
oBnd.AssignPerfectE(["NAME:PerfE_Patches", "Objects:=", ["Patches"], "InfGroundPlane:=", False])

# air box + radiation
oEditor.CreateBox(["NAME:BoxParameters", "XPosition:=", half + "-pad_xy", "YPosition:=", half + "-pad_xy",
    "ZPosition:=", "-pad_bot", "XSize:=", W + "+2*pad_xy", "YSize:=", W + "+2*pad_xy",
    "ZSize:=", "pad_bot+" + H + "+pad_top"], attrs("Air", "air", True, "(143 175 143)", 0.9))
oBnd.AssignRadiation(["NAME:Rad1", "Objects:=", ["Air"], "IsFssReference:=", False, "IsForPML:=", False])

oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven", ["NAME:Setup1", "Frequency:=", "10.54GHz", "MaxDeltaS:=", 0.02,
    "MaximumPasses:=", 15, "MinimumPasses:=", 2, "MinimumConvergedPasses:=", 2, "PercentRefinement:=", 30,
    "IsEnabled:=", True, "BasisOrder:=", 1, "DoLambdaRefine:=", True, "DoMaterialLambda:=", True,
    "SetLambdaTarget:=", False, "Target:=", 0.3333, "UseMaxTetIncrease:=", False])
if ADD_SWEEP:
    oAna.InsertFrequencySweep("Setup1", ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearStep",
        "RangeStart:=", "10GHz", "RangeEnd:=", "11GHz", "RangeStep:=", "0.02GHz", "Type:=", "Interpolating",
        "SaveFields:=", False, "SaveRadFields:=", False, "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 250,
        "InterpMinSolns:=", 0, "InterpMinSubranges:=", 1])
oRad = oDesign.GetModule("RadField")
oRad.InsertFarFieldSphereSetup(["NAME:3D", "UseCustomRadiationSurface:=", False,
    "ThetaStart:=", "0deg", "ThetaStop:=", "90deg", "ThetaStep:=", "0.5deg",
    "PhiStart:=", "0deg", "PhiStop:=", "360deg", "PhiStep:=", "1deg", "UseLocalCS:=", False])
oEditor.FitAll()
oProject.SaveAs(oDesktop.GetProjectDirectory() + "/HMSA_4port_square_gap@@GAP@@.aedt", True)
oDesktop.AddMessage(oProject.GetName(), "HMSA_4port_gap@@GAP@@", 0, "Built: @@N@@ patches, 4 monopoles, 4 lumped ports. Validate, then Analyze All.")
