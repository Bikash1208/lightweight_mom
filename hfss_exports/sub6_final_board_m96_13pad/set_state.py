# -*- coding: utf-8 -*-
# SUB6 final board: switch all 240 diode boundaries between OFF and ON in the ACTIVE project without rebuilding.  MODEL_ONLY / NOT_HFSS_VALIDATED.
# AEDT 2026 R1 IronPython:  Tools > Run Script.  Edit STATE below.  Not tested inside AEDT by the author.
# The 240 RLC boundaries reference the project variables $Rd, $Ld, $Cd; only these change (the mesh is kept).
#   OFF: R 3876 ohm || C 0.1703 pF (L = 100 uH ~ open)      ON: R 162.4 ohm || L 0.6113 nH (C = 10 aF ~ open)
import ScriptEnv

STATE = "ON"            # "OFF" or "ON"

ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oProject = oDesktop.GetActiveProject()
VALS = {"OFF": {"R": "3876ohm", "C": "0.1703pF", "L": "100000nH"},
        "ON": {"R": "162.4ohm", "L": "0.6113nH", "C": "0.00001pF"}}
vals = VALS[STATE]


def setv(name, value):
    oProject.ChangeProperty(["NAME:AllTabs", ["NAME:ProjectVariableTab", ["NAME:PropServers", "ProjectVariables"],
                            ["NAME:ChangedProps", ["NAME:" + name, "Value:=", value]]]])


setv("$Rd", vals["R"])
setv("$Ld", vals["L"])
setv("$Cd", vals["C"])
setv("$STATE_ON", "1" if STATE == "ON" else "0")
print("STATE = %s  ($Rd %s, $Ld %s, $Cd %s). Re-run the analysis (Analyze) for the new state." % (STATE, vals["R"], vals["L"], vals["C"]))
