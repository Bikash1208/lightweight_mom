"""AEDT (IronPython) helper - UNTESTED draft: assigns one parallel-RLC boundary per diode sheet (R = R_A for GA_*, R_B for GB_*, C 28.16 fF) from sheets.csv.
Run in HFSS: Tools > Run Script. Create the project variables R_A and R_B first. Check the first boundary by hand (current-line direction) before trusting the rest."""
import csv
import os

oDesign = oProject.GetActiveDesign()
oMod = oDesign.GetModule("BoundarySetup")
here = os.path.dirname(os.path.abspath(__file__))
for r in csv.DictReader(open(os.path.join(here, "sheets.csv"))):
    var = "R_A" if r["group"] == "GA" else "R_B"
    z = r["z_mm"] + "mm"
    oMod.AssignLumpedRLC(["NAME:" + r["sheet"], "Objects:=", [r["sheet"]],
                          ["NAME:CurrentLine", "Start:=", [r["current_line_start_x"] + "mm", r["current_line_start_y"] + "mm", z], "End:=", [r["current_line_end_x"] + "mm", r["current_line_end_y"] + "mm", z]],
                          "RLC Type:=", "Parallel", "UseResist:=", True, "Resistance:=", var, "UseInduct:=", False, "UseCap:=", True, "Capacitance:=", "28.16fF"])
